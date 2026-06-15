"""
filter.py — Peptide/receptor filtering and fixed-CIF generation.

Main entry point: process_entry(pdb_id, n_mer, pdb_cache, fixed_dir)
Returns a result dict with accepted bool and all metadata fields.
"""

import shutil
from difflib import SequenceMatcher
from pathlib import Path

from Bio.PDB.MMCIF2Dict import MMCIF2Dict
import numpy as np
from scipy.spatial import cKDTree
from pdbfixer import PDBFixer
from openmm.app import PDBxFile, Modeller
import openmm.unit as unit

from .metadata import (
    extract_sequence,
    extract_full_sequence,
    extract_resolution,
    extract_ss_string,
    extract_ss_dssp_chain,
    extract_receptor_description,
    classify_receptor_function,
)

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

AMBER14_AA = {
    "ALA", "ARG", "ASN", "ASP", "CYS",
    "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO",
    "SER", "THR", "TRP", "TYR", "VAL",
    "HID", "HIE", "HIP",
}

KNOWN_PTM = {
    "SEP": "phosphoserine",
    "TPO": "phosphothreonine",
    "PTR": "phosphotyrosine",
    "NEP": "N1-phosphohistidine",
    "MLY": "N-methyllysine",
    "M3L": "N-trimethyllysine",
    "ALY": "N6-acetyllysine",
    "HYP": "hydroxyproline",
    "CGU": "gamma-carboxyglutamate",
    "CME": "S-(2-hydroxyethyl)thiocysteine",
    "CSO": "S-hydroxycysteine",
    "OCS": "cysteinesulfinic acid",
    "CSD": "cysteinesulfinic acid",
    "CSS": "S-mercaptocysteine",
    "MSE": "selenomethionine (crystallographic Met substitute)",
    "SEC": "selenocysteine",
    "LYZ": "5-hydroxylysine",
    "DAL": "D-alanine",
    "DCY": "D-cysteine",
    "DGL": "D-glutamic acid",
    "DPN": "D-phenylalanine",
    "DSN": "D-serine",
    "DTH": "D-threonine",
    "DVA": "D-valine",
}

AMBER14_HETATM = {
    "NA", "K", "CL", "CA", "MG", "ZN",
    "FE", "FE2", "CU", "MN",
}

WATER_RESIDUES = {"HOH", "WAT", "H2O", "TIP", "TIP3", "TIP3P", "SOL"}

# MSE → MET: crystallographic phasing artefact, biologically identical.
# PDBFixer converts these via replaceNonstandardResidues() during build.
PDBFIXER_CONVERTS = {"MSE"}

# Terminal capping groups: not real residues, can be stripped if not at interface.
TERMINAL_CAPS = {"ACE", "NME", "NH2"}

CAP_INTERFACE_CUTOFF_A = 4.5  # caps within this distance of receptor are kept / rejected

NEIGHBOR_CUTOFF_A    = 4.0
MAX_RECEPTOR_GAP     = 10
PEPTIDE_BOND_CUTOFF  = 1.8
PEPTIDE_COPY_THRESHOLD = 0.7
CONTACT_CUTOFF_A     = 4.5   # peptide heavy-atom contacts
HBOND_CUTOFF_A       = 3.5   # N/O/S donor-acceptor proximity proxy
HBOND_ELEMENTS       = {"N", "O", "S"}


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def is_protein_chain(residues):
    return any(r in AMBER14_AA or r in KNOWN_PTM for r in residues)


def peptide_similarity(residues_a, residues_b):
    return SequenceMatcher(None, residues_a, residues_b).ratio()


def check_amber14(residue_names, exempt=frozenset()):
    """Return list of non-AMBER14 residue names (empty = pass).
    exempt: additional residue names to skip (e.g. PDBFIXER_CONVERTS, TERMINAL_CAPS)."""
    bad = []
    for name in residue_names:
        if name in WATER_RESIDUES or name in exempt:
            continue
        if name not in AMBER14_AA and name not in AMBER14_HETATM:
            bad.append(name)
    return bad


def ptm_description(resname):
    return KNOWN_PTM.get(resname, f"unknown non-standard residue {resname}")


def raw_chain_data(cif_path):
    """Parse CIF atom_site and entity tables.
    Returns (observed, full_seq_auth, full_seq_label).
    """
    d = MMCIF2Dict(str(cif_path))
    auth_ids   = d.get("_atom_site.auth_asym_id",    [])
    label_ids  = d.get("_atom_site.label_asym_id",   [])
    res_names  = d.get("_atom_site.auth_comp_id",    [])
    res_nums   = d.get("_atom_site.auth_seq_id",     [])
    entity_ids = d.get("_atom_site.label_entity_id", [])

    observed = {}
    auth_to_entity  = {}
    label_to_entity = {}
    for auth, label, name, num, eid in zip(auth_ids, label_ids, res_names, res_nums, entity_ids):
        if name.upper() not in WATER_RESIDUES:
            observed.setdefault(auth, set()).add(num)
        if auth not in auth_to_entity:
            auth_to_entity[auth] = eid
        if label not in label_to_entity:
            label_to_entity[label] = eid
    observed = {c: len(s) for c, s in observed.items()}

    entity_seq_len = {}
    for eid, seq in zip(
        d.get("_entity_poly.entity_id",                []),
        d.get("_entity_poly.pdbx_seq_one_letter_code", []),
    ):
        entity_seq_len[eid] = sum(1 for ch in seq if ch.isalpha())

    full_seq_auth  = {a: entity_seq_len[e] for a, e in auth_to_entity.items()  if e in entity_seq_len}
    full_seq_label = {l: entity_seq_len[e] for l, e in label_to_entity.items() if e in entity_seq_len}
    return observed, full_seq_auth, full_seq_label


# ---------------------------------------------------------------------------
# Interface contacts and H-bond proxies
# ---------------------------------------------------------------------------

def compute_interface_contacts(fixed_cif_path, pep_chain, rec_chains):
    """Count peptide-receptor contacts and H-bond proxies from a fixed CIF.

    Returns (n_contacts, n_hbonds):
      n_contacts — peptide heavy atoms within CONTACT_CUTOFF_A of any receptor heavy atom.
      n_hbonds   — N/O/S atoms from peptide within HBOND_CUTOFF_A of N/O/S atoms from
                   receptor (distance-only proxy; no explicit H atoms required).
    """
    try:
        d = MMCIF2Dict(str(fixed_cif_path))
    except Exception:
        return 0, 0

    chains   = d.get("_atom_site.auth_asym_id", [])
    xs       = d.get("_atom_site.Cartn_x", [])
    ys       = d.get("_atom_site.Cartn_y", [])
    zs       = d.get("_atom_site.Cartn_z", [])
    elements = d.get("_atom_site.type_symbol", [])

    rec_set = set(rec_chains)
    pep_coords, pep_elems = [], []
    rec_coords, rec_elems = [], []

    for c, x, y, z, el in zip(chains, xs, ys, zs, elements):
        try:
            pt = (float(x), float(y), float(z))
        except (ValueError, TypeError):
            continue
        el_up = el.strip().upper()
        if c == pep_chain:
            pep_coords.append(pt)
            pep_elems.append(el_up)
        elif c in rec_set:
            rec_coords.append(pt)
            rec_elems.append(el_up)

    if not pep_coords or not rec_coords:
        return 0, 0

    pep_arr = np.array(pep_coords, dtype=np.float32)
    rec_arr = np.array(rec_coords, dtype=np.float32)

    dists, _ = cKDTree(rec_arr).query(pep_arr, k=1, workers=1)
    n_contacts = int((dists < CONTACT_CUTOFF_A).sum())

    pep_hb_idx = [i for i, e in enumerate(pep_elems) if e in HBOND_ELEMENTS]
    rec_hb_idx = [i for i, e in enumerate(rec_elems) if e in HBOND_ELEMENTS]
    if not pep_hb_idx or not rec_hb_idx:
        return n_contacts, 0

    dists_hb, _ = cKDTree(rec_arr[rec_hb_idx]).query(pep_arr[pep_hb_idx], k=1, workers=1)
    n_hbonds = int((dists_hb < HBOND_CUTOFF_A).sum())

    return n_contacts, n_hbonds


# ---------------------------------------------------------------------------
# Fixed CIF builder
# ---------------------------------------------------------------------------

def build_fixed_cif(pdb_id, cif_path, peptide_chain_id, receptor_chain_ids,
                    fixed_dir, ion_chain_ids=None, convert_mse=True, caps_to_strip=None):
    """Load raw CIF, strip non-system chains, fill missing residues, write fixed CIF.
    Returns (True, {residues_added, chain_id_map}) or (False, error_msg).
    convert_mse: call replaceNonstandardResidues() to convert MSE→MET.
    caps_to_strip: list of cap residue names (ACE/NME/NH2) to remove from peptide terminus.
    """
    fixed_dir = Path(fixed_dir)
    try:
        fixer = PDBFixer(filename=str(cif_path))
    except Exception as exc:
        return False, f"PDBFixer reload error: {exc}"

    keep_ids = {peptide_chain_id} | set(receptor_chain_ids) | set(ion_chain_ids or [])
    remove_indices = [i for i, c in enumerate(fixer.topology.chains()) if c.id not in keep_ids]
    if remove_indices:
        fixer.removeChains(remove_indices)

    modeller = Modeller(fixer.topology, fixer.positions)
    waters = [r for r in modeller.topology.residues() if r.name.upper() in WATER_RESIDUES]
    if waters:
        modeller.delete(waters)
        fixer.topology = modeller.topology
        fixer.positions = modeller.positions

    # findMissingResidues must run while caps are still present in the topology.
    # If we delete the cap first, PDBFixer sees it as missing in SEQRES and re-adds
    # it during addMissingAtoms — defeating the strip entirely.
    fixer.findMissingResidues()

    chains = list(fixer.topology.chains())
    for key in list(fixer.missingResidues.keys()):
        chain_idx, res_idx = key
        if chains[chain_idx].id == peptide_chain_id:
            continue
        n_res = sum(1 for _ in chains[chain_idx].residues())
        if res_idx == 0 or res_idx >= n_res:
            del fixer.missingResidues[key]

    # Strip terminal caps after findMissingResidues so they are not catalogued as
    # missing residues and therefore not re-added by addMissingAtoms.
    if caps_to_strip:
        cap_names = set(caps_to_strip)
        pep_chain = next((c for c in modeller.topology.chains()
                          if c.id == peptide_chain_id), None)
        if pep_chain:
            pep_residues = list(pep_chain.residues())
            to_delete = []
            if pep_residues and pep_residues[0].name in cap_names:
                to_delete.append(pep_residues[0])
            if pep_residues and pep_residues[-1].name in cap_names:
                to_delete.append(pep_residues[-1])
            if to_delete:
                modeller.delete(to_delete)
                fixer.topology = modeller.topology
                fixer.positions = modeller.positions

    if convert_mse:
        fixer.findNonstandardResidues()
        fixer.replaceNonstandardResidues()

    residues_before = {c.id: sum(1 for _ in c.residues()) for c in fixer.topology.chains()}
    try:
        fixer.findMissingAtoms()
        fixer.addMissingAtoms()
    except Exception as exc:
        return False, f"addMissingAtoms error: {exc}"

    residues_added = {
        c.id: sum(1 for _ in c.residues()) - residues_before.get(c.id, 0)
        for c in fixer.topology.chains()
    }

    pos_arr = np.array(
        [[v.x, v.y, v.z] for v in fixer.positions.value_in_unit(unit.angstrom)],
        dtype=np.float32,
    )
    for chain in fixer.topology.chains():
        ca_positions = []
        for residue in chain.residues():
            ca = next((a for a in residue.atoms() if a.name == "CA"), None)
            if ca is not None:
                ca_positions.append(pos_arr[ca.index])
        for j in range(len(ca_positions) - 1):
            dist = float(np.linalg.norm(ca_positions[j + 1] - ca_positions[j]))
            if dist > 5.0:
                return False, (
                    f"Unfilled backbone gap in chain {chain.id} "
                    f"(consecutive CA distance {dist:.1f} A)"
                )

    letters = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    chain_id_map = {
        chain.id: letters[i]
        for i, chain in enumerate(fixer.topology.chains())
        if i < len(letters)
    }

    fixed_dir.mkdir(exist_ok=True)
    out_path = fixed_dir / f"{pdb_id}.cif"
    try:
        with open(out_path, "w", encoding="ascii") as fh:
            PDBxFile.writeFile(fixer.topology, fixer.positions, fh)
    except Exception as exc:
        return False, f"CIF write error: {exc}"

    return True, {"residues_added": residues_added, "chain_id_map": chain_id_map}


# ---------------------------------------------------------------------------
# Assembly selection
# ---------------------------------------------------------------------------

def get_assemblies(pdb_id, pdb_cache):
    """Return sorted list of (n_chains, asm_num, cif_path) for assemblies
    containing a plausible peptide chain."""
    pdb_cache = Path(pdb_cache)
    result = []
    for cif_path in pdb_cache.glob(f"{pdb_id}-*.cif"):
        try:
            asm_num = int(cif_path.stem.rsplit("-", 1)[1])
        except (ValueError, IndexError):
            continue
        try:
            d = MMCIF2Dict(str(cif_path))
            auth_ids   = d.get("_atom_site.auth_asym_id",    [])
            label_ids  = d.get("_atom_site.label_asym_id",   [])
            res_names  = d.get("_atom_site.auth_comp_id",    [])
            entity_ids = d.get("_atom_site.label_entity_id", [])

            chain_res_names = {}
            label_to_entity = {}
            for auth, label, name, eid in zip(auth_ids, label_ids, res_names, entity_ids):
                if name.upper() not in WATER_RESIDUES:
                    chain_res_names.setdefault(auth, set()).add(name.upper())
                if label not in label_to_entity:
                    label_to_entity[label] = eid

            auth_to_label = {}
            for auth, label in zip(auth_ids, label_ids):
                if auth not in auth_to_label:
                    auth_to_label[auth] = label

            entity_seq_len = {}
            for eid, seq in zip(
                d.get("_entity_poly.entity_id",                []),
                d.get("_entity_poly.pdbx_seq_one_letter_code", []),
            ):
                entity_seq_len[eid] = sum(1 for ch in seq if ch.isalpha())

            def chain_full_len(auth):
                eid = label_to_entity.get(auth_to_label.get(auth))
                if eid and eid in entity_seq_len:
                    return entity_seq_len[eid]
                return len(chain_res_names.get(auth, set()))

            n_mer_est = 50  # conservative; will be refined in _check_assembly
            peptide_max_len = max(50, n_mer_est * 2)
            has_peptide = any(
                is_protein_chain(chain_res_names[c]) and chain_full_len(c) <= peptide_max_len
                for c in chain_res_names
            )
            if has_peptide:
                result.append((len(chain_res_names), asm_num, cif_path))
        except Exception:
            continue
    result.sort()
    return result


# ---------------------------------------------------------------------------
# Assembly checker (core logic)
# ---------------------------------------------------------------------------

def _check_assembly(pdb_id, cif_path, n_mer, fixed_dir, convert_mse=True, strip_caps=False):
    """Evaluate one assembly CIF. Returns result dict with accepted bool and
    all metadata fields (sequences, resolution, secondary structure)."""
    raw_counts, full_seq_auth, full_seq_label = raw_chain_data(cif_path)

    try:
        fixer = PDBFixer(filename=str(cif_path))
    except Exception as exc:
        return dict(accepted=False, reject_reason=f"PDBFixer load error: {exc}")

    try:
        positions_A = np.array(
            [[v.x, v.y, v.z] for v in fixer.positions.value_in_unit(unit.angstrom)],
            dtype=np.float32,
        )
    except Exception as exc:
        return dict(accepted=False, reject_reason=f"Position extraction error: {exc}")

    topology = fixer.topology

    chains = []
    for chain in topology.chains():
        residue_names = [r.name.upper() for r in chain.residues()]
        residue_ids   = [r.id           for r in chain.residues()]
        atom_indices  = [a.index        for a in chain.atoms()]
        if not atom_indices:
            continue
        chains.append(dict(
            chain_id=chain.id,
            residue_names=residue_names,
            residue_ids=residue_ids,
            atom_indices=np.array(atom_indices, dtype=np.int32),
            is_protein=is_protein_chain(residue_names),
        ))

    protein_chains = [c for c in chains if c["is_protein"]]
    if not protein_chains:
        return dict(accepted=False, reject_reason="No protein chains in assembly")

    def effective_len(c):
        cid = c["chain_id"]
        obs = len(c["residue_names"])
        candidates = [v for v in (full_seq_auth.get(cid), full_seq_label.get(cid))
                      if v is not None and v >= obs]
        return min(candidates) if candidates else (full_seq_auth.get(cid) or full_seq_label.get(cid) or obs)

    peptide_max_len = max(50, n_mer * 2)
    peptide_candidates = [c for c in protein_chains if effective_len(c) <= peptide_max_len]
    if not peptide_candidates:
        return dict(accepted=False, reject_reason=(
            f"No protein chain with <={peptide_max_len} residues "
            f"(shortest: {min(effective_len(c) for c in protein_chains)})"
        ))

    peptide_chain = min(peptide_candidates, key=lambda c: abs(effective_len(c) - n_mer))

    if len(peptide_chain["residue_names"]) < 2:
        return dict(accepted=False, reject_reason="Peptide chain has only 1 residue")

    topology_pep = len(peptide_chain["residue_names"])
    raw_pep = raw_counts.get(
        peptide_chain["chain_id"],
        min(raw_counts.values(), key=lambda v: abs(v - topology_pep)) if raw_counts else 0,
    )
    if raw_pep > topology_pep:
        return dict(accepted=False, reject_reason=(
            f"Peptide chain {peptide_chain['chain_id']} has {raw_pep} residues in CIF "
            f"but PDBFixer loaded only {topology_pep} (non-standard residues dropped)"
        ))

    # Build exempt sets: MSE if convert_mse; caps if strip_caps (distance checked later)
    pep_exempt = set()
    if convert_mse:
        pep_exempt |= PDBFIXER_CONVERTS
    if strip_caps:
        pep_exempt |= TERMINAL_CAPS
    rec_exempt = PDBFIXER_CONVERTS if convert_mse else frozenset()

    bad = check_amber14(peptide_chain["residue_names"], exempt=pep_exempt)
    if bad:
        descs = [ptm_description(r) for r in bad]
        return dict(accepted=False, reject_reason=(
            "Non-AMBER14 residues in peptide: "
            + ", ".join(f"{r} ({d})" for r, d in zip(bad, descs))
        ))

    peptide_pos  = positions_A[peptide_chain["atom_indices"]]
    peptide_tree = cKDTree(peptide_pos)
    pep_residues = peptide_chain["residue_names"]

    receptor_chains = []
    for c in protein_chains:
        if c["chain_id"] == peptide_chain["chain_id"]:
            continue
        if peptide_similarity(pep_residues, c["residue_names"]) >= PEPTIDE_COPY_THRESHOLD:
            continue
        rec_pos = positions_A[c["atom_indices"]]
        distances, _ = peptide_tree.query(rec_pos, k=1, workers=1)
        if np.any(distances < NEIGHBOR_CUTOFF_A):
            receptor_chains.append(c)

    if not receptor_chains:
        return dict(accepted=False, reject_reason=f"No receptor chains within {NEIGHBOR_CUTOFF_A} A of peptide")

    for c in receptor_chains:
        bad = check_amber14(c["residue_names"], exempt=rec_exempt)
        if bad:
            descs = [ptm_description(r) for r in bad]
            return dict(accepted=False, reject_reason=(
                f"Non-AMBER14 residues in receptor chain {c['chain_id']}: "
                + ", ".join(f"{r} ({d})" for r, d in zip(bad, descs))
            ))

    for c in receptor_chains:
        nums = sorted(int(r) for r in c["residue_ids"] if r.lstrip("-").isdigit())
        if len(nums) >= 2:
            gap = max(nums[i + 1] - nums[i] for i in range(len(nums) - 1))
            if gap > MAX_RECEPTOR_GAP:
                return dict(accepted=False, reject_reason=(
                    f"Receptor chain {c['chain_id']} has a gap of {gap} residues "
                    f"(max allowed: {MAX_RECEPTOR_GAP})"
                ))

    # --- Terminal cap detection and interface distance ---
    caps_to_strip = []
    cap_min_dist  = None
    if strip_caps:
        pep_names = peptide_chain["residue_names"]
        rec_pos_all = np.concatenate([positions_A[c["atom_indices"]] for c in receptor_chains])
        rec_tree = cKDTree(rec_pos_all)
        topo_pep = next((c for c in topology.chains()
                         if c.id == peptide_chain["chain_id"]), None)
        pep_res_list = list(topo_pep.residues()) if topo_pep else []

        for cap, idx in [("ACE", 0), ("NME", -1), ("NH2", -1)]:
            if not pep_names or pep_names[idx] != cap:
                continue
            if not pep_res_list:
                continue
            cap_res = pep_res_list[idx]
            cap_indices = [a.index for a in cap_res.atoms()]
            if not cap_indices:
                continue
            dists, _ = rec_tree.query(positions_A[cap_indices], k=1)
            min_d = float(dists.min())
            if cap_min_dist is None or min_d < cap_min_dist:
                cap_min_dist = min_d
            if min_d <= CAP_INTERFACE_CUTOFF_A:
                return dict(accepted=False, reject_reason=(
                    f"{cap} cap at peptide interface (min dist {min_d:.1f} A)"
                ))
            caps_to_strip.append(cap)

    topo_chains = {c.id: c for c in topology.chains()}

    def terminal_atom_pos(chain_id, atom_name, end):
        tc = topo_chains.get(chain_id)
        if tc is None:
            return None
        residues = list(tc.residues())
        if not residues:
            return None
        res = residues[-1] if end == "C" else residues[0]
        for atom in res.atoms():
            if atom.name == atom_name:
                return positions_A[atom.index]
        return None

    pep_c = terminal_atom_pos(peptide_chain["chain_id"], "C", "C")
    pep_n = terminal_atom_pos(peptide_chain["chain_id"], "N", "N")
    for c in receptor_chains:
        rec_n = terminal_atom_pos(c["chain_id"], "N", "N")
        rec_c = terminal_atom_pos(c["chain_id"], "C", "C")
        if pep_c is not None and rec_n is not None:
            d = float(np.linalg.norm(pep_c - rec_n))
            if d < PEPTIDE_BOND_CUTOFF:
                return dict(accepted=False, reject_reason=(
                    f"Peptide C-terminus bonded to receptor {c['chain_id']} N-terminus ({d:.2f} A)"
                ))
        if pep_n is not None and rec_c is not None:
            d = float(np.linalg.norm(pep_n - rec_c))
            if d < PEPTIDE_BOND_CUTOFF:
                return dict(accepted=False, reject_reason=(
                    f"Peptide N-terminus bonded to receptor {c['chain_id']} C-terminus ({d:.2f} A)"
                ))

    receptor_chain_ids = [c["chain_id"] for c in receptor_chains]

    system_indices = np.concatenate(
        [peptide_chain["atom_indices"]] + [c["atom_indices"] for c in receptor_chains]
    )
    system_tree = cKDTree(positions_A[system_indices])
    ion_chain_ids = []
    for c in chains:
        if c["chain_id"] in {peptide_chain["chain_id"]} | set(receptor_chain_ids):
            continue
        if all(r in AMBER14_HETATM for r in c["residue_names"]):
            distances, _ = system_tree.query(positions_A[c["atom_indices"]], k=1, workers=1)
            if np.any(distances < NEIGHBOR_CUTOFF_A):
                ion_chain_ids.append(c["chain_id"])

    # --- Extract metadata before building fixed CIF ---
    resolution = extract_resolution(cif_path)

    # Full canonical sequence (includes missing residues PDBFixer will later fill)
    peptide_seq = (
        extract_full_sequence(cif_path, peptide_chain["chain_id"])
        or extract_sequence(peptide_chain["residue_names"])  # fallback to observed
    )
    receptor_seq_raw = {}
    for c in receptor_chains:
        seq = (
            extract_full_sequence(cif_path, c["chain_id"])
            or extract_sequence(c["residue_names"])
        )
        receptor_seq_raw[c["chain_id"]] = seq

    # Per-residue secondary structure string (H/E/C)
    peptide_ss = extract_ss_string(
        cif_path, peptide_chain["chain_id"], peptide_chain["residue_ids"]
    )
    receptor_ss_raw = {
        c["chain_id"]: extract_ss_string(cif_path, c["chain_id"], c["residue_ids"])
        for c in receptor_chains
    }

    # --- Build fixed CIF ---
    ok, build_result = build_fixed_cif(
        pdb_id, cif_path, peptide_chain["chain_id"], receptor_chain_ids, fixed_dir,
        ion_chain_ids, convert_mse=convert_mse, caps_to_strip=caps_to_strip,
    )
    if not ok:
        return dict(accepted=False, reject_reason=f"Fixed CIF build failed: {build_result}")

    fixed_cif_path = Path(fixed_dir) / f"{pdb_id}.cif"
    pep_id_check = build_result["chain_id_map"].get(
        peptide_chain["chain_id"], peptide_chain["chain_id"]
    )
    _d = MMCIF2Dict(str(fixed_cif_path))
    _pep_comps = {
        comp
        for chain, comp in zip(
            _d.get("_atom_site.auth_asym_id", []),
            _d.get("_atom_site.label_comp_id", []),
        )
        if chain == pep_id_check
    }
    _bad = sorted(_pep_comps - AMBER14_AA - WATER_RESIDUES)
    if _bad:
        fixed_cif_path.unlink(missing_ok=True)
        return dict(accepted=False,
                    reject_reason=f"Non-AMBER14 residues added by PDBFixer: {', '.join(_bad)}")

    cmap           = build_result["chain_id_map"]
    residues_added = {cmap.get(k, k): v for k, v in build_result["residues_added"].items()}
    pep_id         = cmap.get(peptide_chain["chain_id"], peptide_chain["chain_id"])
    rec_ids        = [cmap.get(cid, cid) for cid in receptor_chain_ids]

    rec_lengths = {
        cmap.get(c["chain_id"], c["chain_id"]):
            len(c["residue_names"]) + residues_added.get(cmap.get(c["chain_id"], c["chain_id"]), 0)
        for c in receptor_chains
    }
    ion_chains_out = {
        cmap.get(k, k): v
        for k, v in {
            c["chain_id"]: sorted(set(c["residue_names"]))
            for c in chains if c["chain_id"] in ion_chain_ids
        }.items()
    }

    receptor_seq_str = "|".join(
        f"{cmap.get(orig, orig)}:{seq}" for orig, seq in receptor_seq_raw.items()
    )
    receptor_ss_str = "|".join(
        f"{cmap.get(orig, orig)}:{ss}" for orig, ss in receptor_ss_raw.items()
    )

    # --- DSSP on the fixed CIF (complete coordinates, covers NMR too) ---
    fixed_cif_path = Path(fixed_dir) / f"{pdb_id}.cif"
    peptide_ss_dssp = extract_ss_dssp_chain(fixed_cif_path, pep_id)
    receptor_ss_dssp = "|".join(
        f"{cid}:{extract_ss_dssp_chain(fixed_cif_path, cid)}" for cid in rec_ids
    )

    rec_description = extract_receptor_description(cif_path, receptor_chain_ids)
    rec_func_class  = classify_receptor_function(rec_description)

    # --- Interface contacts and H-bond proxies (from fixed CIF) ---
    n_contacts, n_hbonds = compute_interface_contacts(fixed_cif_path, pep_id, rec_ids)

    all_pep_names = peptide_chain["residue_names"]
    all_rec_names = [r for c in receptor_chains for r in c["residue_names"]]
    had_mse = "MSE" in all_pep_names or "MSE" in all_rec_names

    return dict(
        accepted=True,
        peptide_chain=pep_id,
        peptide_length=len(peptide_chain["residue_names"]) + residues_added.get(pep_id, 0),
        peptide_seq=peptide_seq,
        peptide_ss=peptide_ss,
        receptor_chains=rec_ids,
        receptor_chain_lengths=rec_lengths,
        receptor_length=sum(rec_lengths.values()),
        receptor_seq=receptor_seq_str,
        receptor_ss=receptor_ss_str,
        peptide_ss_dssp=peptide_ss_dssp,
        receptor_ss_dssp=receptor_ss_dssp,
        resolution=resolution,
        receptor_description=rec_description,
        receptor_func_class=rec_func_class,
        n_contacts=n_contacts,
        n_hbonds=n_hbonds,
        ion_chains=ion_chains_out,
        residues_added=residues_added,
        mse_converted="Y" if (convert_mse and had_mse) else "",
        caps_stripped=",".join(caps_to_strip) if caps_to_strip else "",
        cap_min_dist_A=f"{cap_min_dist:.1f}" if cap_min_dist is not None else "",
        reject_reason="",
    )


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def process_entry(pdb_id, n_mer, pdb_cache, fixed_dir, convert_mse=True, strip_caps=False):
    """Select best assembly and run full filter pipeline for one entry."""
    assemblies = get_assemblies(pdb_id, pdb_cache)
    if not assemblies:
        return dict(accepted=False, assembly_num=None,
                    reject_reason="No assembly contains a peptide candidate")

    asm_summary = "  ".join(f"asm{n}:{c}" for c, n, _ in assemblies)
    n_chains, asm_num, cif_path = assemblies[0]
    print(f"    assemblies: {asm_summary}  -> using asm{asm_num}")

    result = _check_assembly(pdb_id, cif_path, n_mer, fixed_dir,
                             convert_mse=convert_mse, strip_caps=strip_caps)
    result["assembly_num"] = asm_num
    return result
