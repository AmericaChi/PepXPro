"""
Step 2: For each downloaded CIF, select the best biological assembly,
identify the peptide chain (shortest protein chain ≤ max(50, n_mer*2)
residues and closest to n_mer), verify AMBER14 compatibility, identify
receptor chains within 4 Å of the peptide, and write a fixed CIF with
missing residues filled (all for peptide; non-terminal only for receptor)
and waters removed.

Input:  peptide_entries.tsv, pdb_cache/<pdb_id>-<N>.cif
Output: filtered_systems.tsv, fixed_cif/<pdb_id>.cif

filtered_systems.tsv columns (accepted rows):
  pdb_id, n_mer, accepted, assembly_num,
  peptide_chain, peptide_length, peptide_added,
  receptor_chains, receptor_length, receptor_added,
  reject_reason

  peptide_length  — post-fix residue count (observed + PDBFixer-added)
  peptide_added   — residues PDBFixer added to the peptide chain
  receptor_length — total post-fix residue count across all receptor chains
  receptor_added  — total residues PDBFixer added across all receptor chains

Run:    mamba activate priml && python -u 02_filter_peptide.py 2>&1 | tee step02_run.log
"""

import shutil
import argparse
from difflib import SequenceMatcher
from pathlib import Path

from Bio.PDB.MMCIF2Dict import MMCIF2Dict

import numpy as np
from scipy.spatial import cKDTree
from pdbfixer import PDBFixer
from openmm.app import PDBxFile, Modeller
import openmm.unit as unit

PDB_CACHE    = Path("pdb_cache")
FIXED_DIR    = Path("fixed_cif")
IN_TSV       = Path("peptide_entries.tsv")
OUT_TSV      = Path("filtered_systems.tsv")

# Maximum deviation from target N-mer length to still accept a chain.

# Residue names accepted by CHARMM36 for standard MD (no patches required).
# All-caps 3-letter codes as they appear in CIF atom_site records.
AMBER14_AA = {
    "ALA", "ARG", "ASN", "ASP", "CYS",
    "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO",
    "SER", "THR", "TRP", "TYR", "VAL",
    # AMBER protonation-state variants (rarely appear in raw CIF but harmless to allow)
    "HID", "HIE", "HIP",
}

# Non-standard residues that indicate a PTM or crystallographic substitution.
# Used only for informative rejection messages.
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

NEIGHBOR_CUTOFF_A = 4.0
MAX_RECEPTOR_GAP  = 10   # reject if any receptor chain has an internal gap larger than this
PEPTIDE_BOND_CUTOFF = 1.8  # Å; C-N distance indicating covalent continuity across chain boundary


# Chains with sequence similarity to the peptide above this threshold are
# treated as additional peptide copies, not receptor chains.
PEPTIDE_COPY_THRESHOLD = 0.7

# Water residue names to strip before simulation.
WATER_RESIDUES = {"HOH", "WAT", "H2O", "TIP", "TIP3", "TIP3P", "SOL"}

# Non-AA residues with native AMBER14/TIP3P-FB parameters (common ions).
# Names must match AMBER14 residue templates in amber14/tip3pfb.xml.
AMBER14_HETATM = {
    "NA",                   # sodium
    "K",                    # potassium
    "CL",                   # chloride
    "CA",                   # calcium
    "MG",                   # magnesium
    "ZN",                   # zinc
    "FE", "FE2",            # iron
    "CU",                   # copper
    "MN",                   # manganese
}


def load_entries(tsv_path):
    entries = {}
    with open(tsv_path, encoding="ascii") as fh:
        header = fh.readline().strip().split("\t")
        for line in fh:
            parts = line.strip().split("\t")
            row = dict(zip(header, parts))
            if row.get("step1_reject"):
                continue  # already rejected
            entries[row["pdb_id"]] = int(row["n_mer"])
    return entries  # {pdb_id: n_mer}


def is_protein_chain(residues):
    """True if any residue in the chain is an amino acid.
    A low threshold avoids excluding short or heavily-modified peptides;
    non-protein chains (DNA, RNA, small molecules) have zero AA residues."""
    return any(r in AMBER14_AA or r in KNOWN_PTM for r in residues)


def peptide_similarity(residues_a, residues_b):
    """Sequence similarity ratio [0-1] between two residue name lists."""
    return SequenceMatcher(None, residues_a, residues_b).ratio()


def check_charmm36(residue_names):
    """Return list of non-CHARMM36 residue names found (empty = pass).
    Waters are ignored here; they are stripped in build_fixed_cif."""
    bad = []
    for name in residue_names:
        if name in WATER_RESIDUES:
            continue
        if name not in AMBER14_AA and name not in AMBER14_HETATM:
            bad.append(name)
    return bad


def ptm_description(resname):
    return KNOWN_PTM.get(resname, f"unknown non-standard residue {resname}")


def build_fixed_cif(pdb_id, cif_path, peptide_chain_id, receptor_chain_ids, ion_chain_ids=None):
    """
    Load the raw CIF fresh, strip non-system chains, fill missing residues
    (all for peptide; non-terminal only for receptor), write fixed_cif/<pdb_id>.cif.
    ion_chain_ids: CHARMM36-compatible ion chains to retain alongside the receptor.
    Returns (True, {chain_id: n_added}) on success or (False, error_message) on failure.
    """
    try:
        fixer = PDBFixer(filename=str(cif_path))
    except Exception as exc:
        return False, f"PDBFixer reload error: {exc}"

    keep_ids = {peptide_chain_id} | set(receptor_chain_ids) | set(ion_chain_ids or [])
    remove_indices = [
        i for i, c in enumerate(fixer.topology.chains())
        if c.id not in keep_ids
    ]
    if remove_indices:
        fixer.removeChains(remove_indices)

    # Strip water residues; keep CHARMM36-compatible HETATMs (ions, cofactors).
    modeller = Modeller(fixer.topology, fixer.positions)
    waters = [r for r in modeller.topology.residues()
              if r.name.upper() in WATER_RESIDUES]
    if waters:
        modeller.delete(waters)
        fixer.topology = modeller.topology
        fixer.positions = modeller.positions

    fixer.findMissingResidues()

    # Drop terminal missing residues for receptor chains only.
    # missingResidues keys: (chain_index, insertion_position)
    # insertion_position == 0              → N-terminal gap
    # insertion_position >= len(chain_res) → C-terminal gap
    chains = list(fixer.topology.chains())
    for key in list(fixer.missingResidues.keys()):
        chain_idx, res_idx = key
        if chains[chain_idx].id == peptide_chain_id:
            continue  # keep all missing residues for the peptide
        n_res = sum(1 for _ in chains[chain_idx].residues())
        if res_idx == 0 or res_idx >= n_res:
            del fixer.missingResidues[key]

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

    # Verify no backbone discontinuities remain (unfilled gaps leave consecutive
    # residues in the topology whose CA atoms are not at peptide-bond distance).
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

    # PDBxFile assigns sequential A, B, C... regardless of topology chain IDs.
    # Build the mapping now so callers can report the IDs that appear in the file.
    letters = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    chain_id_map = {
        chain.id: letters[i]
        for i, chain in enumerate(fixer.topology.chains())
        if i < len(letters)
    }

    FIXED_DIR.mkdir(exist_ok=True)
    out_path = FIXED_DIR / f"{pdb_id}.cif"
    try:
        with open(out_path, "w", encoding="ascii") as fh:
            PDBxFile.writeFile(fixer.topology, fixer.positions, fh)
    except Exception as exc:
        return False, f"CIF write error: {exc}"

    return True, {"residues_added": residues_added, "chain_id_map": chain_id_map}



def get_assemblies(pdb_id, n_mer):
    """Return list of (n_chains, asm_num, cif_path) for assemblies that contain
    at least one protein chain short enough to be a peptide candidate
    (full entity sequence length <= max(50, n_mer*2)), sorted by chain count
    (fewest first), then by assembly number. Assemblies with no peptide
    candidate are excluded."""
    peptide_max_len = max(50, n_mer * 2)
    result = []
    for cif_path in PDB_CACHE.glob(f"{pdb_id}-*.cif"):
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

            # chain_res_names keyed by auth_asym_id (for chain counting and is_protein check)
            chain_res_names  = {}
            label_to_entity  = {}
            for auth, label, name, eid in zip(auth_ids, label_ids, res_names, entity_ids):
                if name.upper() not in WATER_RESIDUES:
                    chain_res_names.setdefault(auth, set()).add(name.upper())
                if label not in label_to_entity:
                    label_to_entity[label] = eid

            # auth -> label (first occurrence)
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


def raw_chain_data(cif_path):
    """Parse CIF atom_site and entity tables. Returns:
      observed:        {auth_asym_id: n_non_water_residues}
      full_seq_auth:   {auth_asym_id: entity_sequence_length}
      full_seq_label:  {label_asym_id: entity_sequence_length}

    Both full_seq dicts are provided because PDBFixer inconsistently uses
    auth_asym_id or label_asym_id for chain.id depending on the CIF file.
    Callers should use the min-valid approach: for a PDBFixer chain with
    observed count obs, the correct entity length is the smallest value
    >= obs across both dicts for that chain ID.
    """
    d = MMCIF2Dict(str(cif_path))

    auth_ids   = d.get("_atom_site.auth_asym_id",    [])
    label_ids  = d.get("_atom_site.label_asym_id",   [])
    res_names  = d.get("_atom_site.auth_comp_id",    [])
    res_nums   = d.get("_atom_site.auth_seq_id",     [])
    entity_ids = d.get("_atom_site.label_entity_id", [])

    observed = {}
    auth_to_entity   = {}
    label_to_entity  = {}
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

    full_seq_auth = {auth: entity_seq_len[eid]
                     for auth, eid in auth_to_entity.items() if eid in entity_seq_len}
    full_seq_label = {label: entity_seq_len[eid]
                      for label, eid in label_to_entity.items() if eid in entity_seq_len}

    return observed, full_seq_auth, full_seq_label


def _check_assembly(pdb_id, cif_path, n_mer):
    """
    Try one assembly CIF: find peptide chain, check residues, find receptor chains.
    Returns a result dict with accepted (bool) and associated fields.
    """
    raw_counts, full_seq_auth, full_seq_label = raw_chain_data(cif_path)

    try:
        fixer = PDBFixer(filename=str(cif_path))
    except Exception as exc:
        return dict(accepted=False, reject_reason=f"PDBFixer load error: {exc}")

    try:
        positions_A = np.array(
            [[v.x, v.y, v.z]
             for v in fixer.positions.value_in_unit(unit.angstrom)],
            dtype=np.float32,
        )
    except Exception as exc:
        return dict(accepted=False, reject_reason=f"Position extraction error: {exc}")

    topology = fixer.topology

    # --- Collect chain metadata ---
    chains = []
    for chain in topology.chains():
        residue_names = [r.name.upper() for r in chain.residues()]
        atom_indices  = [a.index for a in chain.atoms()]
        if not atom_indices:
            continue
        chains.append(
            dict(
                chain_id=chain.id,
                residue_names=residue_names,
                residue_ids=[r.id for r in chain.residues()],
                atom_indices=np.array(atom_indices, dtype=np.int32),
                is_protein=is_protein_chain(residue_names),
            )
        )

    protein_chains = [c for c in chains if c["is_protein"]]
    if not protein_chains:
        return dict(accepted=False, reject_reason="No protein chains in assembly")

    # --- Find the peptide chain ---
    # Strategy: filter to chains that could plausibly be a peptide (full sequence
    # length <= max(50, n_mer*2)), then pick the one closest to n_mer.
    # Full sequence length (from _entity_poly_seq via auth->label->entity mapping)
    # is used so that peptides with missing residues are not penalised.
    # PDBFixer inconsistently uses auth or label IDs for chain.id, so we
    # check both lookups and take the smallest entity length that is >= the
    # chain's observed residue count (physically: entity >= observed always).
    def effective_len(c):
        cid = c["chain_id"]
        obs = len(c["residue_names"])
        candidates = [v for v in (full_seq_auth.get(cid), full_seq_label.get(cid))
                      if v is not None and v >= obs]
        if candidates:
            return min(candidates)
        # Fallback: any available entity length, else observed count
        return full_seq_auth.get(cid) or full_seq_label.get(cid) or obs

    peptide_max_len = max(50, n_mer * 2)
    peptide_candidates = [c for c in protein_chains if effective_len(c) <= peptide_max_len]
    if not peptide_candidates:
        return dict(
            accepted=False,
            reject_reason=(
                f"No protein chain with <={peptide_max_len} residues "
                f"(shortest: {min(effective_len(c) for c in protein_chains)})"
            ),
        )

    peptide_chain = min(peptide_candidates, key=lambda c: abs(effective_len(c) - n_mer))

    if len(peptide_chain["residue_names"]) < 2:
        return dict(accepted=False, reject_reason="Peptide chain has only 1 residue (not a peptide)")

    # --- Check PDBFixer didn't silently drop residues from the peptide ---
    # PDBFixer may use auth or label IDs; match by closest observed count.
    topology_pep = len(peptide_chain["residue_names"])
    raw_pep = raw_counts.get(peptide_chain["chain_id"],
                             min(raw_counts.values(),
                                 key=lambda v: abs(v - topology_pep)) if raw_counts else 0)
    if raw_pep > topology_pep:
        return dict(
            accepted=False,
            reject_reason=(
                f"Peptide chain {peptide_chain['chain_id']} has {raw_pep} residues in CIF "
                f"but PDBFixer loaded only {topology_pep} (non-standard residues dropped)"
            ),
        )

    # --- Check CHARMM36 compatibility of peptide chain ---
    bad_residues = check_charmm36(peptide_chain["residue_names"])
    if bad_residues:
        descs = [ptm_description(r) for r in bad_residues]
        return dict(
            accepted=False,
            reject_reason=(
                f"Non-CHARMM36 residues in peptide: "
                + ", ".join(f"{r} ({d})" for r, d in zip(bad_residues, descs))
            ),
        )

    # --- Identify receptor chains (protein chains within 4 A of peptide) ---
    peptide_pos = positions_A[peptide_chain["atom_indices"]]
    peptide_tree = cKDTree(peptide_pos)

    pep_residues = peptide_chain["residue_names"]
    receptor_chains = []
    for c in protein_chains:
        if c["chain_id"] == peptide_chain["chain_id"]:
            continue
        # Skip chains that look like additional copies of the peptide.
        sim = peptide_similarity(pep_residues, c["residue_names"])
        if sim >= PEPTIDE_COPY_THRESHOLD:
            continue
        rec_pos = positions_A[c["atom_indices"]]
        distances, _ = peptide_tree.query(rec_pos, k=1, workers=1)
        if np.any(distances < NEIGHBOR_CUTOFF_A):
            receptor_chains.append(c)

    if not receptor_chains:
        return dict(accepted=False, reject_reason=f"No receptor chains within {NEIGHBOR_CUTOFF_A} A of peptide")

    # --- Check CHARMM36 compatibility of receptor chains ---
    for c in receptor_chains:
        bad_residues = check_charmm36(c["residue_names"])
        if bad_residues:
            descs = [ptm_description(r) for r in bad_residues]
            return dict(
                accepted=False,
                reject_reason=(
                    f"Non-CHARMM36 residues in receptor chain {c['chain_id']}: "
                    + ", ".join(f"{r} ({d})" for r, d in zip(bad_residues, descs))
                ),
            )

    # --- Check for large internal gaps in receptor chains ---
    for c in receptor_chains:
        nums = sorted(int(r) for r in c["residue_ids"] if r.lstrip("-").isdigit())
        if len(nums) >= 2:
            gap = max(nums[i + 1] - nums[i] for i in range(len(nums) - 1))
            if gap > MAX_RECEPTOR_GAP:
                return dict(
                    accepted=False,
                    reject_reason=(
                        f"Receptor chain {c['chain_id']} has a gap of "
                        f"{gap} residues (max allowed: {MAX_RECEPTOR_GAP})"
                    ),
                )

    # --- Reject if peptide and receptor are covalently continuous ---
    # Some entries have a single polypeptide split across two auth chain IDs;
    # the peptide C (or N) terminus will be at peptide-bond distance from the
    # receptor N (or C) terminus.
    topo_chains = {c.id: c for c in topology.chains()}

    def terminal_atom_pos(chain_id, atom_name, end):
        """Position of atom_name in the first ('N') or last ('C') residue of chain_id."""
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
                return dict(
                    accepted=False,
                    reject_reason=(
                        f"Peptide C-terminus bonded to receptor chain {c['chain_id']} "
                        f"N-terminus ({d:.2f} A) - chains are covalently continuous"
                    ),
                )
        if pep_n is not None and rec_c is not None:
            d = float(np.linalg.norm(pep_n - rec_c))
            if d < PEPTIDE_BOND_CUTOFF:
                return dict(
                    accepted=False,
                    reject_reason=(
                        f"Peptide N-terminus bonded to receptor chain {c['chain_id']} "
                        f"C-terminus ({d:.2f} A) - chains are covalently continuous"
                    ),
                )

    receptor_length = sum(len(c["residue_names"]) for c in receptor_chains)
    receptor_chain_ids = [c["chain_id"] for c in receptor_chains]

    # --- Collect CHARMM36-compatible ion chains near the system ---
    # Ions are often deposited in their own single-residue chains and would
    # otherwise be excluded because they fail is_protein_chain.
    system_indices = np.concatenate(
        [peptide_chain["atom_indices"]]
        + [c["atom_indices"] for c in receptor_chains]
    )
    system_tree = cKDTree(positions_A[system_indices])
    ion_chain_ids = []
    for c in chains:
        if c["chain_id"] == peptide_chain["chain_id"]:
            continue
        if c["chain_id"] in receptor_chain_ids:
            continue
        if all(r in AMBER14_HETATM for r in c["residue_names"]):
            distances, _ = system_tree.query(positions_A[c["atom_indices"]], k=1, workers=1)
            if np.any(distances < NEIGHBOR_CUTOFF_A):
                ion_chain_ids.append(c["chain_id"])

    ok, result = build_fixed_cif(pdb_id, cif_path, peptide_chain["chain_id"], receptor_chain_ids, ion_chain_ids)
    if not ok:
        return dict(accepted=False, reject_reason=f"Fixed CIF build failed: {result}")

    cmap            = result["chain_id_map"]
    residues_added  = {cmap.get(k, k): v for k, v in result["residues_added"].items()}
    pep_id          = cmap.get(peptide_chain["chain_id"], peptide_chain["chain_id"])

    rec_ids         = [cmap.get(cid, cid) for cid in receptor_chain_ids]
    # Post-fix lengths (observed + added by PDBFixer).
    rec_lengths     = {cmap.get(c["chain_id"], c["chain_id"]):
                           len(c["residue_names"]) + residues_added.get(cmap.get(c["chain_id"], c["chain_id"]), 0)
                       for c in receptor_chains}
    ion_chains_out  = {cmap.get(k, k): v
                       for k, v in {c["chain_id"]: sorted(set(c["residue_names"]))
                                    for c in chains if c["chain_id"] in ion_chain_ids}.items()}

    return dict(
        accepted=True,
        peptide_chain=pep_id,
        peptide_length=len(peptide_chain["residue_names"]) + residues_added.get(pep_id, 0),
        receptor_chains=rec_ids,
        receptor_chain_lengths=rec_lengths,
        receptor_length=sum(rec_lengths.values()),
        ion_chains=ion_chains_out,
        residues_added=residues_added,
        reject_reason="",
    )


def process_entry(pdb_id, n_mer):
    """Use the assembly with the fewest chains that contains a peptide candidate."""
    assemblies = get_assemblies(pdb_id, n_mer)
    if not assemblies:
        return dict(accepted=False, assembly_num=None, reject_reason="No assembly contains a peptide candidate")

    asm_summary = "  ".join(f"asm{n}:{c}" for c, n, _ in assemblies)
    n_chains, asm_num, cif_path = assemblies[0]
    print(f"    assemblies: {asm_summary}  -> using asm{asm_num}")

    result = _check_assembly(pdb_id, cif_path, n_mer)
    result["assembly_num"] = asm_num
    return result


def main():
    parser = argparse.ArgumentParser(
        description="Step 2: Filter peptide entries, check CHARMM36 compatibility, "
                    "find receptor chains, and write fixed CIFs to fixed_cif/.")
    parser.add_argument("--clean", action="store_true",
                        help="Remove files generated by this step (filtered_systems.tsv, fixed_cif/).")
    args = parser.parse_args()

    if args.clean:
        if OUT_TSV.exists():
            OUT_TSV.unlink()
            print(f"Removed {OUT_TSV}")
        if FIXED_DIR.exists():
            shutil.rmtree(FIXED_DIR)
            print(f"Removed {FIXED_DIR}/")
        return

    print(f"Reading {IN_TSV} ...")
    entries = load_entries(IN_TSV)
    print(f"  {len(entries)} entries to process")

    out_cols = [
        "pdb_id", "n_mer", "accepted", "assembly_num",
        "peptide_chain", "peptide_length", "peptide_added",
        "receptor_chains", "receptor_length", "receptor_added",
        "reject_reason",
    ]

    accepted_count = 0
    rejected_count = 0

    with open(OUT_TSV, "w", encoding="ascii") as fh:
        fh.write("\t".join(out_cols) + "\n")

        n_total = len(entries)
        for i, (pdb_id, n_mer) in enumerate(entries.items(), 1):
            result = process_entry(pdb_id, n_mer)

            if result["accepted"]:
                rec_chains_str = ",".join(result["receptor_chains"])
                accepted_count += 1
                pep_id    = result["peptide_chain"]
                pep_added = result["residues_added"].get(pep_id, 0)
                rec_added = sum(result["residues_added"].get(cid, 0)
                                for cid in result["receptor_chains"])
                row = [
                    pdb_id,
                    str(n_mer),
                    "Y",
                    str(result["assembly_num"]),
                    result["peptide_chain"],
                    str(result["peptide_length"]),
                    str(pep_added),
                    rec_chains_str,
                    str(result["receptor_length"]),
                    str(rec_added),
                    "",
                ]
                pep_str   = f"{pep_id}:{result['peptide_length'] - pep_added}+{pep_added}"
                rec_str   = "  ".join(
                    f"{cid}:{result['receptor_chain_lengths'][cid] - result['residues_added'].get(cid, 0)}"
                    f"+{result['residues_added'].get(cid, 0)}"
                    for cid in result["receptor_chains"]
                )
                ions_all = sorted(set(
                    r for resnames in result["ion_chains"].values() for r in resnames
                ))
                ions_str = f"  ions=[{','.join(ions_all)}]" if ions_all else ""
                print(f"  [{i}/{n_total}] PASS   {pdb_id}  asm={result['assembly_num']}  "
                      f"pep={pep_str}  rec={rec_str}{ions_str}  "
                      f"(accepted={accepted_count} rejected={rejected_count})")
                if result["peptide_length"] != n_mer:
                    print(f"    WARNING: peptide length {result['peptide_length']} != n_mer {n_mer}")
            else:
                row = [pdb_id, str(n_mer), "N", str(result.get("assembly_num") or ""), "", "", "", "", "", "", result["reject_reason"]]
                rejected_count += 1
                print(f"  [{i}/{n_total}] REJECT {pdb_id}: {result['reject_reason']}  "
                      f"(accepted={accepted_count} rejected={rejected_count})")

            fh.write("\t".join(row) + "\n")

    print(
        f"\nDone. accepted={accepted_count} rejected={rejected_count} "
        f"total={accepted_count + rejected_count}"
    )
    print(f"Written {OUT_TSV}")


if __name__ == "__main__":
    main()
