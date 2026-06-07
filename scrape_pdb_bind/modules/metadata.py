"""
metadata.py — Extract structural metadata from CIF files.

Functions:
  extract_sequence(residue_names)                    → 1-letter string (observed only)
  extract_full_sequence(cif_path, chain_id)          → 1-letter string (complete canonical)
  extract_resolution(cif_path)                       → '1.80' or 'NA'
  extract_ss_string(cif_path, chain_id, residue_ids) → per-residue 'CCCHHHHCCEE' or 'NA'
"""

DSSP_3STATE = {
    'H': 'H', 'G': 'H', 'I': 'H',   # helix variants
    'E': 'E', 'B': 'E',               # strand variants
    '-': 'C', 'T': 'C', 'S': 'C', ' ': 'C',  # coil/turns
}

AA3_TO_1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    "HID": "H", "HIE": "H", "HIP": "H",
    "MSE": "M", "SEC": "C",
}


def extract_sequence(residue_names):
    """Convert observed 3-letter residue names to 1-letter string. Unknown → 'X'."""
    return "".join(AA3_TO_1.get(r.upper(), "X") for r in residue_names)


def extract_full_sequence(cif_path, chain_id):
    """Read the complete canonical sequence for chain_id from _entity_poly.
    Includes residues missing from the deposited coordinates.
    Returns 1-letter string, or None if not available.
    """
    from Bio.PDB.MMCIF2Dict import MMCIF2Dict
    try:
        d = MMCIF2Dict(str(cif_path))
    except Exception:
        return None

    auth_ids   = d.get("_atom_site.auth_asym_id",    [])
    label_ids  = d.get("_atom_site.label_asym_id",   [])
    entity_ids = d.get("_atom_site.label_entity_id", [])

    # Build chain_id → entity_id (try auth then label)
    auth_to_entity  = {}
    label_to_entity = {}
    for auth, label, eid in zip(auth_ids, label_ids, entity_ids):
        if auth not in auth_to_entity:
            auth_to_entity[auth] = eid
        if label not in label_to_entity:
            label_to_entity[label] = eid

    eid = auth_to_entity.get(chain_id) or label_to_entity.get(chain_id)
    if not eid:
        return None

    for poly_eid, seq in zip(
        d.get("_entity_poly.entity_id",                []),
        d.get("_entity_poly.pdbx_seq_one_letter_code", []),
    ):
        if poly_eid == eid:
            # Strip whitespace/newlines inserted for formatting
            return "".join(c for c in seq if c.isalpha())

    return None


def extract_resolution(cif_path):
    """Read crystallographic resolution from CIF.
    Returns '1.80', or 'NA' for NMR/missing data."""
    from Bio.PDB.MMCIF2Dict import MMCIF2Dict
    try:
        d = MMCIF2Dict(str(cif_path))
        for key in ("_refine.ls_d_res_high",
                    "_reflns.d_resolution_high",
                    "_refine_hist.d_res_high"):
            val = d.get(key)
            if val:
                v = val[0] if isinstance(val, list) else val
                try:
                    return f"{float(v):.2f}"
                except (ValueError, TypeError):
                    pass
    except Exception:
        pass
    return "NA"


def extract_ss_dssp_chain(fixed_cif_path, chain_id):
    """Compute per-residue secondary structure using DSSP (Biopython + mkdssp).

    Runs on a fixed CIF (post-PDBFixer, complete coordinates).
    mkdssp requires PDB format — the fixed CIF is converted to a temp PDB
    in memory before calling DSSP, then the temp file is removed.
    Returns a per-residue string 'CCCHHHHCCEE' (3-state: H/E/C) sorted by
    residue number, or 'NA' on failure.
    Covers all experimental methods including NMR.
    """
    import os
    import tempfile
    from Bio.PDB import MMCIFParser, MMCIFIO, DSSP as BiopythonDSSP

    try:
        parser = MMCIFParser(QUIET=True)
        structure = parser.get_structure("x", str(fixed_cif_path))
        model = structure[0]

        # mkdssp v4 only reads mmCIF — write a temp CIF (not PDB)
        with tempfile.NamedTemporaryFile(suffix=".cif", mode="w",
                                         delete=False, encoding="ascii") as tmp:
            tmp_path = tmp.name
        io = MMCIFIO()
        io.set_structure(structure)
        io.save(tmp_path)

        dssp = BiopythonDSSP(model, tmp_path, dssp="mkdssp")
    except Exception:
        return "NA"
    finally:
        try:
            os.unlink(tmp_path)
        except Exception:
            pass

    residues = []
    for key in dssp.keys():
        c_id, res_id = key
        if c_id != chain_id:
            continue
        seqnum = res_id[1]
        ss_char = DSSP_3STATE.get(dssp[key][2], "C")
        residues.append((seqnum, ss_char))

    if not residues:
        return "NA"

    residues.sort(key=lambda x: x[0])
    return "".join(ss for _, ss in residues)


def ss_percentages(ss_str):
    """Compute H/E/C fractions from a per-residue SS string.

    Handles multi-chain format 'A:CCCHHH...|B:EEECC...'.
    Returns (pct_helix, pct_beta, pct_coil) as '0.30' strings, or 'NA' if unavailable.
    """
    if not ss_str or ss_str == "NA":
        return "NA", "NA", "NA"
    if "|" in ss_str:
        chars = "".join(p.split(":", 1)[1] for p in ss_str.split("|") if ":" in p)
    else:
        chars = ss_str
    n = len(chars)
    if n == 0:
        return "NA", "NA", "NA"
    return (
        f"{chars.count('H') / n:.2f}",
        f"{chars.count('E') / n:.2f}",
        f"{chars.count('C') / n:.2f}",
    )


def extract_ss_string(cif_path, chain_id, residue_ids):
    """Return a per-residue secondary structure string for chain_id.

    Each character maps to one residue in residue_ids order:
      H = helix   (from _struct_conf, conf_type HELX_P)
      E = strand  (from _struct_sheet_range)
      C = coil    (all others)

    Returns string of the same length as residue_ids, or 'NA' if no
    secondary structure records exist in the CIF.
    residue_ids: ordered list of auth_seq_id strings from the topology.
    """
    from Bio.PDB.MMCIF2Dict import MMCIF2Dict
    try:
        d = MMCIF2Dict(str(cif_path))
    except Exception:
        return "NA"

    n = len(residue_ids)
    if n == 0:
        return "NA"

    # Build {residue_int: ss_char} — default to coil
    # For residues with insertion codes, int() will fail; they stay 'C'
    res_to_int = {}
    for rid in residue_ids:
        try:
            res_to_int[rid] = int(rid)
        except (ValueError, TypeError):
            res_to_int[rid] = None

    ss_map = {}  # int_residue_num → 'H' or 'E'

    # Helices
    conf_chains = d.get("_struct_conf.beg_auth_asym_id", [])
    conf_beg    = d.get("_struct_conf.beg_auth_seq_id",   [])
    conf_end    = d.get("_struct_conf.end_auth_seq_id",   [])
    conf_types  = d.get("_struct_conf.conf_type_id",      [])
    has_conf    = bool(conf_chains)

    for c, b, e, t in zip(conf_chains, conf_beg, conf_end, conf_types):
        if c != chain_id or "HELX" not in t.upper():
            continue
        try:
            for r in range(int(b), int(e) + 1):
                ss_map[r] = "H"
        except (ValueError, TypeError):
            pass

    # Strands
    sh_chains = d.get("_struct_sheet_range.beg_auth_asym_id", [])
    sh_beg    = d.get("_struct_sheet_range.beg_auth_seq_id",   [])
    sh_end    = d.get("_struct_sheet_range.end_auth_seq_id",   [])
    has_sheet = bool(sh_chains)

    for c, b, e in zip(sh_chains, sh_beg, sh_end):
        if c != chain_id:
            continue
        try:
            for r in range(int(b), int(e) + 1):
                if r not in ss_map:  # helices take priority
                    ss_map[r] = "E"
        except (ValueError, TypeError):
            pass

    if not has_conf and not has_sheet:
        return "NA"

    # Build per-residue string in topology order
    result = []
    for rid in residue_ids:
        r_int = res_to_int.get(rid)
        if r_int is not None:
            result.append(ss_map.get(r_int, "C"))
        else:
            result.append("C")  # insertion code residue → coil

    return "".join(result)
