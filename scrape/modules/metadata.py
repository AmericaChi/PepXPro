"""
metadata.py — Extract structural metadata from CIF files.

Functions:
  extract_sequence(residue_names)                    → 1-letter string (observed only)
  extract_full_sequence(cif_path, chain_id)          → 1-letter string (complete canonical)
  extract_resolution(cif_path)                       → '1.80' or 'NA'
  extract_ss_string(cif_path, chain_id, residue_ids) → per-residue 'CCCHHHHCCEE' or 'NA'
  extract_receptor_entity(cif_path, receptor_seqs)   → (description, ec_number)
"""

import re

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


def extract_exptl_method(cif_path):
    """Return experimental method string from CIF, e.g. 'SOLUTION NMR', 'X-RAY DIFFRACTION'.
    Returns 'NA' if not found."""
    from Bio.PDB.MMCIF2Dict import MMCIF2Dict
    try:
        d = MMCIF2Dict(str(cif_path))
        val = d.get("_exptl.method")
        if val:
            v = val[0] if isinstance(val, list) else val
            return str(v).strip()
    except Exception:
        pass
    return "NA"


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

    Runs mkdssp directly on the cached fixed CIF (post-PDBFixer, complete coordinates).
    Returns a per-residue string 'CCCHHHHCCEE' (3-state: H/E/C) sorted by
    residue number, or 'NA' on failure.
    Covers all experimental methods including NMR.
    """
    from Bio.PDB import MMCIFParser, DSSP as BiopythonDSSP

    try:
        parser = MMCIFParser(QUIET=True)
        structure = parser.get_structure("x", str(fixed_cif_path))
        model = structure[0]
        dssp = BiopythonDSSP(model, str(fixed_cif_path), dssp="mkdssp")
    except Exception:
        return "NA"

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


_FUNC_RULES = [
    # Immune first: MHC/HLA/antibody chains are never annotated with an EC number and
    # would otherwise be mis-caught by generic rules ("binding protein", "receptor").
    ("immune",          ["histocompatibility antigen", "mhc class", "hla class",
                         "h-2 class", "major histocompatibility", "antigen-presenting",
                         "t-cell receptor", "t cell receptor", "tcr ",
                         "beta-2-microglobulin", "beta-2 microglobulin",
                         "immunoglobulin", "antibody", "fab ", "fab fragment",
                         "fab' ", "nanobody", "tapasin"]),
    ("enzyme",          ["kinase", "phosphatase", "protease", "peptidase", "proteinase",
                         "transferase", "synthase", "synthetase", "ligase", "lyase",
                         "isomerase", "oxidase", "reductase", "dehydrogenase",
                         "hydrolase", "hydroxylase", "oxygenase", "lipase", "nuclease",
                         "polymerase", "methyltransferase", "demethylase",
                         "acetyltransferase", "deacetylase", "acetylase",
                         "caspase", "trypsin", "chymotrypsin", "elastase", "thrombin",
                         "subtilisin", "papain", "cathepsin", "calpain", "granzyme",
                         "renin", "carboxylase", "decarboxylase", "esterase",
                         "glycosylase", "topoisomerase", "helicase", "gyrase",
                         "cyclase", "phosphorylase", "dismutase", "lactamase",
                         "amidase", "deaminase", "aminotransferase", "racemase",
                         "mutase", "epimerase", "ubiquitin-conjugating",
                         "ubiquitin ligase", "e3 ubiquitin", "sumo ligase",
                         " enzyme", "catalytic subunit"]),
    ("chaperone",       ["chaperone", "chaperonin", "heat shock protein",
                         "heat shock cognate", "hsp70", "hsp90",
                         "dnak", "groel", "groes", "trigger factor", "prefoldin",
                         "peptidyl-prolyl cis-trans isomerase", "calnexin",
                         "calreticulin"]),
    ("transport",       ["importin", "exportin", "karyopherin", "coatomer", "clathrin",
                         "ap-2 complex", "ap-1 complex", "adaptor protein complex",
                         "adapter-related protein complex",
                         "vacuolar protein sorting", "escrt", "sorting nexin",
                         "tumor susceptibility gene 101", "nuclear pore", "nucleoporin",
                         "transport protein", "kinesin", "dynein", "snare", "syntaxin",
                         "nuclear rna export"]),
    ("structural",      ["structural protein", "collagen", "actin", "tubulin",
                         "myosin", "keratin", "laminin", "fibronectin",
                         "cytoskeletal", "filamin", "spectrin", "radixin", "moesin",
                         "ezrin", "capsid", "coat protein", "viral envelope",
                         "titin", "vinculin", "talin", "desmoplakin", "plectin",
                         "septin", "intermediate filament"]),
    ("regulatory",      ["transcription factor", "transcription activator",
                         "transcription repressor", "regulatory protein",
                         "transcriptional regulator", "transcriptional coactivator",
                         "transcriptional corepressor", "transcription elongation",
                         "elongation regulator", "repressor protein",
                         "tumor suppressor", "cyclin",
                         "androgen receptor", "estrogen receptor", "oestrogen receptor",
                         "glucocorticoid receptor", "retinoic acid receptor",
                         "thyroid hormone receptor", "nuclear receptor",
                         "peroxisome proliferator-activated receptor",
                         "vitamin d receptor", "bile acid receptor",
                         "b-cell lymphoma 6", "bcl-6", "mdm2", "mdm4",
                         "retinoblastoma", "telomeric repeat-binding", "e2f",
                         "histone", "chromatin", "nucleosome", "polycomb",
                         "sin3", "steroid receptor coactivator", "splicing factor"]),
    ("signal",          ["signal transduction", "signaling protein", "gtpase",
                         "g-protein", "g protein", "ras protein", "rho protein",
                         "growth factor receptor", "receptor tyrosine kinase",
                         "g-protein coupled receptor", "g protein-coupled receptor",
                         "calmodulin", "s100", "s-100 protein", "calcineurin",
                         "arrestin", "guanine nucleotide", "chemotaxis protein",
                         "response regulator", "substance-p receptor",
                         "opioid receptor", "chemokine receptor",
                         "adrenergic receptor", "muscarinic", "integrin",
                         "troponin"]),
    ("peptide-binding", ["binding protein", "pdz domain", "sh2 domain", "sh3 domain",
                         "bromodomain", "chromodomain", "14-3-3", "recognition module",
                         "wd repeat", "wd-repeat", "wd40", "disks large", "harmonin",
                         "scribble", "syntenin", "nherf",
                         "na(+)/h(+) exchange regulatory",
                         "adaptor protein", "adaptor", "grb2", "scaffold",
                         "ankyrin repeat", "tpr repeat", "tetratricopeptide",
                         "leucine-rich repeat", "speckle-type poz", "spop",
                         "kelch-like", "kelch", "f-box", "baculoviral iap repeat",
                         "inhibitor of apoptosis", "bir domain", "ww domain",
                         "eh domain", "ph domain", "sam domain", "mit domain",
                         "uba domain", "ubiquitin-like", "ubiquitin receptor",
                         "dynein light chain", "proliferating cell nuclear antigen",
                         "sliding clamp", "dna polymerase iii subunit beta",
                         "gephyrin", "streptavidin", "avidin", "profilin"]),
]

_EC_RE = re.compile(r"^\d+\.\d+\.\d+\.[\d-]+$")


def classify_receptor_function(description, ec=None):
    """Classify receptor into one of 8 functional classes from its PDB entity description.

    Classes (priority order): immune, enzyme, chaperone, transport, structural,
    regulatory, signal, peptide-binding, other.

    An EC number (_entity.pdbx_ec) is authoritative for 'enzyme': it catches enzymes
    named only by their proper name (trypsin, caspase-7, ...) that no generic keyword
    would match. It is checked after the immune rule only because MHC/antibody chains
    never carry one.

    Entries with no usable description fall into 'other' (no separate 'NA' class).
    """
    if not description or description == "NA":
        return "other"
    desc = description.lower()
    ec_ok = bool(ec) and any(_EC_RE.match(e.strip()) for e in str(ec).split(";"))
    for cls, keywords in _FUNC_RULES:
        if cls == "enzyme" and ec_ok:
            return "enzyme"
        if any(kw in desc for kw in keywords):
            return cls
    return "other"


def extract_deposition_year(cif_path):
    """Extract publication/deposition year from an assembly CIF.

    Tries _citation.year first (present in all assembly CIFs).
    Falls back to _pdbx_database_status.recvd_initial_deposition_date.
    Returns 4-digit year string ('2010') or 'NA'.
    """
    import re as _re
    from Bio.PDB.MMCIF2Dict import MMCIF2Dict
    try:
        d = MMCIF2Dict(str(cif_path))
        for key in ("_citation.year",
                    "_pdbx_database_status.recvd_initial_deposition_date"):
            val = d.get(key)
            if not val:
                continue
            v = val[0] if isinstance(val, list) else val
            m = _re.search(r'\b(19|20)\d{2}\b', str(v))
            if m:
                return m.group(0)
    except Exception:
        pass
    return "NA"


def fetch_resolution_rcsb(pdb_id):
    """Query RCSB REST API for crystallographic resolution.

    Used as fallback when assembly CIFs lack _refine.ls_d_res_high.
    Returns '1.80' or 'NA'. Makes one HTTP request per call.
    """
    import urllib.request
    import json
    url = f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id.lower()}"
    try:
        with urllib.request.urlopen(url, timeout=10) as resp:
            data = json.loads(resp.read())
        # rcsb_entry_info.resolution_combined is a list; take the first
        combined = data.get("rcsb_entry_info", {}).get("resolution_combined")
        if combined:
            return f"{float(combined[0]):.2f}"
        # Fallback: refine block
        refine = data.get("refine")
        if refine:
            val = refine[0].get("ls_d_res_high")
            if val is not None:
                return f"{float(val):.2f}"
    except Exception:
        pass
    return "NA"


_ENTITY_MATCH_MIN = 0.8   # min identity to accept an entity as a receptor chain


def extract_receptor_entity(cif_path, receptor_seqs):
    """Return (description, ec_number) for the receptor entities of a complex.

    Receptor entities are identified by matching their canonical sequence against the
    receptor chain sequences, rather than by chain ID: PDB author chain IDs (H/L/P for
    a Fab, for instance) do not survive the pipeline's chain renaming, so a chain-ID
    lookup silently misses the entity and reports no description at all.

    receptor_seqs: list of 1-letter receptor chain sequences (order irrelevant).
    Multi-entity receptors yield all distinct descriptions, joined with '; '.
    Returns ('NA', '') when no entity can be matched.
    """
    from Bio.PDB.MMCIF2Dict import MMCIF2Dict
    try:
        d = MMCIF2Dict(str(cif_path))
    except Exception:
        return "NA", ""

    ent_ids   = d.get("_entity.id", [])
    ent_types = d.get("_entity.type", [])
    ent_descs = d.get("_entity.pdbx_description", [])
    ent_ecs   = d.get("_entity.pdbx_ec", []) or ["?"] * len(ent_ids)

    entities = {
        eid: {"desc": desc, "ec": ec, "seq": ""}
        for eid, typ, desc, ec in zip(ent_ids, ent_types, ent_descs, ent_ecs)
        if typ == "polymer"
    }
    if not entities:
        return "NA", ""

    seqs = {}
    for eid, mon in zip(d.get("_entity_poly_seq.entity_id", []),
                        d.get("_entity_poly_seq.mon_id", [])):
        seqs.setdefault(eid, []).append(AA3_TO_1.get(mon, "X"))
    for eid in entities:
        entities[eid]["seq"] = "".join(seqs.get(eid, []))

    matched = []
    for rseq in receptor_seqs:
        if not rseq:
            continue
        best, best_score = None, 0.0
        for eid, ent in entities.items():
            eseq = ent["seq"]
            if not eseq:
                continue
            # PDBFixer may have filled gaps, so the pipeline sequence can differ in
            # length from the deposited one: accept containment either way, else fall
            # back to identity over the common prefix.
            if eseq in rseq or rseq in eseq:
                score = 1.0
            else:
                n = min(len(eseq), len(rseq))
                score = (sum(a == b for a, b in zip(eseq[:n], rseq[:n])) / n) if n else 0.0
            if score > best_score:
                best, best_score = eid, score
        if best is not None and best_score >= _ENTITY_MATCH_MIN and best not in matched:
            matched.append(best)

    descs, ecs = [], []
    for eid in matched:
        desc, ec = entities[eid]["desc"], entities[eid]["ec"]
        if desc and desc not in descs:
            descs.append(desc)
        if ec and ec not in ("?", ".", "") and ec not in ecs:
            ecs.append(ec)

    return ("; ".join(descs) if descs else "NA"), ";".join(ecs)


SS_MIXED_MIN = 0.20   # each of H/E must reach this to call a peptide mixed (HE)
SS_MAJOR_MIN = 0.50   # H/E fraction at or above this is a full helix (H) / sheet (E)


def classify_peptide_ss(pct_helix_str, pct_beta_str, pct_coil_str=""):
    """Classify peptide secondary structure into 7 categories.

    Categories (priority order):
      coil          (C)   only coil residues: no helix and no strand
      mixed         (HE)  pct_helix >= 0.20 AND pct_beta >= 0.20
      helix         (H)   pct_helix >= 0.50, strand below the mixed floor
      beta          (E)   pct_beta  >= 0.50, helix below the mixed floor
      partial-helix (pH)  helix present but < 0.50, and helix-dominant
      partial-beta  (pE)  strand present but < 0.50, and strand-dominant
      other               none of the above

    A trace of the opposite element (below the 20% mixed floor) does not disqualify a
    peptide from H/pH/E/pE: it is treated as noise rather than as mixed character, so a
    60% helix with a 7% strand bridge is still a helix. The partial classes break ties
    by whichever element dominates, so a 5%-helix / 16%-strand peptide is partial-beta
    rather than partial-helix.

    Returns one of: 'helix', 'partial-helix', 'beta', 'partial-beta', 'mixed',
    'coil', 'other', 'NA'.
    """
    try:
        h = float(pct_helix_str)
        e = float(pct_beta_str)
    except (ValueError, TypeError):
        return "NA"

    if h == 0 and e == 0:
        return "coil"
    if h >= SS_MIXED_MIN and e >= SS_MIXED_MIN:
        return "mixed"
    if h >= SS_MAJOR_MIN and e < SS_MIXED_MIN:
        return "helix"
    if e >= SS_MAJOR_MIN and h < SS_MIXED_MIN:
        return "beta"
    if h > 0 and h >= e and e < SS_MIXED_MIN:
        return "partial-helix"
    if e > 0 and e > h and h < SS_MIXED_MIN:
        return "partial-beta"
    return "other"


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
