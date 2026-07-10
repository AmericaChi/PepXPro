from pathlib import Path


def load_tsv(path):
    """Parse a tab-separated file into a list of dicts, one per row.

    Uses rstrip("\\n") rather than strip() when splitting a line: a bare
    strip() also eats trailing tabs on rows whose last column(s) are empty
    (e.g. sparse/merged TSVs with optional trailing metadata), silently
    dropping those columns from the resulting dict.
    """
    path = Path(path)
    rows = []
    with open(path, encoding="ascii") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            if not line.strip():
                continue
            rows.append(dict(zip(header, line.rstrip("\n").split("\t"))))
    return rows


def save_tsv(rows, path, columns):
    path = Path(path)
    with open(path, "w", encoding="ascii") as fh:
        fh.write("\t".join(columns) + "\n")
        for r in rows:
            fh.write("\t".join(str(r.get(c, "")) for c in columns) + "\n")


def parse_receptor(receptor_str):
    """'AB:440' → ('AB', 440)"""
    parts = receptor_str.split(":")
    return parts[0], int(parts[1])


def parse_peptide(peptide_str):
    """'C:8' → ('C', 8)"""
    parts = peptide_str.split(":")
    return parts[0], int(parts[1])


SS_CLASS_ORDER   = ["helix", "partial-helix", "beta", "mixed", "coil", "other"]
FUNC_CLASS_ORDER = ["enzyme", "chaperone", "structural", "signal",
                    "regulatory", "peptide-binding", "other"]

def enrich_rows(rows):
    """Add _dg, _pep_len, _rec_len, _n_rec_chains, _resolution, _ss_class to each row."""
    valid = []
    for r in rows:
        try:
            r["_dg"] = float(r["dG_kcal_mol"])
            _, r["_pep_len"] = parse_peptide(r["peptide"])
            chains, r["_rec_len"] = parse_receptor(r["receptor"])
            r["_n_rec_chains"] = len(chains)
        except (ValueError, KeyError, IndexError):
            continue
        # Resolution: float if numeric, None if NMR/NA/missing
        res_raw = r.get("resolution", "NA")
        try:
            r["_resolution"] = float(res_raw)
        except (ValueError, TypeError):
            r["_resolution"] = None
        # Year: int if numeric, None if NA/missing
        try:
            r["_year"] = int(r.get("year", "NA"))
        except (ValueError, TypeError):
            r["_year"] = None
        # SS class: ordinal int for sampling, None if NA
        ss_class = r.get("peptide_ss_class", "NA")
        try:
            r["_ss_class"] = SS_CLASS_ORDER.index(ss_class)
        except ValueError:
            r["_ss_class"] = None
        # Receptor functional class: ordinal int for sampling, None if NA
        func_class = r.get("receptor_func_class", "NA")
        try:
            r["_func_class"] = FUNC_CLASS_ORDER.index(func_class)
        except ValueError:
            r["_func_class"] = None
        # SS percentages (peptide)
        for src, dst in [("peptide_pct_helix", "_pep_helix"),
                          ("peptide_pct_beta",  "_pep_beta"),
                          ("peptide_pct_coil",  "_pep_coil")]:
            try:
                r[dst] = float(r.get(src) or 0)
            except (ValueError, TypeError):
                r[dst] = 0.0
        # SS percentages (receptor)
        for src, dst in [("receptor_pct_helix", "_rec_helix"),
                          ("receptor_pct_beta",  "_rec_beta"),
                          ("receptor_pct_coil",  "_rec_coil")]:
            try:
                r[dst] = float(r.get(src) or 0)
            except (ValueError, TypeError):
                r[dst] = 0.0
        # Interface features (n_contacts, n_hbonds) — None if absent
        for src, dst in [("n_contacts", "_n_contacts"), ("n_hbonds", "_n_hbonds"),
                         ("pep_sasa_burial", "_sasa_burial")]:
            raw = r.get(src, "")
            try:
                r[dst] = float(raw) if raw != "" else None
            except (ValueError, TypeError):
                r[dst] = None
        # Binding type flag (1.0 = Kd, 0.0 = Ki/IC50)
        r["_binding_kd"] = 1.0 if r.get("binding_type") == "Kd" else 0.0
        valid.append(r)
    return valid
