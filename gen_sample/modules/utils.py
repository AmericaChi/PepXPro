from pathlib import Path


def load_tsv(path):
    path = Path(path)
    rows = []
    with open(path, encoding="ascii") as fh:
        header = fh.readline().strip().split("\t")
        for line in fh:
            if not line.strip():
                continue
            rows.append(dict(zip(header, line.strip().split("\t"))))
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


def enrich_rows(rows):
    """Add _dg, _pep_len, _rec_len, _n_rec_chains, _resolution to each row."""
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
            r["_resolution"] = None  # NMR or NA
        valid.append(r)
    return valid
