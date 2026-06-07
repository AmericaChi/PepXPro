from pathlib import Path


def load_tsv(path):
    rows = []
    with open(path, encoding="ascii") as fh:
        header = fh.readline().strip().split("\t")
        for line in fh:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            rows.append(dict(zip(header, parts)))
    return rows


def save_tsv(rows, path, columns):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="ascii") as fh:
        fh.write("\t".join(columns) + "\n")
        for row in rows:
            fh.write("\t".join(str(row.get(c, "")) for c in columns) + "\n")
