"""
05_external.py — Integrate an external peptide-protein dataset into PepXPro.

Parses an external file (xlsx/tsv/csv), normalizes affinity data,
deduplicates against PDBbind, downloads CIFs, runs the filter+summary
pipeline, and produces a merged summary ready for genSample.

Usage:
  python 05_external.py --input inputs/external/PEPBI.xlsx
  python 05_external.py --input inputs/external/dataset.tsv --name mydata
  python 05_external.py --input inputs/external/PEPBI.xlsx --no-merge

Outputs (in outputs/):
  peptide_entries_<name>.tsv   — normalized entries (input for 02_filter.py)
  filtered_systems_<name>.tsv  — structure filter output
  summary_<name>.tsv           — summary for this dataset
  summary_merged.tsv           — PDBbind + external (ready for genSample)

Column map is saved to inputs/<name>_colmap.json and reused on re-runs.
"""

import argparse
import json
import math
import re
import subprocess
import sys
from pathlib import Path

_here = Path(__file__).resolve().parent
_root = _here.parent
sys.path.insert(0, str(_root))

from modules.utils import load_tsv, save_tsv, setup_logging
from modules.download import download_entries

R_KCAL          = 0.001987   # kcal / (mol * K)
DEFAULT_TEMP_K  = 300.0      # K

OUTPUTS_DIR = _root / "outputs"
LOGS_DIR    = _root / "logs"
PDB_CACHE   = _root / "pdb_cache"
LOG_FILE    = LOGS_DIR / "05_external.log"

# Columns consumed from entries TSV by standard pipeline — not treated as "extra"
_STANDARD_ENTRY_COLS = frozenset({
    "pdb_id", "n_mer", "year", "binding_type",
    "dG_kcal", "dG_approx", "step1_reject", "resolution",
})

# ── Column aliases ────────────────────────────────────────────────────────────

FIELD_ALIASES = {
    "pdb_id": [
        "pdb_id", "PDB ID", "PDB", "pdb", "PDBID", "pdb_code", "entry",
        "accession", "pdb code", "complex", "id", "PDB Code",
    ],
    "dg_exp": [
        "dG_exp", "DG_exp", "deltaG_exp", "dG_experimental",
        "ΔG (kcal/mol)",          # ΔG (kcal/mol) — PEPBI ITC column
        "DeltaG", "exp_dG", "Exp. dG",
        "experimental deltaG", "dG(kcal/mol)", "dG (kcal/mol)",
    ],
    "dg_calc": [
        "dG_calc", "DG_calc", "deltaG_calc", "calculated_dG",
        "Calculated ΔG (kcal/mol)",   # Calculated ΔG (kcal/mol) — PEPBI
        "dG_predicted", "calc_dG", "CalcDG", "predicted dG",
        "Calculated dG (kcal/mol)",
    ],
    "kd": [
        "Kd", "kd", "KD", "k_d", "dissociation_constant",
        "KD (M)", "Kd(M)", "Kd (M)", "Kd(nM)", "kd_exp", "Kd_exp",
        "KD(M)",
    ],
    "ki": [
        "Ki", "ki", "KI", "k_i", "inhibition_constant",
        "Ki(M)", "Ki (M)", "Ki(nM)",
    ],
    "n_mer": [
        "n_mer", "peptide_length", "pep_len", "nmer", "mer",
        "Peptide Length", "peptide_len", "peplen",
    ],
    "peptide_seq": [
        "peptide_seq", "Peptide Sequence", "peptide_sequence",
        "sequence", "pep_seq", "seq", "Sequence",
    ],
    "resolution": [
        "resolution", "Resolution", "res", "resolution_A",
        "resolution(A)", "Resolution (A)",
    ],
    "year": [
        "year", "Year", "release_year", "deposition_year",
        "deposit_year", "PDB_year",
    ],
    "binding_type": [
        "binding_type", "measurement_type", "type", "affinity_type",
    ],
    "temperature": [
        "temperature", "Temperature", "T (K)", "temp_K", "T(K)",
        "temp", "Temp", "Temperature (K)", "T", "temp_celsius",
    ],
}


# ── File loading ──────────────────────────────────────────────────────────────

def _unnamed_fraction(columns):
    """Fraction of columns that look like pandas unnamed placeholders."""
    return sum(1 for c in columns if str(c).startswith("Unnamed:")) / max(len(columns), 1)


def _load_file(path):
    """Load xlsx/tsv/csv → list of dicts with string values."""
    path   = Path(path)
    suffix = path.suffix.lower()

    if suffix == ".xlsx":
        try:
            import pandas as pd
        except ImportError:
            print("  Error: pandas required for xlsx.  conda install pandas")
            sys.exit(1)
        # Detect whether we need to skip merged-cell header rows
        df_sample = pd.read_excel(path, header=0, dtype=str, nrows=3)
        if _unnamed_fraction(df_sample.columns) > 0.3:
            # Multi-row header (e.g. PEPBI): the real column names are in row 1
            df = pd.read_excel(path, header=1, dtype=str)
        else:
            df = pd.read_excel(path, header=0, dtype=str)
        df.fillna("", inplace=True)
        return df.to_dict(orient="records")

    if suffix == ".csv":
        import csv
        with open(path, encoding="utf-8-sig") as fh:
            return [dict(r) for r in csv.DictReader(fh)]

    if suffix in (".tsv", ".txt"):
        return load_tsv(path)

    # Unknown extension — try TSV
    try:
        return load_tsv(path)
    except Exception:
        print(f"  Error: unsupported format '{suffix}'. Use xlsx, csv, or tsv.")
        sys.exit(1)


# ── Column mapping ────────────────────────────────────────────────────────────

def _auto_map(columns):
    """Match column names against FIELD_ALIASES (case-insensitive).
    Returns {field: matched_col_name}.
    """
    col_lower = {c.lower().strip(): c for c in columns}
    found = {}
    for field, aliases in FIELD_ALIASES.items():
        for alias in aliases:
            candidate = col_lower.get(alias.lower().strip())
            if candidate is not None:
                found[field] = candidate
                break
    return found


def _prompt_col(field, columns, required=False):
    """Prompt user to pick a column. Returns column name or None."""
    # Write directly to the real terminal — sys.stdout may be redirected to the log.
    out = sys.__stdout__
    avail = [c for c in columns if c]
    out.write(f"\n  Column for '{field}' not detected.\n")
    for i, c in enumerate(avail, 1):
        out.write(f"    {i:2}. {c}\n")
    out.flush()
    prompt = f"  Enter number, column name, or {'REQUIRED — cannot skip' if required else 'Enter to skip'}: "
    while True:
        out.write(prompt)
        out.flush()
        resp = sys.stdin.readline().strip()
        if not resp:
            if required:
                out.write("  This field is required.\n")
                out.flush()
                continue
            return None
        try:
            idx = int(resp) - 1
            if 0 <= idx < len(avail):
                return avail[idx]
        except ValueError:
            pass
        if resp in columns:
            return resp
        out.write(f"  '{resp}' not found.\n")
        out.flush()


def _resolve_colmap(columns, name, colmap_path):
    """
    Load saved column map if valid; otherwise auto-detect and prompt
    for unresolved fields. Persists the resolved map for future runs.
    """
    if colmap_path.exists():
        with open(colmap_path, encoding="utf-8") as fh:
            saved = json.load(fh)
        invalid = [f for f, c in saved.items() if c and c not in columns]
        if not invalid:
            # Enrich saved colmap with any newly-recognised fields
            new_fields = _auto_map(columns)
            added = {f: c for f, c in new_fields.items() if f not in saved}
            if added:
                saved.update(added)
                with open(colmap_path, "w", encoding="utf-8") as fh:
                    json.dump(saved, fh, indent=2, ensure_ascii=False)
                print(f"  Column map loaded from {colmap_path.name} (added: {list(added.keys())})")
            else:
                print(f"  Column map loaded from {colmap_path.name}")
            return saved
        print(f"  Saved map references unknown columns ({invalid}) — re-detecting.")

    colmap = _auto_map(columns)
    print(f"  Auto-detected fields: {list(colmap.keys())}")

    # pdb_id is always required
    if "pdb_id" not in colmap:
        colmap["pdb_id"] = _prompt_col("pdb_id", columns, required=True)

    # Need at least one affinity source
    if not any(f in colmap for f in ("dg_exp", "dg_calc", "kd", "ki")):
        print("\n  No affinity column detected. Checking manually ...")
        for field in ("dg_exp", "kd", "ki", "dg_calc"):
            col = _prompt_col(field, columns)
            if col:
                colmap[field] = col
                break

    # Optional fields — prompt if not detected
    for field in ("n_mer", "peptide_seq", "resolution", "year", "temperature"):
        if field not in colmap:
            colmap[field] = _prompt_col(field, columns)

    colmap_path.parent.mkdir(parents=True, exist_ok=True)
    with open(colmap_path, "w", encoding="utf-8") as fh:
        json.dump(colmap, fh, indent=2, ensure_ascii=False)
    print(f"  Column map saved -> {colmap_path.name}")
    return colmap


# ── Affinity parsing ──────────────────────────────────────────────────────────

_VAL_UNIT_PAT = re.compile(
    r'^\s*([\d.eE+\-]+)\s*(pM|nM|[uµ]M|mM|M)?\s*$', re.IGNORECASE
)
_UNIT_FACTORS = {
    "pm": 1e-12, "nm": 1e-9, "um": 1e-6, "µm": 1e-6,
    "mm": 1e-3,  "m":  1.0,
}


def _parse_conc(val_str, col_hint=""):
    """Parse a concentration string ('0.43 nM', '1.6e-5', '430').
    Unit guessed from col_hint if not in string. Returns Molar or None.
    """
    val_str = str(val_str).strip()
    if not val_str or val_str in ("X", "NA", "nan", ""):
        return None
    m = _VAL_UNIT_PAT.match(val_str)
    if not m:
        return None
    num = float(m.group(1))
    if num <= 0:
        return None
    unit = m.group(2)
    if unit:
        factor = _UNIT_FACTORS.get(unit.lower(), 1.0)
    else:
        hint = col_hint.lower()
        if "nm" in hint:
            factor = 1e-9
        elif "um" in hint or "µm" in hint:
            factor = 1e-6
        elif "pm" in hint:
            factor = 1e-12
        elif "mm" in hint:
            factor = 1e-3
        else:
            factor = 1.0
    return num * factor


def _row_temperature(row, colmap, fallback_temp=DEFAULT_TEMP_K):
    """Return experimental temperature in K, or fallback_temp if not available."""
    if colmap.get("temperature"):
        raw = str(row.get(colmap["temperature"], "")).strip()
        try:
            t = float(raw)
            if 200.0 < t < 400.0:   # sanity: plausible experimental range in K
                return t
        except (ValueError, TypeError):
            pass
    return fallback_temp


def _compute_dg(row, colmap, fallback_temp=DEFAULT_TEMP_K):
    """
    Derive dG (kcal/mol), binding_type, dG_approx from a row.
    Priority: dg_exp → kd → ki → dg_calc.
    Uses experimental temperature if available, else fallback_temp.
    Returns (dg_str, btype, approx_flag) or ("NA", "NA", "N").
    """
    t = _row_temperature(row, colmap, fallback_temp)

    # 1. Experimental dG (direct kcal/mol — temperature already baked in)
    if colmap.get("dg_exp"):
        raw = str(row.get(colmap["dg_exp"], "")).strip()
        if raw and raw not in ("X", "NA", "nan", ""):
            try:
                return f"{float(raw):.3f}", "kd", "N"
            except ValueError:
                # Try "TYPE=VALUE UNIT" format (e.g. "Kd=1000 nM", "Ki=0.5 uM")
                m = re.match(r'(Kd|Ki|kd|ki)\s*=\s*(.+)', raw, re.IGNORECASE)
                if m:
                    atype = m.group(1).lower()
                    conc  = _parse_conc(m.group(2).strip(), atype)
                    if conc:
                        return f"{R_KCAL * t * math.log(conc):.3f}", atype, "N"

    # 2. Kd → dG = RT ln(Kd)
    if colmap.get("kd"):
        col  = colmap["kd"]
        conc = _parse_conc(row.get(col, ""), col)
        if conc:
            return f"{R_KCAL * t * math.log(conc):.3f}", "kd", "N"

    # 3. Ki → dG = RT ln(Ki)
    if colmap.get("ki"):
        col  = colmap["ki"]
        conc = _parse_conc(row.get(col, ""), col)
        if conc:
            return f"{R_KCAL * t * math.log(conc):.3f}", "ki", "N"

    # 4. Calculated dG (fallback — marked as approx)
    if colmap.get("dg_calc"):
        raw = str(row.get(colmap["dg_calc"], "")).strip()
        if raw and raw not in ("X", "NA", "nan", ""):
            try:
                return f"{float(raw):.3f}", "kd", "Y"
            except ValueError:
                pass

    return "NA", "NA", "N"


# ── Column name sanitization ─────────────────────────────────────────────────

_UNICODE_MAP = {
    'Δ': 'delta', 'δ': 'delta',
    'μ': 'u',     'µ': 'u',
    'α': 'alpha', 'β': 'beta', 'γ': 'gamma',
    'ε': 'eps',   'σ': 'sigma',
}


def _sanitize_col(name):
    """Replace common Unicode characters in column names with ASCII equivalents."""
    for uni, asc in _UNICODE_MAP.items():
        name = name.replace(uni, asc)
    return re.sub(r'[^\x00-\x7F]', '_', name).strip()


def _sanitize_val(val):
    """Transliterate common Unicode chars in values; replace rest with '?'."""
    val = str(val)
    for uni, asc in _UNICODE_MAP.items():
        val = val.replace(uni, asc)
    return re.sub(r'[^\x00-\x7F]', '?', val)


# ── Normalization ─────────────────────────────────────────────────────────────

def _get_n_mer(row, colmap):
    """Extract peptide length. Falls back to len(peptide_seq), then 0."""
    if colmap.get("n_mer"):
        val = str(row.get(colmap["n_mer"], "")).strip()
        try:
            return int(float(val))
        except (ValueError, TypeError):
            pass
    if colmap.get("peptide_seq"):
        seq = re.sub(r'[^A-Za-z]', '', str(row.get(colmap["peptide_seq"], "")))
        if seq:
            return len(seq)
    return 0


def normalize_entries(raw_rows, colmap, dataset_name, fallback_temp=DEFAULT_TEMP_K):
    """
    Convert raw external rows to peptide_entries-format dicts.
    Returns (entries, extra_col_names).
    extra_col_names = dataset columns not consumed by the standard pipeline.
    """
    mapped_cols = {v for v in colmap.values() if v}
    entries       = []
    skipped       = 0
    extra_col_set = set()

    for raw in raw_rows:
        pdb_col = colmap.get("pdb_id")
        pdb_id  = str(raw.get(pdb_col, "")).strip().lower() if pdb_col else ""
        if not re.match(r'^[a-z0-9]{4}$', pdb_id):
            skipped += 1
            continue

        dg_str, btype, dg_approx = _compute_dg(raw, colmap, fallback_temp)
        if dg_str == "NA":
            skipped += 1
            continue

        n_mer = _get_n_mer(raw, colmap)
        year  = (str(raw.get(colmap["year"],       "")).strip() or "NA") \
                if colmap.get("year")       else "NA"
        res   = (str(raw.get(colmap["resolution"], "")).strip() or "NA") \
                if colmap.get("resolution") else "NA"

        # Carry extra columns through with sanitized ASCII names and values
        extra = {_sanitize_col(k): _sanitize_val(v)
                 for k, v in raw.items() if k not in mapped_cols}
        extra_col_set.update(extra.keys())

        entry = {
            "pdb_id":       pdb_id,
            "n_mer":        str(n_mer),
            "year":         year,
            "binding_type": btype,
            "dG_kcal":      dg_str,
            "dG_approx":    dg_approx,
            "step1_reject": "",
            "resolution":   res,
            "ext_source":   dataset_name,
        }
        entry.update(extra)
        entries.append(entry)

    if skipped:
        print(f"  Skipped {skipped} rows (invalid PDB ID or no affinity data).")
    return entries, sorted(extra_col_set)


# ── Deduplication ─────────────────────────────────────────────────────────────

def load_existing_pdb_ids(outputs_dir):
    """Collect all pdb_ids already in PDBbind outputs and any previously processed external datasets."""
    existing = set()
    outputs_dir = Path(outputs_dir)
    skip = {"summary.tsv", "summary_merged.tsv"}
    sources = (
        [outputs_dir / f for f in ("peptide_entries.tsv", "summary.tsv")]
        + [p for p in outputs_dir.glob("summary_*.tsv") if p.name not in skip]
    )
    for path in sources:
        if not path.exists():
            continue
        for row in load_tsv(path):
            pid = row.get("pdb_id", "").strip().lower()
            if pid:
                existing.add(pid)
    return existing


# ── Merge ─────────────────────────────────────────────────────────────────────

def discover_external_summaries(outputs_dir):
    """Return sorted list of summary_<name>.tsv paths, excluding summary.tsv and summary_merged.tsv."""
    outputs_dir = Path(outputs_dir)
    skip = {"summary.tsv", "summary_merged.tsv"}
    found = sorted(
        p for p in outputs_dir.glob("summary_*.tsv")
        if p.name not in skip
    )
    return found


def merge_summaries(std_tsv, ext_tsvs, merged_tsv):
    """
    Concatenate standard PDBbind summary with one or more external summaries.
    Extra columns are appended at the right; PDBbind entries get empty strings for them.
    ext_tsvs: list of Path objects.
    """
    rows_std = load_tsv(std_tsv) if Path(std_tsv).exists() else []
    std_cols = list(rows_std[0].keys()) if rows_std else []

    all_cols = list(std_cols)
    ext_rows_by_source = []
    for tsv in ext_tsvs:
        rows = load_tsv(tsv) if tsv.exists() else []
        for c in (rows[0].keys() if rows else []):
            if c not in all_cols:
                all_cols.append(c)
        ext_rows_by_source.append((tsv.stem.removeprefix("summary_"), rows))

    merged = [{c: r.get(c, "") for c in all_cols} for r in rows_std]
    counts = []
    seen_ids: dict[str, str] = {}  # pdb_id -> source name
    for src_name, rows in ext_rows_by_source:
        for r in rows:
            pid = r.get("pdb_id", "").strip().lower()
            if pid and pid in seen_ids:
                print(f"  WARNING: {pid} appears in both '{seen_ids[pid]}' and '{src_name}' "
                      f"-- remove one summary_*.tsv and re-run --merge-only")
            elif pid:
                seen_ids[pid] = src_name
        merged.extend({c: r.get(c, "") for c in all_cols} for r in rows)
        counts.append(f"{src_name}={len(rows)}")

    save_tsv(merged, Path(merged_tsv), all_cols)
    ext_total = len(merged) - len(rows_std)
    print(f"  {len(rows_std)} PDBbind + {ext_total} external ({', '.join(counts)}) = {len(merged)} total")
    print(f"  Written {merged_tsv}")


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Integrate an external dataset into the PepXPro pipeline.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument("--input", required=False, metavar="FILE",
                        help="External dataset file (xlsx, tsv, csv). "
                             "Not required when --merge-only is used.")
    parser.add_argument("--name", metavar="NAME",
                        help="Short identifier for this dataset "
                             "(default: input filename stem).")
    parser.add_argument("--no-download", action="store_true",
                        help="Skip CIF download (assume CIFs already in pdb_cache/).")
    parser.add_argument("--no-filter", action="store_true",
                        help="Skip 02_filter.py step.")
    parser.add_argument("--no-summary", action="store_true",
                        help="Skip 03_summary.py step.")
    parser.add_argument("--no-merge", action="store_true",
                        help="Skip merge with PDBbind summary.tsv.")
    parser.add_argument("--merge-only", action="store_true",
                        help="Skip all processing steps; only re-run the merge of all "
                             "existing summary_*.tsv files into summary_merged.tsv.")
    parser.add_argument("--force-colmap", action="store_true",
                        help="Ignore saved column map and re-detect interactively.")
    parser.add_argument("--temperature", default=DEFAULT_TEMP_K, type=float, metavar="K",
                        help="Fallback temperature in K when not found in the dataset (default: 300.0).")
    parser.add_argument("--no-convert-mse", dest="convert_mse",
                        action="store_false", default=True,
                        help="Disable MSE→MET conversion in the filter step (on by default).")
    parser.add_argument("--strip-caps", action="store_true",
                        help="Remove ACE/NME/NH2 terminal caps if not at interface (>4.5 Å).")
    parser.add_argument("-v", "--verbose", action="store_true",
                        help="Print a brief start/done line to stdout.")
    parser.add_argument("--log", metavar="FILE", default=str(LOG_FILE),
                        help=f"Log file path (default: {LOG_FILE.name}).")
    args = parser.parse_args()

    log_file = Path(args.log)
    _logger, stdout = setup_logging(log_file, args.verbose)

    # ── Merge-only shortcut ───────────────────────────────────────────────────
    if args.merge_only:
        if args.verbose:
            stdout.write("[scrape/05-external] merge-only ...\n")
            stdout.flush()
        std_summary = OUTPUTS_DIR / "summary.tsv"
        summary_merged = OUTPUTS_DIR / "summary_merged.tsv"
        ext_tsvs = discover_external_summaries(OUTPUTS_DIR)
        if not ext_tsvs:
            print("No external summary_*.tsv files found in outputs/.")
            sys.exit(0)
        print(f"[Merge] Found {len(ext_tsvs)} external summary file(s):")
        for p in ext_tsvs:
            print(f"  {p.name}")
        merge_summaries(std_summary, ext_tsvs, summary_merged)
        print(f"  Use in genSample: python 01_genSample.py --summary-tsv {summary_merged}")
        if args.verbose:
            stdout.write(f"[scrape/05-external] done  —  merged\n")
            stdout.flush()
        sys.exit(0)

    if not hasattr(args, 'input') or not args.input:
        sys.stderr.write("Error: --input is required unless --merge-only is used.\n")
        sys.exit(1)

    input_path = Path(args.input)
    if not input_path.exists():
        input_path = _root / args.input  # try relative to scrape/
    if not input_path.exists():
        sys.stderr.write(f"Error: {Path(args.input)} not found.\n")
        sys.exit(1)

    name        = args.name or re.sub(r'[^a-z0-9_]', '_', input_path.stem.lower())
    colmap_path = _root / "inputs" / "external" / f"{name}_colmap.json"
    if args.force_colmap and colmap_path.exists():
        colmap_path.unlink()

    OUTPUTS_DIR.mkdir(parents=True, exist_ok=True)

    if args.verbose:
        stdout.write(f"[scrape/05-external] starting ... ({name})\n")
        stdout.flush()

    entries_ext_tsv  = OUTPUTS_DIR / f"peptide_entries_{name}.tsv"
    filtered_ext_tsv = OUTPUTS_DIR / f"filtered_systems_{name}.tsv"
    summary_ext_tsv  = OUTPUTS_DIR / f"summary_{name}.tsv"
    summary_merged   = OUTPUTS_DIR / "summary_merged.tsv"

    # ── 1. Load file ──────────────────────────────────────────────────────────
    print(f"\n[1/6] Loading {input_path} ...")
    raw_rows = _load_file(input_path)
    print(f"  {len(raw_rows)} rows")
    if not raw_rows:
        print("  Empty file. Aborting.")
        sys.exit(1)

    # ── 2. Column map ─────────────────────────────────────────────────────────
    print(f"\n[2/6] Resolving column map ...")
    columns = list(raw_rows[0].keys())
    colmap  = _resolve_colmap(columns, name, colmap_path)

    # ── 3. Normalize ──────────────────────────────────────────────────────────
    print(f"\n[3/6] Normalizing entries ...")
    entries, extra_col_names = normalize_entries(raw_rows, colmap, name, fallback_temp=args.temperature)
    print(f"  {len(entries)} valid entries")
    if not entries:
        print("  No valid entries. Aborting.")
        sys.exit(1)

    # ── 4. Deduplicate ────────────────────────────────────────────────────────
    print(f"\n[4/6] Deduplicating against PDBbind ...")
    existing = load_existing_pdb_ids(OUTPUTS_DIR)
    before   = len(entries)
    entries  = [e for e in entries if e["pdb_id"] not in existing]
    print(f"  {before - len(entries)} duplicates removed -> {len(entries)} new entries")
    if not entries:
        print("  All entries already in PDBbind. Nothing to do.")
        sys.exit(0)

    # ── 5. Download CIFs ──────────────────────────────────────────────────────
    if not args.no_download:
        print(f"\n[5/6] Downloading CIFs ...")
        pdb_ids = [e["pdb_id"] for e in entries]
        counts, errors = download_entries(pdb_ids, PDB_CACHE)
        print(f"  ok={counts['ok']}  cached={counts['cached']}  errors={counts['error']}")
        if errors:
            for pid, msg in errors[:10]:
                print(f"    {pid}: {msg}")
    else:
        print(f"\n[5/6] Skipping CIF download.")

    # Fill year and resolution from CIF now that CIFs are downloaded
    from modules.metadata import (extract_deposition_year, extract_resolution,
                                   fetch_resolution_rcsb)
    for field, extractor, label in (
        ("year",       extract_deposition_year, "year"),
        ("resolution", extract_resolution,      "resolution"),
    ):
        missing = [e for e in entries if e[field] == "NA"]
        if not missing:
            continue
        print(f"  Filling {label} from CIF for {len(missing)} entries ...")
        filled = 0
        for e in missing:
            cif_path = PDB_CACHE / f"{e['pdb_id']}-1.cif"
            if cif_path.exists():
                val = extractor(cif_path)
                if val != "NA":
                    e[field] = val
                    filled += 1
        print(f"  {label.capitalize()} filled for {filled}/{len(missing)} entries "
              f"({len(missing) - filled} remain NA).")

    # Resolution fallback: RCSB REST API (assembly CIFs lack _refine.ls_d_res_high)
    still_missing = [e for e in entries if e["resolution"] == "NA"]
    if still_missing:
        from modules.metadata import extract_exptl_method
        print(f"  Filling resolution from RCSB API for {len(still_missing)} entries ...")
        filled = 0
        for e in still_missing:
            val = fetch_resolution_rcsb(e["pdb_id"])
            if val != "NA":
                e["resolution"] = val
                filled += 1
        remaining = [e for e in still_missing if e["resolution"] == "NA"]
        print(f"  Resolution filled for {filled}/{len(still_missing)} entries "
              f"({len(remaining)} remain NA).")
        if remaining:
            method_counts: dict[str, list[str]] = {}
            for e in remaining:
                cif = PDB_CACHE / f"{e['pdb_id']}-1.cif"
                method = extract_exptl_method(cif) if cif.exists() else "NA"
                method_counts.setdefault(method, []).append(e["pdb_id"])
            for method, pids in sorted(method_counts.items()):
                print(f"    {method}: {len(pids)} entries (e.g. {', '.join(pids[:5])})")

    # Write peptide_entries_<name>.tsv (after CIF fill)
    entry_cols = ["pdb_id", "n_mer", "year", "binding_type",
                  "dG_kcal", "dG_approx", "step1_reject", "resolution",
                  "ext_source"] + extra_col_names
    save_tsv(entries, entries_ext_tsv, entry_cols)
    print(f"  Written {entries_ext_tsv}")

    # ── 6a. Filter ────────────────────────────────────────────────────────────
    if not args.no_filter:
        print(f"[6a/6] Running 02_filter.py ...")
        cmd = [sys.executable, str(_here / "02_filter.py"),
               "--input",  str(entries_ext_tsv),
               "--output", str(filtered_ext_tsv),
               "--log",    str(LOGS_DIR / f"02_filter_{name}.log")]
        if args.verbose:
            cmd.append("-v")
        if not args.convert_mse:
            cmd.append("--no-convert-mse")
        if args.strip_caps:
            cmd.append("--strip-caps")
        result = subprocess.run(cmd, check=False, stderr=subprocess.PIPE)
        if result.stderr:
            _logger.warning(result.stderr.decode(errors="replace").strip())
        if result.returncode != 0:
            sys.stderr.write(f"02_filter.py failed (exit {result.returncode})\n")
            sys.exit(result.returncode)
    else:
        print("[6a/6] Skipping filter step.")

    # ── 6b. Summary ───────────────────────────────────────────────────────────
    if not args.no_summary and filtered_ext_tsv.exists():
        print(f"[6b/6] Running 03_summary.py ...")
        cmd = [sys.executable, str(_here / "03_summary.py"),
               "--entries",  str(entries_ext_tsv),
               "--filtered", str(filtered_ext_tsv),
               "--output",   str(summary_ext_tsv),
               "--log",      str(LOGS_DIR / f"03_summary_{name}.log")]
        if args.verbose:
            cmd.append("-v")
        result = subprocess.run(cmd, check=False, stderr=subprocess.PIPE)
        if result.stderr:
            _logger.warning(result.stderr.decode(errors="replace").strip())
        if result.returncode != 0:
            sys.stderr.write(f"03_summary.py failed (exit {result.returncode})\n")
            sys.exit(result.returncode)
    else:
        print("[6b/6] Skipping summary step.")

    # ── Merge ─────────────────────────────────────────────────────────────────
    if not args.no_merge and summary_ext_tsv.exists():
        std_summary = OUTPUTS_DIR / "summary.tsv"
        ext_tsvs = discover_external_summaries(OUTPUTS_DIR)
        print(f"[Merge] Merging {len(ext_tsvs)} external source(s) with PDBbind summary ...")
        if not std_summary.exists():
            print(f"  Warning: {std_summary} not found — merged will contain external only.")
        merge_summaries(std_summary, ext_tsvs, summary_merged)
    else:
        print("[Merge] Skipping merge.")

    print(f"Done. External summary: {summary_ext_tsv}")
    if not args.no_merge:
        print(f"  Merged: {summary_merged}")

    if args.verbose:
        stdout.write(f"[scrape/05-external] done  —  {name}\n")
        stdout.flush()


if __name__ == "__main__":
    main()
