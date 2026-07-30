"""
PepXPro — peptide-protein benchmark pipeline.

Usage:
  python pepxpro.py scrape -v
  python pepxpro.py scrape --steps filter summary
  python pepxpro.py scrape --force
  python pepxpro.py scrape -v --external-input inputs/external/PEPBI.xlsx
  python pepxpro.py scrape -v --index scrape/inputs/INDEX_general_PL.2020
  python pepxpro.py config --print
"""

import argparse
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent

SCRAPE_SCRIPTS = {
    "download":      ROOT / "scrape/scripts/01_download.py",
    "filter":        ROOT / "scrape/scripts/02_filter.py",
    "summary":       ROOT / "scrape/scripts/03_summary.py",
    "external":      ROOT / "scrape/scripts/04_external.py",
    "distributions": ROOT / "scrape/scripts/05_distributions.py",
    "tm_matrix":     ROOT / "scrape/scripts/06_tm_matrix.py",
}
SCRAPE_ORDER         = ["download", "filter", "summary", "distributions"]
SCRAPE_ORDER_EXT     = ["download", "filter", "summary", "external", "distributions"]
SCRAPE_STEP_CHOICES  = ["download", "filter", "summary", "external", "distributions", "tm_matrix"]


def run_scrape(args):
    if args.steps:
        steps = args.steps
    elif getattr(args, "external_input", None):
        steps = SCRAPE_ORDER_EXT + (["tm_matrix"] if getattr(args, "precompute_tm", False) else [])
    else:
        steps = SCRAPE_ORDER + (["tm_matrix"] if getattr(args, "precompute_tm", False) else [])

    for step in steps:
        cmd = [sys.executable, str(SCRAPE_SCRIPTS[step])]
        if args.verbose:
            cmd.append("-v")
        if step == "download" and getattr(args, "index", None):
            cmd.extend(["--index", str(args.index)])
        if step == "filter" and getattr(args, "force", False):
            cmd.append("--force")
        if step in ("filter", "external") and not getattr(args, "convert_mse", True):
            cmd.append("--no-convert-mse")
        if step in ("filter", "external") and getattr(args, "strip_caps", False):
            cmd.append("--strip-caps")
        if step == "external":
            if not getattr(args, "external_input", None):
                sys.stderr.write("[scrape/external] --external-input FILE is required.\n")
                sys.exit(1)
            cmd.extend(["--input", str(args.external_input)])
            if getattr(args, "name", None):
                cmd.extend(["--name", args.name])
        if step == "distributions":
            summary_dir = ROOT / "scrape/outputs"
            if getattr(args, "summary_tsv", None):
                cmd.extend(["--summary-tsv", str(args.summary_tsv)])
            elif "external" in steps:
                cmd.extend(["--summary-tsv", str(summary_dir / "summary_merged.tsv")])
            else:
                cmd.extend(["--summary-tsv", str(summary_dir / "summary.tsv")])

        t0 = time.perf_counter()
        result = subprocess.run(
            cmd,
            stdout=None if args.verbose else subprocess.DEVNULL,
            stderr=subprocess.PIPE,
        )
        elapsed = time.perf_counter() - t0

        if result.returncode != 0:
            sys.stderr.write(f"\n[scrape/{step}] ERROR (exit {result.returncode})\n")
            if result.stderr:
                sys.stderr.write(result.stderr.decode())
            sys.exit(result.returncode)

        if args.verbose:
            mins, secs = divmod(int(elapsed), 60)
            t_str = f"{mins}m {secs}s" if mins else f"{secs}s"
            sys.stdout.write(f"[scrape/{step}] elapsed: {t_str}\n")
            sys.stdout.flush()


def run_config(_args):
    import re

    def _const(filepath, varname):
        """Read one constant from a Python source file by regex, stripping inline comments."""
        src = Path(filepath).read_text(encoding="utf-8")
        m = re.search(rf'^{varname}\s*=\s*(.+)', src, re.MULTILINE)
        if not m:
            return "?"
        val = m.group(1)
        val = re.sub(r'\s*#.*$', '', val)  # strip inline comment
        return val.strip()

    FILTER_PY = ROOT / "scrape/modules/filter.py"
    DL_PY     = ROOT / "scrape/modules/download.py"
    SR_PY     = ROOT / "genSample/modules/sequence_redundancy.py"
    BENCH_PY  = ROOT / "pepxpro_benchmark/scripts/02_make_benchmark.py"

    W = 48  # column width for alignment

    def row(key, val, note="", flag=""):
        label = f"[{flag}]" if flag else ""
        line  = f"{key:<28} = {str(val):<14} {label}"
        if note:
            line += f"  # {note}"
        print(line)

    def section(title):
        bar = "─" * (W - len(title) - 3)
        print(f"\n# ── {title} {bar}")

    print("# PepXPro — default parameters")
    print("# Usage: python pepxpro.py config --print")
    print("#")
    print("# [CLI]        change via command-line flag (see --help)")
    print("# [hardcoded]  edit the source file directly to change")

    # ── scrape / download ─────────────────────────────────────────────────────
    section("scrape / 01_download.py")
    row("index",         "inputs/INDEX_general_PL.2020R1.lst", "PDBbind index file",     "CLI")
    row("temperature_K", 300.0,  "fallback when not in index comment field",             "CLI")
    row("n_workers",     _const(DL_PY, "MAX_WORKERS"), "concurrent CIF downloads",        "hardcoded")

    # ── scrape / filter ───────────────────────────────────────────────────────
    section("scrape / 02_filter.py")
    row("convert_mse",            "true",  "MSE→MET via PDBFixer  (--no-convert-mse)",  "CLI")
    row("strip_caps",             "false", "remove ACE/NME/NH2 if not at interface",    "CLI")
    row("neighbor_cutoff_A",      _const(FILTER_PY, "NEIGHBOR_CUTOFF_A"),      "receptor identification radius","hardcoded")
    row("max_receptor_gap",       _const(FILTER_PY, "MAX_RECEPTOR_GAP"),       "max gap in receptor numbering", "hardcoded")
    row("peptide_bond_cutoff_A",  _const(FILTER_PY, "PEPTIDE_BOND_CUTOFF"),    "covalent terminus detection",   "hardcoded")
    row("peptide_copy_threshold", _const(FILTER_PY, "PEPTIDE_COPY_THRESHOLD"), "seq similarity → reject peptide copy", "hardcoded")
    row("contact_cutoff_A",       _const(FILTER_PY, "CONTACT_CUTOFF_A"),       "n_contacts heavy-atom cutoff",  "hardcoded")
    row("hbond_cutoff_A",         _const(FILTER_PY, "HBOND_CUTOFF_A"),         "n_hbonds N/O/S proximity proxy","hardcoded")
    row("cap_interface_cutoff_A", _const(FILTER_PY, "CAP_INTERFACE_CUTOFF_A"), "reject cap if closer than this to receptor", "hardcoded")

    # ── genSample / quality ──────────────────────────────────────────────────
    section("genSample / quality filters  (01_genSample.py)")
    row("single_chain",    "true",  "receptor must have exactly 1 chain (--no-single-chain)", "CLI")
    row("min_peptide",     "",      "--min-peptide N",           "CLI")
    row("max_peptide",     "",      "--max-peptide N",           "CLI")
    row("min_receptor",    "",      "--min-receptor N",          "CLI")
    row("max_receptor",    "",      "--max-receptor N",          "CLI")
    row("min_dg",          "",      "kcal/mol  --min-dg F",      "CLI")
    row("max_dg",          "",      "kcal/mol  --max-dg F",      "CLI")
    row("min_resolution",  "",      "Å  --min-resolution F",     "CLI")
    row("max_resolution",  "",      "Å  --max-resolution F",     "CLI")
    row("binding_type",    "both",  "kd | ki | both",            "CLI")
    row("max_added",       "",      "max residues added by PDBFixer to peptide", "CLI")
    row("max_receptor_added", "",   "max residues added by PDBFixer to receptor","CLI")

    # ── genSample / redundancy ───────────────────────────────────────────────
    section("genSample / redundancy filters  (01_genSample.py)")
    row("seq_id_receptor",  "strict", "preset → CD-HIT ≥ 40%   (strict/medium/lenient/custom)", "CLI")
    row("seq_id_peptide",   "strict", "preset → pairwise ≥ 70%  (strict/medium/lenient/custom)", "CLI")
    row("tm_score_receptor","",       "TM-score threshold  --tm-score-receptor F  (disabled by default)", "CLI")

    # ── benchmark ─────────────────────────────────────────────────────────────
    section("benchmark / 02_make_benchmark.py  (frozen profiles, edit file directly)")
    if BENCH_PY.exists():
        import importlib.util
        spec     = importlib.util.spec_from_file_location("_bench_config_mod", BENCH_PY)
        bench_mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(bench_mod)
        for name, profile in bench_mod.PROFILES.items():
            row(f"PROFILES['{name}']", f"-> {profile['version']}", "frozen release", "hardcoded")
            row("  seq_id_receptor",  profile["seq_id_receptor"], flag="hardcoded")
            row("  seq_id_peptide",   profile["seq_id_peptide"],  flag="hardcoded")
            row("  tm_score",         profile["tm_score"],        flag="hardcoded")
            row("  sampling_method",  profile["sampling_method"], flag="hardcoded")
            row("  dims",             "+".join(profile["dims"]),  flag="hardcoded")
            row("  sample_n",         profile["sample_n"],        flag="hardcoded")
            row("  seed",             profile["seed"],            flag="hardcoded")
        print("  Also available: --from-config auto|GCn  (derives params live from "
              "pepxpro_benchmark/optimization/results.tsv instead of a frozen profile)")
    else:
        print(f"  (script not found: {BENCH_PY})")

    print()


def main():
    parser = argparse.ArgumentParser(
        prog="pepxpro",
        description="PepXPro — peptide-protein benchmark pipeline.",
    )
    sub = parser.add_subparsers(dest="command", required=True)

    p_scrape = sub.add_parser("scrape", help="Run the scrape pipeline.")
    p_scrape.add_argument("-v", "--verbose", action="store_true",
                          help="Print one brief line per step to stdout.")
    p_scrape.add_argument(
        "--steps", nargs="+", choices=SCRAPE_STEP_CHOICES, metavar="STEP",
        help=f"Steps to run (default: all standard). Choices: {', '.join(SCRAPE_STEP_CHOICES)}.",
    )
    p_scrape.add_argument("--index", metavar="FILE",
                          help="PDBbind index file (default: inputs/INDEX_general_PL.2020R1.lst). "
                               "Passed to 01_download.py.")
    p_scrape.add_argument("--force", action="store_true",
                          help="Force reprocess in filter step (ignore cache).")
    p_scrape.add_argument("--no-convert-mse", dest="convert_mse",
                          action="store_false", default=True,
                          help="Disable MSE→MET conversion (enabled by default).")
    p_scrape.add_argument("--strip-caps", action="store_true",
                          help="Strip ACE/NME/NH2 terminal caps if not at interface (>4.5 Å).")
    p_scrape.add_argument("--precompute-tm", action="store_true",
                          help="Pre-compute receptor TM-score matrix after filtering "
                               "(06_tm_matrix.py). Saves outputs/tm_score_matrix.npz "
                               "for fast reuse in genSample --tm-score-receptor.")
    p_scrape.add_argument("--external-input", metavar="FILE",
                          help="External dataset file (xlsx/tsv/csv). "
                               "Activates the external step (04_external.py) automatically.")
    p_scrape.add_argument("--name", metavar="NAME",
                          help="Short identifier for external dataset output files "
                               "(e.g. --name PEPBI → summary_PEPBI.tsv). "
                               "Defaults to the input filename stem.")
    p_scrape.add_argument("--summary-tsv", metavar="FILE", dest="summary_tsv",
                          help="Summary TSV to use for the distributions step "
                               "(e.g. summary_merged.tsv). Inferred automatically when omitted.")

    p_config = sub.add_parser("config", help="Show pipeline default parameters.")
    p_config.add_argument("--print", dest="print_params", action="store_true",
                          help="Print all default parameters for every pipeline stage.")

    args = parser.parse_args()
    if args.command == "scrape":
        run_scrape(args)
    elif args.command == "config":
        if args.print_params:
            run_config(args)
        else:
            p_config.print_help()


if __name__ == "__main__":
    main()
