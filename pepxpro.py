"""
PepXPro — peptide-protein benchmark pipeline.

Usage:
  python pepxpro.py scrape -v
  python pepxpro.py scrape --steps filter summary
  python pepxpro.py scrape --force
  python pepxpro.py scrape -v --external-input inputs/external/PEPBI.xlsx
"""

import argparse
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent

SCRAPE_SCRIPTS = {
    "download":      ROOT / "scrape_pdb_bind/scripts/01_download.py",
    "filter":        ROOT / "scrape_pdb_bind/scripts/02_filter.py",
    "summary":       ROOT / "scrape_pdb_bind/scripts/03_summary.py",
    "external":      ROOT / "scrape_pdb_bind/scripts/04_external.py",
    "distributions": ROOT / "scrape_pdb_bind/scripts/05_distributions.py",
}
SCRAPE_ORDER         = ["download", "filter", "summary", "distributions"]
SCRAPE_ORDER_EXT     = ["download", "filter", "summary", "external", "distributions"]
SCRAPE_STEP_CHOICES  = ["download", "filter", "summary", "external", "distributions"]


def run_scrape(args):
    if args.steps:
        steps = args.steps
    elif getattr(args, "external_input", None):
        steps = SCRAPE_ORDER_EXT
    else:
        steps = SCRAPE_ORDER

    for step in steps:
        cmd = [sys.executable, str(SCRAPE_SCRIPTS[step])]
        if args.verbose:
            cmd.append("-v")
        if step == "filter" and getattr(args, "force", False):
            cmd.append("--force")
        if step == "external":
            if not getattr(args, "external_input", None):
                sys.stderr.write("[scrape/external] --external-input FILE is required.\n")
                sys.exit(1)
            cmd.extend(["--input", str(args.external_input)])
        if step == "distributions" and "external" in steps:
            merged = ROOT / "scrape_pdb_bind/outputs/summary_merged.tsv"
            cmd.extend(["--summary-tsv", str(merged)])

        result = subprocess.run(
            cmd,
            stdout=None if args.verbose else subprocess.DEVNULL,
            stderr=subprocess.PIPE,
        )

        if result.returncode != 0:
            sys.stderr.write(f"\n[scrape/{step}] ERROR (exit {result.returncode})\n")
            if result.stderr:
                sys.stderr.write(result.stderr.decode())
            sys.exit(result.returncode)


def main():
    parser = argparse.ArgumentParser(
        prog="pepxpro",
        description="PepXPro — peptide-protein benchmark pipeline.",
    )
    sub = parser.add_subparsers(dest="command", required=True)

    p_scrape = sub.add_parser("scrape", help="Run the scrape_pdb_bind pipeline.")
    p_scrape.add_argument("-v", "--verbose", action="store_true",
                          help="Print one brief line per step to stdout.")
    p_scrape.add_argument(
        "--steps", nargs="+", choices=SCRAPE_STEP_CHOICES, metavar="STEP",
        help=f"Steps to run (default: all standard). Choices: {', '.join(SCRAPE_STEP_CHOICES)}.",
    )
    p_scrape.add_argument("--force", action="store_true",
                          help="Force reprocess in filter step (ignore cache).")
    p_scrape.add_argument("--external-input", metavar="FILE",
                          help="External dataset file (xlsx/tsv/csv). "
                               "Activates the external step (04_external.py) automatically.")

    args = parser.parse_args()
    if args.command == "scrape":
        run_scrape(args)


if __name__ == "__main__":
    main()
