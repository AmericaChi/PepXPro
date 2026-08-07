"""
PepXPro — peptide-protein benchmark pipeline.

Usage:
  python pepxpro.py scrape -v
  python pepxpro.py scrape --steps filter summary
  python pepxpro.py scrape --force
  python pepxpro.py scrape -v --external-input inputs/external/PEPBI.xlsx
  python pepxpro.py scrape -v --index scrape/inputs/INDEX_general_PL.2020
  python pepxpro.py genSample gen --sample 20 --sampling kennard-stone --seed 42
  python pepxpro.py genSample plot --run run_01
  python pepxpro.py benchmark optimize --full
  python pepxpro.py benchmark make --from-config auto --version v1
  python pepxpro.py config --print

Every module/step above is also runnable directly, e.g.:
  python genSample/scripts/01_gen_sample.py --sample 20 ...
  python pepxpro_benchmark/scripts/02_make_benchmark.py --from-config auto ...

After `pip install -e .` (see README), the same commands work without the
"python " / ".py" via the bare `pepxpro` entry point, e.g. `pepxpro scrape -v`.
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

GENSAMPLE_SCRIPTS = {
    "methods": ROOT / "genSample/scripts/00_plot_samplingMethods.py",
    "gen":     ROOT / "genSample/scripts/01_gen_sample.py",
    "plot":    ROOT / "genSample/scripts/02_plot_sampleDist.py",
    "view":    ROOT / "genSample/scripts/03_view_sample.py",
}

BENCHMARK_SCRIPTS = {
    "optimize": ROOT / "pepxpro_benchmark/scripts/01_optimize_benchmark.py",
    "make":     ROOT / "pepxpro_benchmark/scripts/02_make_benchmark.py",
}


def _run_step(cmd, label, verbose):
    """Run one subprocess step, printing an error on failure and (if verbose)
    a one-line elapsed-time summary on success."""
    t0 = time.perf_counter()
    result = subprocess.run(
        cmd,
        stdout=None if verbose else subprocess.DEVNULL,
        stderr=subprocess.PIPE,
    )
    elapsed = time.perf_counter() - t0

    if result.returncode != 0:
        sys.stderr.write(f"\n[{label}] ERROR (exit {result.returncode})\n")
        if result.stderr:
            sys.stderr.write(result.stderr.decode())
        sys.exit(result.returncode)

    if verbose:
        mins, secs = divmod(int(elapsed), 60)
        t_str = f"{mins}m {secs}s" if mins else f"{secs}s"
        sys.stdout.write(f"[{label}] elapsed: {t_str}\n")
        sys.stdout.flush()


def _run_passthrough(cmd):
    """Run one subprocess step with stdout/stderr streamed live (no -v gate —
    used for genSample/benchmark scripts, which print their own progress
    unconditionally and have no quiet mode)."""
    result = subprocess.run(cmd)
    if result.returncode != 0:
        sys.exit(result.returncode)


def _build_argv(args, spec):
    """Reconstruct a subprocess argv fragment from a parsed argparse Namespace
    and a declarative flag spec: a list of (dest, kind) tuples. The flag name
    is derived from dest (e.g. "seq_id_receptor" -> "--seq-id-receptor"),
    which always matches argparse's own auto-derived dest for a flag added
    without an explicit dest= override.

    Every mirrored flag in this file is declared with default=None, so only
    flags the user actually typed are forwarded — the real script underneath
    always supplies its own defaults for anything omitted here. This means
    pepxpro.py's copy of each flag can never silently drift from the real
    script's default behavior, only from its flag *names*.

    kind:
      "store_true" -> append flag if True
      "bool_opt"   -> BooleanOptionalAction; append flag or --no-<flag> based on value
      "value"      -> append flag, str(value)
      "nargs"      -> append flag, *values
      "append"     -> append (flag, item) for each item in the list
    """
    argv = []
    for dest, kind in spec:
        val = getattr(args, dest, None)
        if val is None:
            continue
        flag = "--" + dest.replace("_", "-")
        if kind == "store_true":
            if val:
                argv.append(flag)
        elif kind == "bool_opt":
            argv.append(flag if val else flag.replace("--", "--no-", 1))
        elif kind == "value":
            argv += [flag, str(val)]
        elif kind == "nargs":
            argv += [flag, *[str(v) for v in val]]
        elif kind == "append":
            for v in val:
                argv += [flag, str(v)]
        else:
            raise ValueError(f"unknown _build_argv kind: {kind}")
    return argv


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

        _run_step(cmd, f"scrape/{step}", args.verbose)


# ---------------------------------------------------------------------------
# genSample
# ---------------------------------------------------------------------------

def _add_gensample_gen_parser(gs_sub):
    p = gs_sub.add_parser(
        "gen", help="Filter, deduplicate, and sample a benchmark set (01_gen_sample.py).")

    mode = p.add_mutually_exclusive_group(required=True)
    mode.add_argument("--sample", metavar="N|PCT%", default=None,
                      help="Sample N entries, or a percentage of the filtered pool (e.g. --sample 20 or --sample 30%%).")
    mode.add_argument("--no-sample", action="store_true", default=None,
                      help="Take the full filtered pool without sampling.")
    mode.add_argument("--dry-run", action="store_true", default=None,
                      help="Print pool statistics after filtering and exit (no output written).")
    mode.add_argument("--pick", metavar="PDB_ID", default=None,
                      help="Extract a single system by PDB ID, bypassing all filters and sampling.")

    samp = p.add_argument_group("sampling (only with --sample)")
    samp.add_argument("--sampling", choices=["stratified", "kennard-stone", "kde", "isolation-forest", "cluster"],
                      metavar="METHOD", default=None, help="Algorithm (default: stratified).")
    samp.add_argument("--dims", nargs="+", metavar="DIM", default=None,
                      help="Feature dimensions (default: dg pep_len). See genSample/CLAUDE.md for the full list.")
    samp.add_argument("--seed", type=int, metavar="N", default=None, help="Random seed for reproducibility.")

    filt = p.add_argument_group("filters (applied before sampling; all optional)")
    filt.add_argument("--single-chain", action=argparse.BooleanOptionalAction, default=None,
                      help="Keep only single receptor-chain entries (script default: on).")
    filt.add_argument("--min-receptor", type=int, metavar="N", default=None)
    filt.add_argument("--max-receptor", type=int, metavar="N", default=None)
    filt.add_argument("--min-peptide", type=int, metavar="N", default=None)
    filt.add_argument("--max-peptide", type=int, metavar="N", default=None)
    filt.add_argument("--min-dg", type=float, metavar="F", default=None)
    filt.add_argument("--max-dg", type=float, metavar="F", default=None)
    filt.add_argument("--binding-type", choices=["kd", "ki", "both"], default=None,
                      help="Measurement type (script default: both).")
    filt.add_argument("--max-added", type=int, metavar="N", default=None,
                      help="Max residues PDBFixer added to peptide.")
    filt.add_argument("--max-receptor-added", type=int, metavar="N", default=None,
                      help="Max residues PDBFixer added to receptor.")
    filt.add_argument("--max-resolution", type=float, metavar="F", default=None)
    filt.add_argument("--min-resolution", type=float, metavar="F", default=None)
    filt.add_argument("--exclude-nmr", action="store_true", default=None,
                      help="Exclude entries with no numeric resolution (NMR, NA).")
    filt.add_argument("--min-year", type=int, metavar="N", default=None)
    filt.add_argument("--max-year", type=int, metavar="N", default=None)
    filt.add_argument("--exclude", metavar="FILE", default=None, help="File with one pdb_id per line to exclude.")
    filt.add_argument("--peptide-topology", choices=["linear", "cyclic_disulfide", "both"], default=None,
                      help="Filter by peptide topology (script default: both).")

    red = p.add_argument_group("redundancy (applied after quality filters, before sampling)")
    red.add_argument("--seq-id-receptor", choices=["strict", "medium", "lenient", "custom"],
                     metavar="PRESET", default=None, help="CD-HIT receptor redundancy preset.")
    red.add_argument("--seq-id-receptor-threshold", type=float, metavar="F", default=None,
                     help="Receptor identity threshold for --seq-id-receptor custom.")
    red.add_argument("--seq-id-peptide", choices=["strict", "medium", "lenient", "custom"],
                     metavar="PRESET", default=None, help="Biopython pairwise peptide redundancy preset.")
    red.add_argument("--seq-id-peptide-threshold", type=float, metavar="F", default=None,
                     help="Peptide identity threshold for --seq-id-peptide custom.")

    struct = p.add_argument_group("structural redundancy (requires fixed_cif/)")
    struct.add_argument("--tm-score-receptor", type=float, metavar="F", default=None,
                        help="Receptor TM-score threshold (e.g. 0.5 = same fold). Requires tmtools.")

    inp = p.add_argument_group("input")
    inp.add_argument("--summary-tsv", metavar="FILE", default=None,
                     help="Path to summary TSV (script default: scrape/outputs/summary.tsv).")

    out = p.add_argument_group("output")
    out.add_argument("--no-save-pools", action="store_true", default=None,
                     help="Do not save pool_quality.tsv / pool_redundancy.tsv.")
    out.add_argument("--full-metadata", action="store_true", default=None,
                     help="Write all columns from summary.tsv to sample.tsv.")

    exp = p.add_argument_group("export (all optional; can be combined)")
    exp.add_argument("--link-cifs", action="store_true", default=None, help="Symlink structure files.")
    exp.add_argument("--copy-cifs", action="store_true", default=None, help="Copy structure files.")
    exp.add_argument("--protonate", action="store_true", default=None, help="Add missing hydrogens (implies copy).")
    exp.add_argument("--ph", type=float, metavar="F", default=None, help="pH for protonation (script default: 7.0).")
    exp.add_argument("--format", choices=["cif", "pdb"], default=None, help="Output format (script default: cif).")
    exp.add_argument("--export-fasta", action="store_true", default=None,
                     help="Write peptide/receptor FASTA files to the run directory.")


GENSAMPLE_GEN_SPEC = [
    ("sample", "value"), ("no_sample", "store_true"), ("dry_run", "store_true"), ("pick", "value"),
    ("sampling", "value"), ("dims", "nargs"), ("seed", "value"),
    ("single_chain", "bool_opt"),
    ("min_receptor", "value"), ("max_receptor", "value"),
    ("min_peptide", "value"), ("max_peptide", "value"),
    ("min_dg", "value"), ("max_dg", "value"), ("binding_type", "value"),
    ("max_added", "value"), ("max_receptor_added", "value"),
    ("max_resolution", "value"), ("min_resolution", "value"), ("exclude_nmr", "store_true"),
    ("min_year", "value"), ("max_year", "value"), ("exclude", "value"), ("peptide_topology", "value"),
    ("seq_id_receptor", "value"), ("seq_id_receptor_threshold", "value"),
    ("seq_id_peptide", "value"), ("seq_id_peptide_threshold", "value"),
    ("tm_score_receptor", "value"),
    ("summary_tsv", "value"),
    ("no_save_pools", "store_true"), ("full_metadata", "store_true"),
    ("link_cifs", "store_true"), ("copy_cifs", "store_true"), ("protonate", "store_true"),
    ("ph", "value"), ("format", "value"), ("export_fasta", "store_true"),
]


def _add_gensample_methods_parser(gs_sub):
    p = gs_sub.add_parser(
        "methods", help="Benchmark and compare sampling methods (00_plot_samplingMethods.py).")
    p.add_argument("--n", type=int, metavar="N", default=None, help="Entries to sample (script default: 20).")
    p.add_argument("--seed", type=int, metavar="N", default=None, help="Random seed for stochastic methods.")
    p.add_argument("--dims", nargs="+", metavar="DIM", default=None,
                  help="Two+ dimensions for the plot axes (script default: dg pep_len).")


GENSAMPLE_METHODS_SPEC = [("n", "value"), ("seed", "value"), ("dims", "nargs")]


def _add_gensample_plot_parser(gs_sub):
    p = gs_sub.add_parser(
        "plot", help="Per-run filtering/sampling distribution plot (02_plot_sampleDist.py).")
    p.add_argument("--run", metavar="RUN_DIR", default=None,
                  help="Run directory name inside genSample/outputs/, e.g. run_01.")
    p.add_argument("--run-dir", metavar="PATH", default=None,
                  help="Full path to any directory with pool_quality.tsv, pool_redundancy.tsv, sample.tsv.")
    p.add_argument("--dims", nargs=2, metavar=("DIM_X", "DIM_Y"), default=None,
                  help="Override sampling dims from run_config.")
    p.add_argument("--output", metavar="FILE", default=None, help="Override output path.")


GENSAMPLE_PLOT_SPEC = [("run", "value"), ("run_dir", "value"), ("dims", "nargs"), ("output", "value")]


def _add_gensample_view_parser(gs_sub):
    p = gs_sub.add_parser(
        "view", help="Browse a run's structures in PyMOL (03_view_sample.py).")
    p.add_argument("--run", type=int, metavar="N", default=None, help="Run number (default: latest run).")
    p.add_argument("--run-dir", metavar="PATH", default=None,
                  help="Full path to any directory with sample.tsv. Overrides --run.")
    p.add_argument("--cif-dir", metavar="DIR", default=None, help="Override CIF directory.")


GENSAMPLE_VIEW_SPEC = [("run", "value"), ("run_dir", "value"), ("cif_dir", "value")]

GENSAMPLE_SPECS = {
    "gen": GENSAMPLE_GEN_SPEC,
    "methods": GENSAMPLE_METHODS_SPEC,
    "plot": GENSAMPLE_PLOT_SPEC,
    "view": GENSAMPLE_VIEW_SPEC,
}


def run_gensample(args):
    cmd = [sys.executable, str(GENSAMPLE_SCRIPTS[args.gs_step])]
    cmd += _build_argv(args, GENSAMPLE_SPECS[args.gs_step])
    _run_passthrough(cmd)


# ---------------------------------------------------------------------------
# benchmark
# ---------------------------------------------------------------------------

def _add_benchmark_optimize_parser(b_sub):
    p = b_sub.add_parser(
        "optimize", help="Grid-search redundancy/sampling configs for diversity (01_optimize_benchmark.py).")
    p.add_argument("--full", action="store_true", default=None,
                  help="Consider the full 27-config redundancy grid for Step 1 winner selection.")
    p.add_argument("--full-diagnostic", action="store_true", default=None,
                  help="Grid-search N x method x dims for ALL 27 configs, not just the winner. Implies --full.")
    p.add_argument("--plot-only", action="store_true", default=None,
                  help="Skip grid search; regenerate plot from existing results.tsv.")
    p.add_argument("--summary-tsv", metavar="PATH", default=None,
                  help="Input summary TSV, resolved against scrape/outputs/ (script default: summary_merged.tsv).")
    p.add_argument("--n-fractions", metavar="F1,F2,...", default=None,
                  help="Comma-separated sample sizes N as a fraction of each config's own pool size.")
    p.add_argument("--methods", metavar="M1,M2,...", default=None,
                  help="Comma-separated sampling methods (subset of stratified,kennard-stone,kde,cluster).")
    p.add_argument("--dims-combo", metavar="DIM1,DIM2,...", action="append", default=None,
                  help="Feature-dimension combo, e.g. 'dg,pep_len'. Repeat for multiple combos.")
    p.add_argument("--tm-thresholds", metavar="T1,T2,...", default=None,
                  help="Comma-separated receptor TM-score cutoffs.")
    p.add_argument("--rec-presets", metavar="P1,P2,...", default=None,
                  help="Comma-separated receptor sequence-identity presets.")
    p.add_argument("--pep-presets", metavar="P1,P2,...", default=None,
                  help="Comma-separated peptide sequence-identity presets.")


BENCHMARK_OPTIMIZE_SPEC = [
    ("full", "store_true"), ("full_diagnostic", "store_true"), ("plot_only", "store_true"),
    ("summary_tsv", "value"), ("n_fractions", "value"), ("methods", "value"),
    ("dims_combo", "append"), ("tm_thresholds", "value"), ("rec_presets", "value"), ("pep_presets", "value"),
]


def _load_benchmark_profile_names():
    """Read PROFILES.keys() from 02_make_benchmark.py without hardcoding them
    here, so --profile's choices can never drift from the real script. Returns
    None (free-form flag, no choices validation) if the script is unreadable."""
    make_py = BENCHMARK_SCRIPTS["make"]
    if not make_py.exists():
        return None
    try:
        import importlib.util
        spec = importlib.util.spec_from_file_location("_bench_profiles_mod", make_py)
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)
        return sorted(mod.PROFILES)
    except Exception:
        return None


def _add_benchmark_make_parser(b_sub):
    p = b_sub.add_parser(
        "make", help="Generate a versioned pepxpro_benchmark release (02_make_benchmark.py).")
    group = p.add_mutually_exclusive_group(required=True)
    group.add_argument("--profile", metavar="NAME", choices=_load_benchmark_profile_names(), default=None,
                       help="Frozen, citation-safe benchmark profile (e.g. GC2 -> v1).")
    group.add_argument("--from-config", metavar="auto|GCn", default=None,
                       help="Derive parameters live from optimization/results.tsv.")
    p.add_argument("--version", metavar="VX", default=None,
                  help="Version label, e.g. v3 (required for --from-config).")
    p.add_argument("--seed", type=int, metavar="N", default=None,
                  help="Override the sampling seed.")
    p.add_argument("--dry-run", action="store_true", default=None,
                  help="Show filtered pool size without generating files.")


BENCHMARK_MAKE_SPEC = [
    ("profile", "value"), ("from_config", "value"), ("version", "value"),
    ("seed", "value"), ("dry_run", "store_true"),
]

BENCHMARK_SPECS = {
    "optimize": BENCHMARK_OPTIMIZE_SPEC,
    "make": BENCHMARK_MAKE_SPEC,
}


def run_benchmark(args):
    cmd = [sys.executable, str(BENCHMARK_SCRIPTS[args.bench_step])]
    cmd += _build_argv(args, BENCHMARK_SPECS[args.bench_step])
    _run_passthrough(cmd)


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

    p_gensample = sub.add_parser("genSample", help="Filter, sample, and export a benchmark set.")
    gs_sub = p_gensample.add_subparsers(dest="gs_step", required=True)
    _add_gensample_gen_parser(gs_sub)
    _add_gensample_methods_parser(gs_sub)
    _add_gensample_plot_parser(gs_sub)
    _add_gensample_view_parser(gs_sub)

    p_benchmark = sub.add_parser("benchmark", help="Grid-search and build a versioned benchmark release.")
    b_sub = p_benchmark.add_subparsers(dest="bench_step", required=True)
    _add_benchmark_optimize_parser(b_sub)
    _add_benchmark_make_parser(b_sub)

    p_config = sub.add_parser("config", help="Show pipeline default parameters.")
    p_config.add_argument("--print", dest="print_params", action="store_true",
                          help="Print all default parameters for every pipeline stage.")

    args = parser.parse_args()
    if args.command == "scrape":
        run_scrape(args)
    elif args.command == "genSample":
        run_gensample(args)
    elif args.command == "benchmark":
        run_benchmark(args)
    elif args.command == "config":
        if args.print_params:
            run_config(args)
        else:
            p_config.print_help()


if __name__ == "__main__":
    main()
