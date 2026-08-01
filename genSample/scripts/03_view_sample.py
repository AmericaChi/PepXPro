"""
Step 5: Browse the fixed CIFs of a genSample run in PyMOL, with named
peptide/receptor selections.

Each structure gets two selections built from its sample.tsv row — "peptide"
and "receptor" — so that a single chain-agnostic preset works on every entry
regardless of which literal chain letter happens to be the peptide for that
PDB ID. The PepXPro preset (white background, cartoon only, gray85 receptor /
orange peptide) is baked into the loader itself and re-applied on every
navigation step, since each step deletes and reloads the whole scene.
Navigate with left/right arrow keys. Structures are shown in ΔG order
(tightest binder first).

Usage:
  python 03_view_sample.py                          # latest run
  python 03_view_sample.py --run 3                   # run_03/sample.tsv
  python 03_view_sample.py --run 3 --cif-dir /path/to/cif
  python 03_view_sample.py --run-dir ../../pepxpro_benchmark/v1   # any sample.tsv dir

  CIF lookup order: <run dir>/cif/  →  <run dir>/pdb/  →  scrape/fixed_cif/
  Override all with --cif-dir.
"""

import argparse
import json
import os
import sys
import tempfile
from pathlib import Path

_here   = Path(__file__).resolve().parent
_root   = _here.parent                          # genSample/
_scrape = _root.parent / "scrape"
sys.path.insert(0, str(_root))

from modules.utils import load_tsv, parse_receptor, parse_peptide

OUTPUTS_DIR  = _root / "outputs"
FALLBACK_CIF = _scrape / "fixed_cif"


# ── helpers ───────────────────────────────────────────────────────────────────

def latest_run_dir():
    runs = sorted(OUTPUTS_DIR.glob("run_*"))
    return runs[-1] if runs else None


def resolve_run_dir(run_n):
    if run_n is None:
        d = latest_run_dir()
        if d is None:
            print("No run directories found in genSample/outputs/.")
            sys.exit(1)
        return d
    return OUTPUTS_DIR / f"run_{run_n:02d}"


def resolve_cif(pdb_id, run_dir, fallback):
    for sub in ("cif", "pdb"):
        p = run_dir / sub / f"{pdb_id}.cif"
        if p.exists():
            return p
    p = fallback / f"{pdb_id}.cif"
    return p if p.exists() else None


def find_pymol():
    """Return path to the pymol executable in the active conda env, or just 'pymol'."""
    conda_prefix = os.environ.get("CONDA_PREFIX", "")
    if conda_prefix:
        for candidate in (
            Path(conda_prefix) / "bin" / "pymol",
            Path(conda_prefix) / "Scripts" / "pymol.exe",  # Windows
        ):
            if candidate.exists():
                return str(candidate)
    return "pymol"


def chain_selection(chain_str):
    """'AB' -> 'A+B' (PyMOL multi-chain selection syntax). Empty on bad input."""
    if not chain_str:
        return ""
    return "+".join(list(chain_str))


# ── PyMOL script generator ────────────────────────────────────────────────────

_VIEWER_TEMPLATE = """\
import json
from pymol import cmd

_data     = __DATA__
_cifs     = _data["cifs"]
_labels   = _data["labels"]
_pep_sels = _data["pep_sels"]
_rec_sels = _data["rec_sels"]
_idx      = [0]

def _load():
    path   = _cifs[_idx[0]]
    pdb_id = path.rsplit("/", 1)[-1].rsplit(".", 1)[0]
    cmd.delete("all")
    cmd.load(path, pdb_id)

    # Named, chain-agnostic selections: same preset works for every entry
    # regardless of which literal chain letter is the peptide this time.
    cmd.delete("peptide")
    cmd.delete("receptor")
    pep_sel = _pep_sels[_idx[0]]
    rec_sel = _rec_sels[_idx[0]]
    if pep_sel:
        cmd.select("peptide", "chain %s" % pep_sel)
    if rec_sel:
        cmd.select("receptor", "chain %s" % rec_sel)
    cmd.deselect()

    # PepXPro preset — re-applied on every navigation step, since _load()
    # deletes and reloads the whole scene each time.
    cmd.bg_color("white")
    cmd.hide("everything", "all")

    cmd.show("cartoon", "receptor")
    cmd.show("cartoon", "peptide")

    cmd.color("gray85", "receptor")
    cmd.color("0x44bf70", "peptide")  # green, matching the peptide-length-distribution color (viridis t=0.7)

    cmd.set("cartoon_fancy_helices", 1)
    # smoothing off: with short peptides (few CA control points) the loop
    # spline over-smooths coil into a near-straight line — keep the raw trace.
    cmd.set("cartoon_smooth_loops", 0)
    cmd.set("orthoscopic", 1)
    cmd.set("depth_cue", 0)
    cmd.set("ray_shadows", 0)
    cmd.set("ray_opaque_background", 0)

    cmd.set("ambient", 0.3)
    cmd.set("direct", 0.8)
    cmd.set("specular", 0.2)
    cmd.set("shininess", 10)

    cmd.set("surface_quality", 1)
    cmd.set("surface_smooth_edges", 1)
    cmd.set("solvent_radius", 1.4)

    cmd.set("ray_trace_mode", 1)
    cmd.set("ray_trace_gain", 0.1)

    cmd.orient("all")
    cmd.zoom("all")

    print("[%d/%d]  %s" % (_idx[0] + 1, len(_cifs), _labels[_idx[0]]))

def _next():
    if _idx[0] < len(_cifs) - 1:
        _idx[0] += 1
    _load()

def _prev():
    if _idx[0] > 0:
        _idx[0] -= 1
    _load()

cmd.set_key("right", _next)
cmd.set_key("left",  _prev)
print("  left/right arrows to navigate  (%d structures)" % len(_cifs))
print("  'peptide' and 'receptor' selections are rebuilt on every navigation step")
_load()
"""


def write_viewer_script(cif_files, labels, pep_sels, rec_sels):
    data = {
        "cifs":     [str(p) for p in cif_files],
        "labels":   labels,
        "pep_sels": pep_sels,
        "rec_sels": rec_sels,
    }
    content = _VIEWER_TEMPLATE.replace("__DATA__", json.dumps(data))
    tmp = tempfile.NamedTemporaryFile(
        suffix=".py", mode="w", delete=False, encoding="utf-8"
    )
    tmp.write(content)
    tmp.close()
    return tmp.name


# ── main ──────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Browse genSample run structures in PyMOL, with named "
                     "peptide/receptor selections.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument("--run", type=int, metavar="N",
                        help="Run number (default: latest run).")
    parser.add_argument("--run-dir", metavar="PATH",
                        help="Full path to any directory containing sample.tsv "
                             "(e.g. a pepxpro_benchmark version). Overrides --run.")
    parser.add_argument("--cif-dir", metavar="DIR",
                        help="Override CIF directory.")
    args = parser.parse_args()

    run_dir    = Path(args.run_dir) if args.run_dir else resolve_run_dir(args.run)
    sample_tsv = run_dir / "sample.tsv"

    if not run_dir.exists():
        print(f"Run directory not found: {run_dir}")
        sys.exit(1)
    if not sample_tsv.exists():
        print(f"sample.tsv not found in {run_dir}")
        sys.exit(1)

    fallback = Path(args.cif_dir) if args.cif_dir else FALLBACK_CIF
    # If --cif-dir is relative, resolve from the run directory
    if args.cif_dir and not fallback.is_absolute():
        fallback = run_dir / fallback

    rows = load_tsv(sample_tsv)
    for r in rows:
        try:
            r["_dg"] = float(r.get("dG_kcal_mol", "0"))
        except (ValueError, TypeError):
            r["_dg"] = 0.0
    rows.sort(key=lambda r: r["_dg"])

    cif_files, labels, pep_sels, rec_sels, missing = [], [], [], [], []
    for r in rows:
        cif = resolve_cif(r["pdb_id"], run_dir, fallback)
        if cif:
            cif_files.append(cif)
            dg  = r.get("dG_kcal_mol", "?")
            pep = r.get("peptide", "?")
            rec = r.get("receptor", "?")
            labels.append(f"{r['pdb_id']}  dG={dg}  pep={pep}  rec={rec}")

            try:
                pep_chain, _ = parse_peptide(r.get("peptide", ""))
            except (ValueError, IndexError):
                pep_chain = ""
            try:
                rec_chain, _ = parse_receptor(r.get("receptor", ""))
            except (ValueError, IndexError):
                rec_chain = ""
            pep_sels.append(chain_selection(pep_chain))
            rec_sels.append(chain_selection(rec_chain))
        else:
            missing.append(r["pdb_id"])

    if missing:
        print(f"  Warning: CIF not found for {len(missing)} entries: {', '.join(missing)}")

    print(f"  Run:    {run_dir.name}")
    print(f"  Sample: {len(rows)} entries  ({len(cif_files)} CIFs resolved)")

    if not cif_files:
        print("No CIF files to display.")
        sys.exit(0)

    script_path = write_viewer_script(cif_files, labels, pep_sels, rec_sels)
    pymol_bin   = find_pymol()

    import subprocess
    try:
        subprocess.run([pymol_bin, "-r", script_path], check=False)
    finally:
        try:
            os.unlink(script_path)
        except OSError:
            pass


if __name__ == "__main__":
    main()
