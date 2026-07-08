import json
import shutil
from pathlib import Path

from .utils import save_tsv

SUMMARY_COLS = ["pdb_id", "receptor", "peptide", "dG_kcal_mol", "binding_type",
                "peptide_ss_class", "receptor_func_class"]


def next_run_dir(outputs_dir):
    """Return the next available outputs/run_XX/ path."""
    outputs_dir = Path(outputs_dir)
    outputs_dir.mkdir(parents=True, exist_ok=True)
    existing = sorted(outputs_dir.glob("run_*"))
    if not existing:
        return outputs_dir / "run_01"
    last = existing[-1].name  # e.g. "run_03"
    try:
        num = int(last.split("_")[1]) + 1
    except (IndexError, ValueError):
        num = len(existing) + 1
    return outputs_dir / f"run_{num:02d}"


def write_sample(rows, run_dir, columns=None):
    run_dir = Path(run_dir)
    run_dir.mkdir(parents=True, exist_ok=True)
    out = run_dir / "sample.tsv"
    save_tsv(rows, out, columns or SUMMARY_COLS)
    print(f"  Written {out}")
    return out


def save_run_config(config, run_dir):
    run_dir = Path(run_dir)
    run_dir.mkdir(parents=True, exist_ok=True)
    out = run_dir / "run_config.json"
    with open(out, "w", encoding="ascii") as fh:
        json.dump(config, fh, indent=2)
    print(f"  Written {out}")
    return out


def _write_structure(topology, positions, dst, fmt):
    """Write topology+positions to dst in 'cif' or 'pdb' format."""
    if fmt == "pdb":
        from openmm.app import PDBFile
        with open(dst, "w", encoding="ascii") as fh:
            PDBFile.writeFile(topology, positions, fh)
    else:
        from openmm.app import PDBxFile
        with open(dst, "w", encoding="ascii") as fh:
            PDBxFile.writeFile(topology, positions, fh)


def link_cifs(rows, run_dir, fixed_cif_dir, copy=False, fmt="cif"):
    """Symlink (or copy) structure files for selected rows into run_dir/cif/ or run_dir/pdb/."""
    run_dir = Path(run_dir)
    fixed_cif_dir = Path(fixed_cif_dir)
    out_dir = run_dir / fmt
    out_dir.mkdir(parents=True, exist_ok=True)

    missing = []
    for r in rows:
        src = fixed_cif_dir / f"{r['pdb_id']}.cif"
        dst = out_dir / f"{r['pdb_id']}.{fmt}"
        if not src.exists():
            missing.append(r["pdb_id"])
            continue
        if dst.exists() or dst.is_symlink():
            dst.unlink()
        if fmt == "pdb":
            # CIF → PDB conversion via OpenMM
            from pdbfixer import PDBFixer
            try:
                fixer = PDBFixer(filename=str(src))
                _write_structure(fixer.topology, fixer.positions, dst, "pdb")
            except Exception as exc:
                print(f"  Warning: CIF→PDB conversion failed for {r['pdb_id']}: {exc}")
                missing.append(r["pdb_id"])
                continue
        elif copy:
            shutil.copy2(src, dst)
        else:
            dst.symlink_to(src.resolve())

    action = "Converted to PDB" if fmt == "pdb" else ("Copied" if copy else "Linked")
    n_ok = len(rows) - len(missing)
    print(f"  {action} {n_ok} files → {out_dir}")
    if missing:
        print(f"  Warning: {len(missing)} file(s) not found or failed: {', '.join(missing)}")


def write_fasta(rows, run_dir):
    """Write peptide_sequences.fasta and receptor_sequences.fasta for selected rows."""
    run_dir = Path(run_dir)
    run_dir.mkdir(parents=True, exist_ok=True)

    for seq_col, fname in [("peptide_seq", "peptide_sequences.fasta"),
                            ("receptor_seq", "receptor_sequences.fasta")]:
        out = run_dir / fname
        written = 0
        missing = 0
        with open(out, "w", encoding="ascii") as fh:
            for r in rows:
                seq = r.get(seq_col, "").strip()
                if not seq or seq == "NA":
                    missing += 1
                    continue
                fh.write(f">{r['pdb_id']}\n{seq}\n")
                written += 1
        print(f"  Written {out}  ({written} sequences)", end="")
        if missing:
            print(f"  [{missing} missing seq skipped]", end="")
        print()


def protonate_cifs(rows, run_dir, fixed_cif_dir, ph=7.0, fmt="cif"):
    """Copy structure files into run_dir/cif/ (or pdb/) and add missing hydrogens at the given pH."""
    from pdbfixer import PDBFixer

    run_dir = Path(run_dir)
    fixed_cif_dir = Path(fixed_cif_dir)
    out_dir = run_dir / fmt
    out_dir.mkdir(parents=True, exist_ok=True)

    missing = []
    done = 0
    for r in rows:
        src = fixed_cif_dir / f"{r['pdb_id']}.cif"
        dst = out_dir / f"{r['pdb_id']}.{fmt}"
        if not src.exists():
            missing.append(r["pdb_id"])
            continue
        try:
            fixer = PDBFixer(filename=str(src))
            fixer.addMissingHydrogens(ph)
            _write_structure(fixer.topology, fixer.positions, dst, fmt)
            done += 1
        except Exception as exc:
            print(f"  Warning: protonation failed for {r['pdb_id']}: {exc}")
            shutil.copy2(src, out_dir / f"{r['pdb_id']}.cif")

    label = fmt.upper()
    print(f"  Protonated {done} {label} files at pH {ph} → {out_dir}")
    if missing:
        print(f"  Warning: {len(missing)} file(s) not found: {', '.join(missing)}")
