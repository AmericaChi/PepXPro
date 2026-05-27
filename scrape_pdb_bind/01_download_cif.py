"""
Step 1: Parse INDEX_general_PL.2020, identify N-mer peptide binders,
compute DeltaG from binding data, and download biological assembly CIFs
(tries assemblies 1, 2, 3, ... stopping at the first 404) to pdb_cache/.
Pre-rejects covalent complexes and entries without Kd or Ki data.

Output: peptide_entries.tsv, pdb_cache/<pdb_id>-<N>.cif
Run:    mamba activate mdtools && python 01_download_cif.py 2>&1 | tee step1_run.log
"""

import re
import math
import gzip
import shutil
import argparse
import urllib.request
import urllib.error
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

INDEX_FILE = Path("INDEX_general_PL.2020")
PDB_CACHE  = Path("pdb_cache")
OUT_TSV    = Path("peptide_entries.tsv")
MAX_WORKERS = 8

R_KCAL = 0.001987   # kcal / (mol * K)
T      = 300.0      # K

UNIT_FACTORS = {"pM": 1e-12, "nM": 1e-9, "uM": 1e-6, "mM": 1e-3, "M": 1.0}

# Entries with these flags in the comment cannot be simulated with standard MD.
REJECT_FLAGS = ["covalent complex"]

_MER_PAT     = re.compile(r'\((\d+)-mer\)')
_BINDING_PAT = re.compile(r'(Kd|Ki|IC50)([<>=~]*)([\d.]+)(pM|nM|uM|mM|M)')


def parse_binding(token):
    """Parse e.g. 'Ki=0.43uM'.
    Returns (binding_type, conc_in_M) or None.
    Only exact measurements (operator '=' or absent) are accepted;
    inequalities (<, <=, >, >=) and approximations (~) return None.
    """
    m = _BINDING_PAT.match(token)
    if not m:
        return None
    btype, op, val, unit = m.group(1), m.group(2), float(m.group(3)), m.group(4)
    if op not in ("", "="):
        return None
    return btype, val * UNIT_FACTORS[unit]


def calc_dG(conc_M):
    """DeltaG = RT * ln(Kd)  [kcal/mol, negative for tight binders]."""
    return R_KCAL * T * math.log(conc_M)


def parse_index(path):
    """Return list of dicts for every N-mer peptide entry in the INDEX."""
    entries = []
    with open(path, encoding="ascii") as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue

            # Must contain '(N-mer)'
            mer_m = _MER_PAT.search(line)
            if not mer_m:
                continue

            n_mer = int(mer_m.group(1))
            parts = line.split()
            pdb_id = parts[0].lower()
            binding_str = parts[3]  # e.g. 'Ki=0.43uM'

            # Reject entries that cannot be simulated with standard MD
            comment = line[line.find("//"):].lower() if "//" in line else ""
            reject_reason = next(
                (f for f in REJECT_FLAGS if f in comment), None
            )

            parsed = parse_binding(binding_str)
            if parsed is not None:
                btype, conc_M = parsed
                dG = calc_dG(conc_M)
                # IC50 is not a true equilibrium constant; flag as approximate
                dG_approx = btype == "IC50"
            else:
                btype, dG, dG_approx = None, None, False

            entries.append(
                dict(
                    pdb_id=pdb_id,
                    n_mer=n_mer,
                    binding_type=btype or "NA",
                    dG_kcal=f"{dG:.3f}" if dG is not None else "NA",
                    dG_approx="Y" if dG_approx else "N",
                    step1_reject=reject_reason or "",
                )
            )
    return entries


def download_all(pdb_id, cache_dir):
    """Download all available biological assembly CIFs for pdb_id.
    Tries assemblies 1, 2, 3, ... stopping on the first 404.
    Files are saved as <pdb_id>-<n>.cif.
    Returns (pdb_id, {asm_num: status_str}).
    """
    results = {}
    # Collect already-cached assemblies first.
    for cif_path in cache_dir.glob(f"{pdb_id}-*.cif"):
        try:
            results[int(cif_path.stem.rsplit("-", 1)[1])] = "cached"
        except (ValueError, IndexError):
            pass
    # Only attempt to download from the next assembly after the highest cached.
    # This avoids a network round-trip on every rerun for already-complete entries.
    start_n = max(results.keys()) + 1 if results else 1
    for n in range(start_n, 10):
        out = cache_dir / f"{pdb_id}-{n}.cif"
        url = f"https://files.rcsb.org/download/{pdb_id.upper()}-assembly{n}.cif.gz"
        gz_path = cache_dir / f"{pdb_id}-{n}.cif.gz"
        try:
            urllib.request.urlretrieve(url, gz_path)
            with gzip.open(gz_path, "rb") as gz, open(out, "wb") as f:
                shutil.copyfileobj(gz, f)
            gz_path.unlink()
            results[n] = "ok"
        except urllib.error.HTTPError as exc:
            if gz_path.exists():
                gz_path.unlink()
            if exc.code == 404:
                break  # no more assemblies for this entry
            results[n] = f"error: {exc}"
        except Exception as exc:
            if gz_path.exists():
                gz_path.unlink()
            results[n] = f"error: {exc}"
    return pdb_id, results


def main():
    parser = argparse.ArgumentParser(
        description="Step 1: Parse PDBbind INDEX for N-mer peptide entries, compute ΔG, "
                    "and download all available biological assembly CIFs to pdb_cache/ "
                    "(saved as <pdb_id>-<N>.cif).")
    parser.add_argument("--clean", action="store_true",
                        help="Remove files generated by this step (peptide_entries.tsv).")
    args = parser.parse_args()

    if args.clean:
        if OUT_TSV.exists():
            OUT_TSV.unlink()
            print(f"Removed {OUT_TSV}")
        return

    PDB_CACHE.mkdir(exist_ok=True)

    print(f"Parsing {INDEX_FILE} ...")
    entries = parse_index(INDEX_FILE)
    print(f"  {len(entries)} peptide-mer entries found")

    # Keep only entries with a real equilibrium binding constant (Kd or Ki).
    # IC50 cannot be reliably converted to ΔG without assay-specific parameters,
    # and NA entries have no usable binding data.
    before = len(entries)
    entries = [e for e in entries if e["binding_type"] in ("Kd", "Ki")]
    print(f"  {before - len(entries)} dropped (IC50 or no binding data)")
    print(f"  {len(entries)} entries remaining (Kd/Ki only)")

    # Count pre-rejected (covalent etc.)
    rejected = [e for e in entries if e["step1_reject"]]
    print(f"  {len(rejected)} rejected at step 1 (e.g. covalent complex)")
    downloadable = [e for e in entries if not e["step1_reject"]]
    print(f"  {len(downloadable)} to download")

    # Write entries TSV (all entries, including step1 rejects)
    cols = ["pdb_id", "n_mer", "binding_type", "dG_kcal", "dG_approx", "step1_reject"]
    with open(OUT_TSV, "w", encoding="ascii") as fh:
        fh.write("\t".join(cols) + "\n")
        for e in entries:
            fh.write("\t".join(str(e[c]) for c in cols) + "\n")
    print(f"  Written {OUT_TSV}")

    # Download all available assemblies concurrently.
    pdb_ids = [e["pdb_id"] for e in downloadable]
    counts = {"ok": 0, "cached": 0, "error": 0}
    errors = []

    with ThreadPoolExecutor(max_workers=MAX_WORKERS) as pool:
        futures = {pool.submit(download_all, pid, PDB_CACHE): pid for pid in pdb_ids}
        for i, fut in enumerate(as_completed(futures), 1):
            pid, asm_results = fut.result()
            statuses = list(asm_results.values())
            if not statuses or all(s.startswith("error") for s in statuses):
                counts["error"] += 1
                msg = statuses[0] if statuses else "no assemblies found"
                errors.append((pid, msg))
                status_str = msg
            elif any(s == "ok" for s in statuses):
                counts["ok"] += 1
                n_asm = len(asm_results)
                status_str = f"ok ({n_asm} {'assembly' if n_asm == 1 else 'assemblies'})"
            else:
                counts["cached"] += 1
                n_asm = len(asm_results)
                status_str = f"cached ({n_asm} {'assembly' if n_asm == 1 else 'assemblies'})"
            if i % 100 == 0 or status_str.startswith("error"):
                print(
                    f"  [{i}/{len(pdb_ids)}] "
                    f"ok={counts['ok']} cached={counts['cached']} "
                    f"errors={counts['error']}  ({pid}: {status_str})"
                )

    print(
        f"\nDone. ok={counts['ok']} cached={counts['cached']} "
        f"errors={counts['error']}"
    )
    if errors:
        print("Download errors:")
        for pid, msg in errors:
            print(f"  {pid}: {msg}")


if __name__ == "__main__":
    main()
