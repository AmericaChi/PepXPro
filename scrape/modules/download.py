import re
import math
import gzip
import shutil
import urllib.request
import urllib.error
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

R_KCAL = 0.001987   # kcal / (mol * K)

UNIT_FACTORS = {"pM": 1e-12, "nM": 1e-9, "uM": 1e-6, "mM": 1e-3, "M": 1.0}

REJECT_FLAGS = ["covalent complex"]

_MER_PAT     = re.compile(r'\((\d+)-mer\)')
_BINDING_PAT = re.compile(r'(Kd|Ki|IC50)([<>=~]*)([\d.]+)(pM|nM|uM|mM|M)')
_TEMP_PAT    = re.compile(r'temp(?:erature)?\s*(?:is\s*)?[=:]?\s*([\d.]+)\s*(?:degree|deg|°|c\b|celsius)?', re.IGNORECASE)

MAX_WORKERS = 8


def parse_binding(token):
    """Parse e.g. 'Ki=0.43uM'. Returns (binding_type, conc_in_M) or None.
    Only exact measurements ('=' or absent) are accepted."""
    m = _BINDING_PAT.match(token)
    if not m:
        return None
    btype, op, val, unit = m.group(1), m.group(2), float(m.group(3)), m.group(4)
    if op not in ("", "="):
        return None
    return btype, val * UNIT_FACTORS[unit]


def parse_temperature(comment):
    """Extract temperature in K from a PDBbind index comment field.

    Handles values in Celsius (< 200) or Kelvin (200-400).
    Returns temperature in K, or None if not found.
    """
    m = _TEMP_PAT.search(comment)
    if not m:
        return None
    try:
        val = float(m.group(1))
    except ValueError:
        return None
    if val < 200:        # assume Celsius
        val += 273.15
    if 200.0 < val < 400.0:
        return val
    return None


def calc_dG(conc_M, t):
    return R_KCAL * t * math.log(conc_M)


def parse_index(path, fallback_temp=300.0):
    """Return list of dicts for every N-mer peptide entry in the INDEX.
    Extracts: pdb_id, n_mer, year, binding_type, dG_kcal, dG_approx,
              step1_reject, temperature.

    Temperature is parsed from the index comment when available;
    falls back to fallback_temp (default 300 K).
    """
    entries = []
    with open(path, encoding="ascii") as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            mer_m = _MER_PAT.search(line)
            if not mer_m:
                continue
            n_mer = int(mer_m.group(1))
            parts = line.split()
            pdb_id = parts[0].lower()
            resolution  = parts[1] if len(parts) > 1 else "NA"
            year        = parts[2] if len(parts) > 2 else "NA"
            binding_str = parts[3] if len(parts) > 3 else ""

            comment = line[line.find("//"):] if "//" in line else ""
            reject_reason = next((f for f in REJECT_FLAGS if f in comment.lower()), None)

            t = parse_temperature(comment) or fallback_temp

            parsed = parse_binding(binding_str)
            if parsed is not None:
                btype, conc_M = parsed
                dG = calc_dG(conc_M, t)
                dG_approx = btype == "IC50"
            else:
                btype, dG, dG_approx = None, None, False

            entries.append(dict(
                pdb_id=pdb_id,
                n_mer=n_mer,
                resolution=resolution,
                year=year,
                binding_type=btype or "NA",
                dG_kcal=f"{dG:.3f}" if dG is not None else "NA",
                dG_approx="Y" if dG_approx else "N",
                step1_reject=reject_reason or "",
                temperature=f"{t:.2f}",
            ))
    return entries


def download_all(pdb_id, cache_dir):
    """Download all available biological assembly CIFs for pdb_id.
    Tries assemblies 1, 2, 3, ... stopping on the first 404.
    Returns (pdb_id, {asm_num: status_str}).
    """
    cache_dir = Path(cache_dir)
    results = {}
    for cif_path in cache_dir.glob(f"{pdb_id}-*.cif"):
        try:
            results[int(cif_path.stem.rsplit("-", 1)[1])] = "cached"
        except (ValueError, IndexError):
            pass
    start_n = max(results.keys()) + 1 if results else 1
    for n in range(start_n, 10):
        out    = cache_dir / f"{pdb_id}-{n}.cif"
        gz_path = cache_dir / f"{pdb_id}-{n}.cif.gz"
        url    = f"https://files.rcsb.org/download/{pdb_id.upper()}-assembly{n}.cif.gz"
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
                break
            results[n] = f"error: {exc}"
        except Exception as exc:
            if gz_path.exists():
                gz_path.unlink()
            results[n] = f"error: {exc}"
    return pdb_id, results


def download_entries(pdb_ids, cache_dir, log=print):
    """Download assemblies for all pdb_ids concurrently. Returns count dict."""
    cache_dir = Path(cache_dir)
    cache_dir.mkdir(exist_ok=True)
    counts = {"ok": 0, "cached": 0, "error": 0}
    errors = []
    with ThreadPoolExecutor(max_workers=MAX_WORKERS) as pool:
        futures = {pool.submit(download_all, pid, cache_dir): pid for pid in pdb_ids}
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
                status_str = f"ok ({len(asm_results)} assemblies)"
            else:
                counts["cached"] += 1
                status_str = f"cached ({len(asm_results)} assemblies)"
            if i % 100 == 0 or status_str.startswith("error"):
                log(f"  [{i}/{len(pdb_ids)}] ok={counts['ok']} cached={counts['cached']} "
                    f"errors={counts['error']}  ({pid}: {status_str})")
    return counts, errors
