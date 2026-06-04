# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Environment

Activate the conda environment before running any scripts:
```bash
conda  activate peppro
```

To create the environment from scratch:
```bash
conda create -n peppro -c conda-forge python=3.12 openmm matplotlib numpy pandas biopython pdbfixer scipy scikit-learn pymol-open-source
```

Key packages and their roles:
- `biopython` — `MMCIF2Dict` for parsing raw CIF files in step 2
- `openmm` — `PDBxFile`, `Modeller` for writing fixed CIFs
- `pdbfixer` — fills missing residues and atoms
- `scipy` — `cKDTree` for spatial distance queries
- `numpy` — array operations
- `pymol-open-source` / `pymol2` — headless RMSD alignment (step 4 `--sample`) and interactive viewer

## Running the Pipeline

Run steps sequentially, logging output in real time (`-u` disables Python's output buffering):
```bash
python -u 01_download_cif.py 2>&1 | tee step01_run.log
python -u 02_filter_peptide.py 2>&1 | tee step02_run.log
python -u 03_summary.py 2>&1 | tee step03_run.log
python -u 04_sample_and_view.py --sample 20 2>&1 | tee step04_run.log
python -u 04_sample_and_view.py --sample 20 --single-chain 2>&1 | tee step04_run.log
python -u 04_sample_and_view.py --expand 10 2>&1 | tee step04_run.log
python 04_sample_and_view.py --view-sample                 # browse only the sampled structures
python 04_sample_and_view.py                               # browse all of fixed_cif/
```

Each script accepts `--clean` to remove its output files (never touches `pdb_cache/`).

## Architecture

Four-step pipeline that extracts peptide-protein binding data from PDBbind:

```
INDEX_general_PL.2020 (19,449 entries)
    → 01_download_cif.py → peptide_entries.tsv + pdb_cache/<pdb_id>-<N>.cif
    → 02_filter_peptide.py → filtered_systems.tsv + fixed_cif/*.cif
    → 03_summary.py → summary.tsv
    → [04_sample_and_view.py] → diversity sampler + PyMOL viewer (optional)
```

**Step 1 (`01_download_cif.py`):** Parses INDEX for `(N-mer)` entries, computes ΔG = RT·ln(Kd) at 300 K, downloads all available biological assembly CIFs from RCSB (tries assemblies 1, 2, 3, ... stopping at the first 404; 8 concurrent workers), pre-rejects covalent complexes and entries without Kd or Ki data (IC50 and NA are dropped). Files saved as `pdb_cache/<pdb_id>-<N>.cif`.

**Step 2 (`02_filter_peptide.py`):** The main filtering and structure preparation step. For each entry:
1. Scans all available assembly CIFs via MMCIF2Dict (no PDBFixer at this stage), filters to those containing at least one protein chain with full entity sequence length ≤ `max(50, n_mer*2)` (i.e., a plausible peptide), then selects the qualifying assembly with the fewest chains (ties broken by assembly number); the chosen assembly number is recorded in `filtered_systems.tsv`
2. Identifies the peptide chain: filters protein chains to those with full entity sequence length ≤ `max(50, n_mer*2)` (structural prior: peptide is always the shortest chain), then picks the one closest to n_mer. No hard tolerance — the upper bound is the meaningful filter. Full entity sequence length (from `_entity_poly.pdbx_seq_one_letter_code` via `_atom_site.label_entity_id`) is used so chains with missing residues are not penalised; both `auth_asym_id`- and `label_asym_id`-keyed dicts are built and the minimum entity length ≥ observed count is used to handle PDBFixer's inconsistent choice of chain ID scheme across files
3. Checks that PDBFixer loaded the same residue count as the raw CIF for the peptide chain; if fewer, non-standard residues were silently dropped → rejected
4. Checks peptide residues against the AMBER14 whitelist (20 standard AAs + HID/HIE/HIP)
5. Excludes protein chains with ≥70% sequence similarity to the peptide (prevents symmetry-related peptide copies from being classified as receptor)
6. Finds receptor chains with ≥1 atom within 4 Å of the peptide; rejects if any receptor chain has an observed residue numbering gap > 10
6a. Rejects if the peptide C-terminus (backbone C atom) is within 1.8 Å of any receptor N-terminus, or peptide N-terminus within 1.8 Å of any receptor C-terminus — indicates the two auth chains are actually one continuous polypeptide
7. Checks receptor residues against the AMBER14 whitelist
8. Finds AMBER14/TIP3P-FB-compatible ion chains (NA, K, CL, CA, MG, ZN, FE/FE2, CU, MN) within 4 Å of the system; these are kept in the fixed CIF but not listed in the receptor column of the TSV
9. Writes `fixed_cif/<pdb_id>.cif`: peptide + receptor chains with all waters removed; missing residues filled (all for peptide, non-terminal only for receptor); AMBER14 ions retained. `PDBxFile` assigns sequential chain IDs (A, B, C, ...) regardless of assembly chain IDs — the screen output and TSV reflect these final IDs.
10. Rejects if any backbone CA-CA distance > 5 Å after fixing (indicates an unfilled gap PDBFixer could not repair)

Output lines: `PASS` includes `asm=N  pep=chainID:obs+added  rec=chainID:obs+added  ions=[CA,ZN,...]`; `REJECT` includes the rejection reason. Both include running `(accepted=N rejected=N)` counters. A `WARNING` line is printed below any PASS where the post-fix peptide length differs from the INDEX n_mer.

`filtered_systems.tsv` columns (accepted rows): `pdb_id`, `n_mer`, `accepted`, `assembly_num`, `peptide_chain`, `peptide_length`, `peptide_added`, `receptor_chains`, `receptor_length`, `receptor_added`, `reject_reason`. `peptide_added`/`receptor_added` are the residue counts PDBFixer added; rejected rows have empty strings for all chain/length/added fields.

**Step 3 (`03_summary.py`):** Joins `filtered_systems.tsv` with `peptide_entries.tsv`, formats receptor as `"AB:300"` (chain IDs concatenated + `:length`) and peptide as `"C:12"` (actual post-fix residue count). Rows are sorted strong-to-weak by ΔG. Columns: `pdb_id`, `receptor`, `peptide`, `dG_kcal_mol`, `binding_type`. n_mer is omitted — the peptide column already reflects the true built length.

**Step 4 (`04_sample_and_view.py`):** Four modes:
- No flags: PyMOL viewer over all of `fixed_cif/`; directory re-globbed on each keypress so new files from a running step 2 appear automatically.
- `--sample N`: selects N structures from `summary.tsv` using diversity sampling across ΔG and peptide length. Both axes are divided into 6 quantile-based bins (6×6 = 36 cells); one random entry is drawn from each occupied cell in random order, with remaining slots filled randomly. Then iteratively replaces the most receptor-similar pair (CA RMSD via headless `pymol2` `cmd.align`) until all comparable pairs exceed `--receptor-rmsd` (default 0, disabled). Writes `sample.tsv`; each run produces a different selection. Does not open PyMOL. Add `--single-chain` to restrict the pool to entries with a single receptor chain before sampling. Use `--max-receptor N`, `--max-peptide N`, or `--max-added N` to cap pool size before sampling (`--max-added` reads `peptide_added` from `filtered_systems.tsv` directly). Aborts if `sample.tsv` already exists unless `--overwrite` is also passed. After sampling, always computes the full pairwise receptor RMSD matrix and prints a cluster report (see below).
- `--expand N`: adds N new structures to an existing `sample.tsv`. Candidates are drawn from `summary.tsv` (excluding already-sampled entries) using the same diversity sampling, then the RMSD diversity filter is applied with the existing sample entries as fixed anchors — new candidates are replaced if too similar to each other or to any existing entry, but existing entries are never evicted. Overwrites `sample.tsv` with the combined list sorted by ΔG. Requires `--overwrite` to proceed (since it writes to `sample.tsv`). `--max-receptor`, `--max-peptide`, and `--max-added` apply to the candidate pool only. `--max-added` reads `peptide_added` from `filtered_systems.tsv` directly. Cluster report is printed for the new entries only.
- `--view-sample`: opens the PyMOL viewer restricted to the systems in `sample.tsv`, ordered by ΔG from tightest to weakest binder.

All viewer modes colour chains with `util.cbc` and use left/right arrow keys to navigate.

**Cluster report** (printed after `--sample` and `--expand`): computes all N*(N-1)/2 pairwise receptor CA-RMSDs via headless PyMOL `cmd.align(cycles=0)` with progress output. Prints overall min/mean/max, then a table at 0.5/1/2/3/5 Å thresholds showing number of clusters, largest cluster size, and singleton count. Clustering uses connected components (single-linkage): structures linked by any path of pairs within the cutoff are one cluster. `stdout` is duplicated before PyMOL's fd-level suppressor so progress lines reach the terminal.

## Key Conventions

- All file I/O uses `encoding='ascii'`
- ΔG calculated with `R = 0.001987 kcal/(mol·K)`, `T = 300.0 K`
- Spatial distances are in Ångströms throughout
- The related expanded pipeline lives at `../pdbbind_peptide/` and references `../scrape_pdb_bind/INDEX_general_PL.2020`
