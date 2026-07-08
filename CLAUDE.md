# PepXPro

PepXPro is a modular framework for generating peptide–protein benchmark datasets from PDBbind.

## Repository structure

- `scrape/`
  Modular pipeline for downloading, filtering, and preparing peptide–protein systems from PDBbind.
  Modules in `modules/`, scripts in `scripts/`, outputs in `outputs/`, logs in `logs/`.

- `genSample/`
  Modular benchmark/sample generation framework. Reads from `scrape/outputs/`.

- `pepxpro_benchmark/`
  Curated benchmark releases. Scripts in `scripts/`, versioned outputs in `v1/`, `v2/`, etc.
  Uses genSample modules directly — no run_XX created in genSample/outputs/.
  Workflow: tune parameters in genSample first, then run `python scripts/01_make_benchmark.py`.

- `manual/`
  Comprehensive user manual: `pepxpro_manual.txt`.

## Development philosophy

- Reproducible workflows
- Modular analysis scripts
- Diversity-aware sampling
- Sequence and structural redundancy reduction
- Benchmark-oriented dataset generation
- ML-friendly dataset splits

## pepxpro.py — top-level runner

```bash
python pepxpro.py scrape -v                          # full pipeline
python pepxpro.py scrape -v --strip-caps             # strip terminal caps not at interface
python pepxpro.py scrape -v --no-convert-mse         # disable MSE→MET (on by default)
python pepxpro.py config --print                     # show all default parameters
```

## Implemented genSample modules

- Quality filtering (length, dG, resolution, binding type, added residues, exclusion list)
- Sequence redundancy filtering (CD-HIT for receptor at ≥40%, Biopython pairwise for peptide)
- Structural redundancy filtering (TM-score receptor via tmtools, parallelized)
- Stratified, Kennard-Stone, KDE, isolation-forest, cluster sampling
- Structure export (CIF/PDB, symlink/copy, protonation)
- Sampling benchmark plots, per-run comparison plots, FASTA export

## Planned genSample modules

- Functional redundancy filtering
- Peptide feature analysis (charge, hydrophobicity, aromaticity)
- Benchmark export with ML-friendly splits

## ML benchmark — ml_benchmark/notebooks/ml_benchmark.ipynb

**Prerequisites:** the full scrape pipeline must have been run first:
```bash
cd scrape/scripts
python 01_download.py
python 02_filter.py   # adds n_contacts and n_hbonds via compute_interface_contacts()
python 03_summary.py  # produces summary.tsv — the notebook's input
```

Goal: show that PepXPro curation choices affect ΔG prediction quality.
Uses Random Forest with 12 structural features (all from scrape pipeline).

**12 features:** `pep_len`, `rec_len`, `peptide_pct_helix/beta/coil`,
`receptor_pct_helix/beta/coil`, `n_contacts` (4.5 Å), `n_hbonds` (N/O/S 3.5 Å),
`binding_kd`, `resolution`.

**6 strategies per scenario:** Random, Seq-70%+KS-2D/4D/8D/12D, Seq-70%+Rec-40%+KS-12D.

**Scenario A — same-model, varied test curation:** random 50/50 split; all strategies train
on the same half; the other half is curated with each strategy to select N test entries.
**Scenario B — partial-helix held out:** ~96 partial-helix entries as fixed test set;
each strategy curates the remaining entries for training.

**N_TEST (A):** 150, 75. **N_TRAIN (B):** 300, 150, 75.
**Metrics:** RMSE, MAE, Pearson r, Spearman ρ (±std over 25 seeds).

**Interface features in pipeline:** `n_contacts` and `n_hbonds` computed in
`scrape/modules/filter.py::compute_interface_contacts()` from fixed CIFs.
They flow through `filtered_systems.tsv` → `summary.tsv` → genSample `enrich_rows()`.
Available as KS dims: `n_contacts`, `n_hbonds` (see genSample/CLAUDE.md).

## Output organization

scrape/outputs/
├── peptide_entries.tsv
├── filtered_systems.tsv
├── summary.tsv
└── plots/
    ├── summary_distributions.png
    └── summary_statistics.tsv

genSample/outputs/run_XX/
├── sample.tsv
├── run_config.json
├── cif/  or  pdb/
└── plots/
    └── sampling_benchmark.png

Large generated files should not be tracked by git.

## Environment

Activate before running any scripts:
```bash
conda activate pepxpro
```

To create from scratch:
```bash
conda create -n pepxpro -c conda-forge python=3.12 openmm matplotlib seaborn numpy pandas biopython pdbfixer scipy scikit-learn dssp pymol-open-source ipykernel jupyter
```

Then install cd-hit (bioconda) and tmtools (pip-only):
```bash
conda install -c conda-forge -c bioconda cd-hit
python -m pip install tmtools
```

Optional — only needed if using external datasets in `.xlsx` format (`04_external.py`):
```bash
conda install -c conda-forge openpyxl
```

To use the tutorial notebook, register the kernel once after creating the env:
```bash
python -m ipykernel install --user --name pepxpro --display-name "pepxpro"
```

Note: always use `python -m pip install` (not plain `pip install`) inside the conda
environment to avoid installing into the system Python.

After creating the environment, download DSSP dictionary files (one-time).
Run with the environment activated so `$CONDA_PREFIX` resolves correctly:
```bash
conda activate pepxpro
mkdir -p $CONDA_PREFIX/var/cache/libcifpp/
curl -o $CONDA_PREFIX/var/cache/libcifpp/mmcif_pdbx.dic https://mmcif.wwpdb.org/dictionaries/ascii/mmcif_pdbx_v50.dic
curl -o $CONDA_PREFIX/var/cache/libcifpp/mmcif_ma.dic https://mmcif.wwpdb.org/dictionaries/ascii/mmcif_ma.dic
curl -o $CONDA_PREFIX/var/cache/libcifpp/components.cif https://files.wwpdb.org/pub/pdb/data/monomers/components.cif
```
Note: `components.cif` is ~400 MB — the last download takes a few minutes.

## Sibling project

`PDBbindVersions/` (one level up, at `PREPPRO/PDBbindVersions/`) contains a Jupyter notebook
comparing v2020 vs v2020R1 index files. Key finding: R1 is a strict subset of v2020
(406 entries removed, 0 new, 314 changed binding data / resolution / year).

## Known gotcha — DSSP

If the libcifpp dictionary files are not downloaded, DSSP returns "NA" for all residues.
genSample then drops all entries when `--dims ss_class` is used (ss_class = "NA" → filtered).
Quick Start and Marty's setup: use `--dims dg pep_len` (does not require DSSP).
Verify DSSP setup: `ls -lh $CONDA_PREFIX/var/cache/libcifpp/`

## Commit convention

**Never add Co-Authored-By lines** to commit messages in this repo.

## Coding guidelines

- Prefer modular scripts over monolithic scripts
- Keep functions small and reusable
- Use pathlib for filesystem operations
- Use reproducible random seeds
- Keep metadata TSVs lightweight and human-readable
