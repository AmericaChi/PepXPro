# PepXPro

PepXPro is a modular framework for generating peptide–protein benchmark datasets from PDBbind.

## Repository structure

- `scrape_pdb_bind/`
  Modular pipeline for downloading, filtering, and preparing peptide–protein systems from PDBbind.
  Modules in `modules/`, scripts in `scripts/`, outputs in `outputs/`, logs in `logs/`.

- `gen_sample/`
  Modular benchmark/sample generation framework. Reads from `scrape_pdb_bind/outputs/`.

- `pepxpro_benchmark/`
  Curated benchmark releases. Scripts in `scripts/`, versioned outputs in `v1/`, `v2/`, etc.
  Uses gen_sample modules directly — no run_XX created in gen_sample/outputs/.
  Workflow: tune parameters in gen_sample first, then run `python scripts/01_make_benchmark.py`.

- `manual/`
  Comprehensive user manual: `pepxpro_manual.txt`.

## Development philosophy

- Reproducible workflows
- Modular analysis scripts
- Diversity-aware sampling
- Sequence and structural redundancy reduction
- Benchmark-oriented dataset generation
- ML-friendly dataset splits

## Implemented gen_sample modules

- Quality filtering (length, dG, resolution, binding type, added residues, exclusion list)
- Sequence redundancy filtering (CD-HIT for receptor at ≥40%, Biopython pairwise for peptide)
- Structural redundancy filtering (TM-score receptor via tmtools, parallelized)
- Stratified, Kennard-Stone, KDE, isolation-forest, cluster sampling
- Structure export (CIF/PDB, symlink/copy, protonation)
- Sampling benchmark plots, per-run comparison plots, FASTA export

## Planned gen_sample modules

- Functional redundancy filtering
- Peptide feature analysis (charge, hydrophobicity, aromaticity)
- Benchmark export with ML-friendly splits

## Pending work — ML benchmark

Design agreed but not implemented. Goal: show that PepXPro curation choices
(redundancy removal, sampling strategy, train/test split) affect ΔG prediction quality.

Axes of comparison:
1. Redundancy strategy — with/without seq-id-receptor, with/without tm-score-receptor
2. Sampling method — stratified vs kennard-stone vs kde vs cluster
3. Train/test split — random vs by-receptor vs by-deposition-year

Structural features to compute per entry (not yet in pipeline):
- Number of peptide–receptor contacts (atoms within X Å)
- Number of hydrogen bonds

iRMSD was discussed but deferred (requires reference structure).

Open question not yet decided: dataset-quality benchmark (measure leakage/distribution
properties) vs model-predictive benchmark (train a model, measure RMSE/Pearson on test set).

## Output organization

scrape_pdb_bind/outputs/
├── peptide_entries.tsv
├── filtered_systems.tsv
├── summary.tsv
└── plots/
    ├── summary_distributions.png
    └── summary_statistics.tsv

gen_sample/outputs/run_XX/
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
conda create -n pepxpro -c conda-forge python=3.12 openmm matplotlib numpy pandas biopython pdbfixer scipy scikit-learn dssp pymol-open-source ipykernel jupyter
```

Then install cd-hit (bioconda) and tmtools (pip-only):
```bash
conda install -c conda-forge -c bioconda cd-hit
python -m pip install tmtools
```

Optional — only needed if using external datasets in `.xlsx` format (`05_external.py`):
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
gen_sample then drops all entries when `--dims ss_class` is used (ss_class = "NA" → filtered).
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
