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
  Workflow: tune parameters in genSample first, run `python scripts/01_optimize_benchmark.py --full`
  to grid-search redundancy/sampling configs (writes `optimization/results.tsv`), then
  `python scripts/02_make_benchmark.py --from-config auto` (or `--from-config GCn` for a specific
  grid config, or `--profile GC2` for a frozen citation-safe release) to generate the versioned release.

- `tutorial/`
  `tutorial.ipynb` (pipeline walkthrough) and `ml_benchmark.ipynb` (curation-vs-ΔG-quality
  demo — writes `tutorial/outputs/`).

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

Every module script is directly runnable on its own (e.g. `python genSample/scripts/01_gen_sample.py
--sample 20`) AND reachable through `pepxpro.py` as a thin dispatcher — each flag on the `pepxpro.py`
side is a 1:1 mirror of the real script's own argparse (declared with `default=None`, so unset flags
never override the real script's defaults; only what you actually type gets forwarded).

```bash
python pepxpro.py scrape -v                                            # full scrape pipeline
python pepxpro.py scrape -v --strip-caps                               # strip terminal caps not at interface
python pepxpro.py scrape -v --no-convert-mse                           # disable MSE→MET (on by default)

python pepxpro.py genSample gen --sample 20 --sampling kennard-stone --seed 42
python pepxpro.py genSample gen --dry-run
python pepxpro.py genSample methods --n 20                             # 00_plot_samplingMethods.py
python pepxpro.py genSample plot --run run_01                          # 02_plot_sampleDist.py
python pepxpro.py genSample view --run 1                               # 03_view_sample.py

python pepxpro.py benchmark optimize --full                            # 01_optimize_benchmark.py
python pepxpro.py benchmark make --from-config auto --version v1       # 02_make_benchmark.py
python pepxpro.py benchmark make --profile GC2                         # frozen citation-safe release

python pepxpro.py config --print                                       # show all default parameters
```

After `python -m pip install -e .` (once, inside the `pepxpro` env — see Environment below), the same
commands work via the bare `pepxpro` entry point, e.g. `pepxpro scrape -v`, `pepxpro genSample gen ...`.

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

## ML benchmark — tutorial/ml_benchmark.ipynb

Lives in `tutorial/` alongside the pipeline walkthrough notebook; writes to `tutorial/outputs/`
(gitignored). Figures are displayed inline (`plt.show()`), not saved to disk.

**Prerequisites:** the full scrape pipeline must have been run first:
```bash
pepxpro scrape -v --steps download
pepxpro scrape -v --steps filter   # adds n_contacts and n_hbonds via compute_interface_contacts()
pepxpro scrape -v --steps summary  # produces summary.tsv — the notebook's input
```

Goal: show that PepXPro curation choices affect ΔG prediction quality.
Uses Random Forest with 12 structural features (all from scrape pipeline).

**12 features:** `pep_len`, `rec_len`, `peptide_pct_helix/beta/coil`,
`receptor_pct_helix/beta/coil`, `n_contacts` (4.5 Å), `n_hbonds` (N/O/S 3.5 Å),
`binding_kd`, `resolution`.

**6 strategies per scenario:** Random, Seq-70%+KS-2D/4D/8D/12D, Seq-70%+Rec-40%+KS-12D.

**Scenario A — same-model, varied test curation:** random 50/50 split; all strategies train
on the same half; the other half is curated with each strategy to select N test entries.
**Scenario B — partial-helix held out:** the partial-helix entries as fixed test set
(74 after the 2026-07 SS reclassification; was ~96 under the old scheme);
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

Optional — registers the bare `pepxpro` command (see "pepxpro.py — top-level runner" above)
in the env's own `bin/`, editable so it always reflects the current source tree:
```bash
python -m pip install -e .
```
This works the same whether the env was created with `conda create` or `mamba create` — both
bundle `pip`, and editable console-script installs are a pip/setuptools mechanism with no
conda/mamba equivalent (same reason `tmtools` above is pip-only).

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
