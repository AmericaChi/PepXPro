# PepPro

PepPro is a modular framework for generating peptide–protein benchmark datasets from PDBbind.

## Repository structure

- `scrape_pdb_bind/`
  Modular pipeline for downloading, filtering, and preparing peptide–protein systems from PDBbind.
  Modules in `modules/`, scripts in `scripts/`, outputs in `outputs/`, logs in `logs/`.

- `gen_sample/`
  Modular benchmark/sample generation framework. Reads from `scrape_pdb_bind/outputs/`.

- `pep_pro_benchmark/`
  Curated benchmark releases. Scripts in `scripts/`, versioned outputs in `v1/`, `v2/`, etc.
  Uses gen_sample modules directly — no run_XX created in gen_sample/outputs/.
  Workflow: tune parameters in gen_sample first, then run `python scripts/01_make_benchmark.py`.

- `manual/`
  Comprehensive user manual: `peppro_manual.txt`.

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
conda activate peppro
```

To create from scratch:
```bash
conda create -n peppro -c conda-forge python=3.12 openmm matplotlib numpy pandas biopython pdbfixer scipy scikit-learn dssp pymol-open-source ipykernel jupyter cd-hit
```

Then install pip-only packages (not available on conda):
```bash
python -m pip install tmtools
```

To use the tutorial notebook, register the kernel once after creating the env:
```bash
python -m ipykernel install --user --name peppro --display-name "peppro"
```

Note: always use `python -m pip install` (not plain `pip install`) inside the conda
environment to avoid installing into the system Python.

After creating the environment, download DSSP dictionary files (one-time):
```bash
mkdir -p /opt/miniconda3/envs/peppro/var/cache/libcifpp/
curl -o /opt/miniconda3/envs/peppro/var/cache/libcifpp/mmcif_pdbx.dic https://mmcif.wwpdb.org/dictionaries/ascii/mmcif_pdbx_v50.dic
curl -o /opt/miniconda3/envs/peppro/var/cache/libcifpp/mmcif_ma.dic https://mmcif.wwpdb.org/dictionaries/ascii/mmcif_ma.dic
curl -o /opt/miniconda3/envs/peppro/var/cache/libcifpp/components.cif https://files.wwpdb.org/pub/pdb/data/monomers/components.cif
```

## Coding guidelines

- Prefer modular scripts over monolithic scripts
- Keep functions small and reusable
- Use pathlib for filesystem operations
- Use reproducible random seeds
- Keep metadata TSVs lightweight and human-readable
