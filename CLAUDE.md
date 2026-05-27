# PepPro

PepPro is a modular framework for generating peptide–protein benchmark datasets from PDBbind.

## Repository structure

- `scrape_pdb_bind/`
  Original pipeline for scraping, filtering, and preparing peptide–protein systems from PDBbind. Do not touch for now.

- `gen_sample/`
  Modular benchmark/sample generation framework. Currently working on this.

- `pep-pro-benchmark/`
  Curated benchmark datasets and metadata releases. Will work on these in the future.

## Development philosophy

- Reproducible workflows
- Modular analysis scripts
- Diversity-aware sampling
- Sequence and structural redundancy reduction
- Benchmark-oriented dataset generation
- ML-friendly dataset splits

## Planned gen_sample modules

- simple filtering 
- metadata extraction
- peptide feature analysis
- sequence similarity
- structural similarity
- secondary structure analysis
- stratified sampling
- benchmark export

## Output organization

Outputs should be organized as:

outputs/
└── run_XX/
    ├── summaries/
    ├── datasets/
    ├── logs/
    ├── cif/
    ├── pdb/
    ├── reports/
    └── run_config.json

Large generated files should not be tracked by git.

## Coding guidelines

- Prefer modular scripts over monolithic scripts
- Keep functions small and reusable
- Use pathlib for filesystem operations
- Use reproducible random seeds
- Keep metadata TSVs lightweight and human-readable
