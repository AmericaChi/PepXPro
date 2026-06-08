# PepPro

A modular framework for scraping, generating, and sharing peptide–protein benchmark datasets (structure + affinity) from PDBbind or external datasets.

PepPro downloads and curates peptide–protein co-crystal structures, removes sequence and structural redundancy, and produces reproducible, diversity-aware benchmark sets ready for computational studies.

---

## Features

- Downloads and filters peptide–protein systems from PDBbind and external affinity datasets
- Removes sequence redundancy (CD-HIT for receptors, Biopython pairwise for peptides)
- Removes structural redundancy (TM-score receptor clustering, parallelized)
- Five sampling methods: stratified, Kennard-Stone, KDE, isolation forest, cluster
- Exports CIF and PDB structures, FASTA sequences, and full metadata TSVs
- Versioned benchmark releases with full provenance and reproducibility

---

## Repository structure

| Directory | Description |
|---|---|
| `scrape_pdb_bind/` | Download, filter, and prepare CIFs from PDBbind and external sources |
| `gen_sample/` | Filter, deduplicate, sample, and export dataset subsets |
| `pep_pro_benchmark/` | Curated, versioned benchmark releases (`v1/`, `v2/`, …) |
| `manual/` | User manual and interactive tutorial notebook |

---

## Requirements

```bash
conda create -n peppro -c conda-forge python=3.12 \
    openmm matplotlib numpy pandas biopython pdbfixer \
    scipy scikit-learn dssp pymol-open-source ipykernel jupyter

conda activate peppro
conda install -c conda-forge -c bioconda cd-hit
python -m pip install tmtools
```

If you plan to use external datasets in `.xlsx` format (`00_external.py`):
```bash
conda install -c conda-forge openpyxl
```

---

## Quick start

```bash
conda activate peppro

# 1. Build the dataset from PDBbind (one-time)
cd scrape_pdb_bind/scripts
python -u 01_download.py 2>&1 | tee ../logs/01_download.log
python -u 02_filter.py   2>&1 | tee ../logs/02_filter.log
python    03_summary.py
python    04_distributions.py

# 2. Generate a sample
cd ../../gen_sample/scripts
python 01_gen_sample.py --sample 50 --sampling kennard-stone \
    --seq-id-receptor strict --seq-id-peptide strict \
    --tm-score-receptor 0.5 --dims dg ss_class func_class \
    --export-fasta --copy-cifs --seed 42

# 3. Generate a versioned benchmark release
cd ../../pep_pro_benchmark/scripts
python 01_make_benchmark.py
```

Or explore interactively with the tutorial notebook:

```bash
jupyter notebook manual/tutorial.ipynb
```

---

## Documentation

Full pipeline documentation, flag reference, and examples: [`manual/peppro_manual.txt`](manual/peppro_manual.txt)

Interactive walkthrough: [`manual/tutorial.ipynb`](manual/tutorial.ipynb)

---

## License

MIT License. See [LICENSE](LICENSE) for details.

---

## Authors

L. America Chi and F. Marty Ytreberg
