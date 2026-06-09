# PepPro

A modular framework for scraping, generating, and sharing diverse peptide–protein benchmark datasets (structure + affinity) from PDBbind with the option of adding an external dataset (curated by the user).

PepPro downloads and curates peptide–protein co-crystal structures, removes sequence and structural redundancy, and produces user-customizable, reproducible, diversity-aware benchmark sets ready for MD, Docking, AF-monomer benchmark, ML studies, or other potential studies.

PepPro also provides a diverse version (v1), curated by the authors, to maximize structural, sequence, and functional diversity within the current dataset.

---

## Features

- Downloads and filters peptide–protein systems from PDBbind and external affinity datasets
- Removes sequence redundancy (CD-HIT for receptors, Biopython pairwise for peptides)
- Removes structural redundancy (TM-score receptor clustering, parallelized)
- Five sampling methods: stratified, Kennard-Stone, KDE, isolation forest, cluster
- Exports CIF and PDB structures, FASTA sequences, and full metadata TSVs
- Versioned diverse benchmark releases with full provenance and reproducibility

---

## Repository structure

| Directory | Description |
|---|---|
| `scrape_pdb_bind/` | Download, filter, and prepare CIFs from PDBbind and external sources |
| `gen_sample/` | Filter, deduplicate, sample, and export dataset subsets |
| `pep_pro_benchmark/` | Curated, versioned diverse benchmark releases (`v1/`, `v2/`, …) |
| `manual/` | User manual and interactive tutorial notebook |

---

## Requirements

Create the environment:

```bash
conda create -n peppro -c conda-forge python=3.12 \
    openmm matplotlib numpy pandas biopython pdbfixer \
    scipy scikit-learn dssp pymol-open-source ipykernel jupyter openpyxl
```

Activate the environment:

```bash
conda activate peppro
```

Install additional dependencies:

```bash
conda install -c conda-forge -c bioconda cd-hit
python -m pip install tmtools
```

After creating the environment, download DSSP dictionary files (one-time).
Run with the environment activated so `$CONDA_PREFIX` resolves correctly:
```bash
conda activate peppro
mkdir -p $CONDA_PREFIX/var/cache/libcifpp/
curl -o $CONDA_PREFIX/var/cache/libcifpp/mmcif_pdbx.dic https://mmcif.wwpdb.org/dictionaries/ascii/mmcif_pdbx_v50.dic
curl -o $CONDA_PREFIX/var/cache/libcifpp/mmcif_ma.dic https://mmcif.wwpdb.org/dictionaries/ascii/mmcif_ma.dic
curl -o $CONDA_PREFIX/var/cache/libcifpp/components.cif https://files.wwpdb.org/pub/pdb/data/monomers/components.cif
```
Note: `components.cif` is ~400 MB — the last download takes a few minutes.

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
python 01_gen_sample.py --sample 20 --sampling kennard-stone \
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
