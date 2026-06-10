# PepXPro

A modular framework for scraping, generating, and sharing diverse peptide–protein benchmark datasets (structure + affinity) from PDBbind with the option of adding an external dataset (curated by the user).

PepXPro downloads and curates peptide–protein co-crystal structures, removes sequence and structural redundancy, and produces user-customizable, reproducible, diversity-aware benchmark sets ready for MD, Docking, AF-monomer benchmark, ML studies, or other potential studies.

PepXPro also provides a diverse version (v1), curated by the authors, to maximize structural, sequence, and functional diversity within the current dataset.

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
| `pepxpro_benchmark/` | Curated, versioned diverse benchmark releases (`v1/`, `v2/`, …) |
| `ml_benchmark/` | Notebook validating that curation choices affect ΔG prediction quality |
| `manual/` | User manual and interactive tutorial notebook |

---

## Requirements

Create the environment:

```bash
conda create -n pepxpro -c conda-forge python=3.12 \
    openmm matplotlib numpy pandas biopython pdbfixer \
    scipy scikit-learn dssp pymol-open-source ipykernel jupyter openpyxl 
```

Activate the environment:

```bash
conda activate pepxpro
```

Install additional dependencies:

```bash
conda install -c conda-forge -c bioconda cd-hit
python -m pip install tmtools
```

After creating the environment, download DSSP dictionary files (one-time).
Run with the environment activated so `$CONDA_PREFIX` resolves correctly:
```bash
mkdir -p $CONDA_PREFIX/var/cache/libcifpp/
curl -o $CONDA_PREFIX/var/cache/libcifpp/mmcif_pdbx.dic https://mmcif.wwpdb.org/dictionaries/ascii/mmcif_pdbx_v50.dic
curl -o $CONDA_PREFIX/var/cache/libcifpp/mmcif_ma.dic https://mmcif.wwpdb.org/dictionaries/ascii/mmcif_ma.dic
curl -o $CONDA_PREFIX/var/cache/libcifpp/components.cif https://files.wwpdb.org/pub/pdb/data/monomers/components.cif
```
Note: `components.cif` is ~400 MB — the last download takes a few minutes.
To verify the files downloaded correctly: `ls -lh $CONDA_PREFIX/var/cache/libcifpp/`

---

## Quick start

Activate the environment:

```bash
conda activate pepxpro
```

### 1. Build the dataset from PDBbind

```bash
cd scrape_pdb_bind/scripts

python 01_download.py
python 02_filter.py
python 03_summary.py
python 04_distributions.py
```

This creates `outputs/summary.tsv`, the main curated dataset used by the sampling tools.

To use the re-processed v2020R1 index instead (better quality structures, updated binding data):

```bash
python 01_download.py --index ../inputs/INDEX_general_PL.2020R1.lst
```

### 2. Generate a diverse sample

Generate a reproducible sample of 20 protein–peptide complexes. The dataset is first filtered to remove sequence and structural redundancy, then sampled using Kennard–Stone to maximize diversity in binding affinity (ΔG) and peptide length.

```bash
cd ../../gen_sample/scripts

python 01_gen_sample.py \
    --sample 20 \
    --seed 42 \
    --copy-cifs \
    --seq-id-receptor strict \
    --seq-id-peptide strict \
    --tm-score-receptor 0.5 \
    --sampling kennard-stone \
    --dims dg pep_len
```

Outputs:

- Sample summary dataset (`.tsv`)
- Selected CIF structures (`--copy-cifs`)
- Run metadata and logs

### Optional: Explore interactively

```bash
jupyter notebook tutorial/tutorial.ipynb
```

---

## Reproducing the PepXPro Benchmark Release

This step is mainly intended for benchmark maintainers or users who want to reproduce the versioned benchmark release exactly.
First, build the PDBbind-derived dataset and external-dataset merge in case you didn't do it before:

```bash
cd scrape_pdb_bind/scripts

python 01_download.py
python 02_filter.py
python 03_summary.py

python 05_external.py \
    --input ../inputs/external/PEPBI.xlsx \
    --name pepbi
```

This creates `outputs/summary_merged.tsv`, which is used as the input for the versioned benchmark release.

Then create the benchmark release:

```bash
cd ../../pepxpro_benchmark/scripts

python 01_make_benchmark.py
```

---

## ML Benchmark

`ml_benchmark/notebooks/ml_benchmark.ipynb` demonstrates that PepXPro curation choices
affect ΔG prediction quality using a Random Forest regressor on 12 structural features.

**Prerequisites:** complete the scrape_pdb_bind pipeline first (steps 1–3 of Quick start above).

**Two scenarios:**

- **Scenario A — same-model, varied test curation:** random 50/50 split; all strategies train
  on the same half; the remaining half is curated with each strategy to select N test entries.
- **Scenario B — partial-helix held out:** ~96 partial-helix entries form a fixed test set;
  each strategy curates the remaining entries for training.

**Six curation strategies per scenario:** Random, Seq-70%+KS 2D/4D/8D/12D, Seq-70%+Rec-40%+KS 12D.

**Metrics:** RMSE, MAE, Pearson r, Spearman ρ (mean ± std over 25 seeds).

Open and run with the `pepxpro` kernel:

```bash
jupyter notebook ml_benchmark/notebooks/ml_benchmark.ipynb
```

---

## Documentation

Full pipeline documentation, flag reference, and examples: [`manual/pepxpro_manual.txt`](manual/pepxpro_manual.txt)

Interactive walkthrough: [`tutorial/tutorial.ipynb`](tutorial/tutorial.ipynb)

---

## License

MIT License. See [LICENSE](LICENSE) for details.

---

## Authors

L. America Chi and F. Marty Ytreberg
