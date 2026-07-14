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
| `scrape/` | Download, filter, and prepare CIFs from PDBbind and external sources |
| `genSample/` | Generate, analyze, and visualize curated dataset subsets |
| `pepxpro_benchmark/` | Curated, versioned diverse benchmark releases (`v1/`, `v2/`, …) |
| `ml_benchmark/` | Notebook validating that curation choices affect ΔG prediction quality |
| `manual/` | User manual and interactive tutorial notebook |

---

## Requirements

Create the environment:

```bash
conda create -n pepxpro -c conda-forge python=3.12 \
    openmm matplotlib seaborn numpy pandas biopython pdbfixer \
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

Use the top-level runner from the project root:

```bash
python pepxpro.py scrape -v
```
This runs steps 01 → 02 → 03 → 05 from scrape module and create `outputs/summary.tsv`, the main PDBbind curated dataset used by the sampling tools. It also produces summary plots showing the distributions of the curated dataset features in `outputs/plots/To `. Use `-h` with any module to view all available parameters and options.

Alternatively, you can run one step at the time:

```bash
python pepxpro.py scrape -v --steps download
python pepxpro.py scrape -v --steps filter
python pepxpro.py scrape -v --steps summary
python pepxpro.py scrape -v --steps distributions
```

Or,

```bash
cd scrape/scripts
python 01_download.py
python 02_filter.py
python 03_summary.py
python 05_distributions.py
```

#### OPTIONAL

To use the re-processed v2020R1 index instead (better quality structures, updated binding data) add the flag to step 01:

```bash
--index scrape/inputs/INDEX_general_PL.2020R1.lst
```

To include MSE (selenomethionine) entries by converting them to MET, or to recover entries with terminal capping groups (ACE/NME/NH2) that are not at the binding interface:

```bash
python pepxpro.py scrape -v --strip-caps          # strip terminal caps not at interface
python pepxpro.py scrape -v --no-convert-mse      # disable MSE→MET (on by default)
```

To process an external-dataset and merge with the processed PDBbind dataset, use:

```bash
python pepxpro.py scrape -v --external-input scrape/inputs/external/PEPBI.xlsx
```

Where `--external-input` is the path for the desired external dataset. This creates `outputs/summary_merged.tsv`.

To view all default parameters for every pipeline stage:

```bash
python pepxpro.py config --print
```

### 2. Generate a sample (CHECK FROM HERE)

Run from genSample/scripts/. Reads from scrape/outputs/. Writes to genSample/outputs/run_XX/.

```bash
cd ../../genSample/scripts

  python 01_genSample.py --sample 20
  python 01_genSample.py --dry-run
  python 00_plot_samplingMethods.py
  python 02_plot_sampleDist.py --run run_01
  python 03_view_sample.py --run 1

```
OUTPUT (genSample/outputs/run_XX/)

------------------------------------
```bash
sample.tsv                 Selected entries.
run_config.json            All parameters used (fully reproducible).
pool_quality.tsv           Entries after quality filters.
pool_redundancy.tsv        Entries after redundancy filters.
peptide_sequences.fasta    (only with --export-fasta)
receptor_sequences.fasta   (only with --export-fasta)
cif/  or  pdb/             Structure files (only with --link-cifs / --copy-cifs / --protonate)
plots/run_comparison.png   (only with 02_plot_sampleDist.py)
```

### Optional: Explore interactively

```bash
jupyter notebook tutorial/tutorial.ipynb
```

---

## Performance reference

Approximate runtimes on first run (~478 accepted entries from PDBbind 2020 peptide subset).
Subsequent runs use the incremental cache and take seconds.

| Step | Command | Time (first run) | Notes |
|---|---|---|---|
| 01 download | `pepxpro.py scrape --steps download` | — | Depends on bandwidth (~1000 CIFs) |
| 02 filter | `pepxpro.py scrape --steps filter` | — | Incremental; re-run is instant |
| 06 tm_matrix | `pepxpro.py scrape --steps tm_matrix` | ~48 min (706 entries, 10 cores) | One-time; reused by all genSample runs |
| genSample TM (matrix) | `--tm-score-receptor 0.5` | — | With pre-computed matrix |
| genSample TM (on-the-fly) | `--tm-score-receptor 0.5` | — | Without matrix |

---

## Reproducing the PepXPro Benchmark Release

This step is mainly intended for benchmark maintainers or users who want to reproduce the versioned benchmark release exactly.
First, build the PDBbind-derived dataset and external-dataset merge in case you didn't do it before:

```bash
python pepxpro.py scrape -v 
python pepxpro.py scrape -v --external-input scrape/inputs/external/PEPBI.xlsx
```

This creates `outputs/summary_merged.tsv`, which is used as the input for the versioned benchmark release.

Then grid-search redundancy/sampling configurations for diversity, and generate the release:

```bash
cd ../../pepxpro_benchmark/scripts

python 01_optimize_benchmark.py --full           # writes optimization/results.tsv
python 02_make_benchmark.py --profile GC2        # frozen v1 release (citation-safe)
# or: python 02_make_benchmark.py --from-config auto --version v1     # re-derive the winner live
# or: python 02_make_benchmark.py --from-config GC11 --version v2     # try a specific grid config
```

---

## ML Benchmark

`ml_benchmark/notebooks/ml_benchmark.ipynb` demonstrates that PepXPro curation choices
affect ΔG prediction quality using a Random Forest regressor on 12 structural features.

**Prerequisites:** complete the scrape pipeline first (steps 1–3 of Quick start above).

**Two scenarios:**

- **Scenario A — same-model, varied test curation:** random 50/50 split; all strategies train
  on the same half; the remaining half is curated with each strategy to select N test entries.
- **Scenario B — partial-helix held out:** the partial-helix entries form a fixed test set (74);
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
