# PepXPro

A modular framework for scraping, generating, and sharing diverse peptide–protein benchmark datasets (structure + affinity) from PDBbind with the option of adding an external dataset (curated by the user).
PepXPro downloads and curates peptide–protein co-crystal structures, removes sequence and structural redundancy, and produces user-customizable, reproducible, representative benchmark sets ready for MD, Docking, AF-monomer benchmark, ML studies, or other potential studies.
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

## Install

Create and activate the environment:

```bash
conda create -n pepxpro -c conda-forge python=3.12 \
    openmm matplotlib seaborn numpy pandas biopython pdbfixer \
    scipy scikit-learn dssp pymol-open-source ipykernel jupyter openpyxl

conda activate pepxpro
```

Install the remaining dependencies:

```bash
conda install -c conda-forge -c bioconda cd-hit
python -m pip install tmtools
```

Register the `pepxpro` command — the easiest way to run everything below without typing
`python pepxpro.py` every time (editable install, always reflects the current source tree):

```bash
python -m pip install -e .
```

**Required** — download the DSSP dictionary files (one-time, ~400 MB, takes a few minutes).
Without these, secondary-structure detection silently returns `"NA"` for every residue instead
of raising an error, which quietly corrupts anything downstream that depends on it — don't skip
this step. Run with the environment activated so `$CONDA_PREFIX` resolves correctly:

```bash
mkdir -p $CONDA_PREFIX/var/cache/libcifpp/
curl -o $CONDA_PREFIX/var/cache/libcifpp/mmcif_pdbx.dic https://mmcif.wwpdb.org/dictionaries/ascii/mmcif_pdbx_v50.dic
curl -o $CONDA_PREFIX/var/cache/libcifpp/mmcif_ma.dic https://mmcif.wwpdb.org/dictionaries/ascii/mmcif_ma.dic
curl -o $CONDA_PREFIX/var/cache/libcifpp/components.cif https://files.wwpdb.org/pub/pdb/data/monomers/components.cif
```

Verify: `ls -lh $CONDA_PREFIX/var/cache/libcifpp/` should show all three files.

---

## Quick start

The simplest possible path: build the PDBbind dataset, then generate a sample from it.

### 1. Build the dataset from PDBbind

```bash
pepxpro scrape -v
```

Downloads and curates the PDBbind peptide–protein subset, producing `scrape/outputs/summary.tsv`
— the curated dataset used by every sampling tool below — plus distribution plots in
`scrape/outputs/plots/`. First run downloads ~1000 structures and takes a while;
re-runs are incremental and fast.

### 2. Generate a sample

```bash
pepxpro genSample gen --sample 20
```

Filters and samples 20 entries from `scrape/outputs/summary.tsv`, writing to a new
`genSample/outputs/run_XX/`:

```
sample.tsv                 Selected entries.
run_config.json            All parameters used (fully reproducible).
pool_quality.tsv           Entries after quality filters.
pool_redundancy.tsv        Entries after redundancy filters.
```

### 3. (optional) Explore interactively

```bash
jupyter notebook tutorial/tutorial.ipynb
```

> Skipped `pip install -e .`? Every `pepxpro ...` command on this page also works as
> `python pepxpro.py ...` — no functional difference.

---

## Repository structure

| Directory | Description |
|---|---|
| `scrape/` | Download, filter, and prepare CIFs from PDBbind and external sources |
| `genSample/` | Generate, analyze, and visualize curated dataset subsets |
| `pepxpro_benchmark/` | Curated, versioned diverse benchmark releases (`v1/`, `v2/`, …) |
| `tutorial/` | Interactive pipeline walkthrough and ML benchmark notebooks |
| `manual/` | User manual (`pepxpro_manual.txt`) |

---

## Next steps

Everything past the Quick start example. Full flag reference and examples live in
[`manual/pepxpro_manual.txt`](manual/pepxpro_manual.txt); every command below also accepts `-h`/`--help`.

**Run individual scrape steps**, or the underlying scripts directly (equivalent to `pepxpro scrape`):

```bash
pepxpro scrape -v --steps download     # ~1000 CIFs; time depends on bandwidth
pepxpro scrape -v --steps filter       # incremental cache; re-run is instant
pepxpro scrape -v --steps summary
pepxpro scrape -v --steps distributions
# or, from scrape/scripts/: python 01_download.py / 02_filter.py / 03_summary.py / 05_distributions.py
```

Useful `scrape` flags: `--index scrape/inputs/INDEX_general_PL.2020` (use the original v2020
index instead of the default, re-processed v2020R1 — better quality structures, updated
binding data), `--strip-caps` (recover entries with terminal caps not at the binding
interface), `--no-convert-mse` (disable MSE→MET, on by default).

**Add an external dataset**, merged with the PDBbind pool into `scrape/outputs/summary_merged.tsv`:

```bash
pepxpro scrape -v --external-input scrape/inputs/external/PEPBI.xlsx
```

**Precompute the TM-score matrix** — speeds up `genSample gen --tm-score-receptor`, reused by
every future run (one-time, ~48 min for ~700 entries on 10 cores; ~2.7 MB output file):

```bash
pepxpro scrape -v --steps tm_matrix
```

**Other `genSample` subcommands:**

```bash
pepxpro genSample methods              # compare the 5 sampling methods (00_plot_samplingMethods.py)
pepxpro genSample plot --run run_01    # per-run distribution plot (02_plot_sampleDist.py)
pepxpro genSample view --run 1         # browse structures in PyMOL (03_view_sample.py)
```

**Reproduce the PepXPro Benchmark Release** — for benchmark maintainers, or exact reproduction
of a citation-safe release. External datasets must be added one at a time, PDBbind first
(never run two `--external-input` calls in parallel — they dedupe against each other's output):

```bash
pepxpro scrape -v
pepxpro scrape -v --external-input scrape/inputs/external/PPIKB/Affinity_Dataset_branch.xlsx --name PPIKB
pepxpro scrape -v --external-input scrape/inputs/external/PEPBI.xlsx --name PEPBI
```

This produces `scrape/outputs/summary_merged.tsv` (PDBbind + PPIKB + PEPBI). Precompute the
TM-score matrix over the full merged pool next — without it, `benchmark optimize --full` falls
back to computing TM-scores on the fly for all 27 grid configs, which can take hours instead of minutes:

```bash
pepxpro scrape -v --steps tm_matrix
```

Then grid-search and generate the release:

```bash
pepxpro benchmark optimize --full                             # writes optimization/results.tsv
pepxpro benchmark make --profile GC2                           # frozen v1 release (citation-safe)
# or: pepxpro benchmark make --from-config auto --version v1   # re-derive the winner live
# or: pepxpro benchmark make --from-config GC11 --version v2   # try a specific grid config
```

**ML benchmark** — `tutorial/ml_benchmark.ipynb` demonstrates that PepXPro curation choices
affect ΔG prediction quality using a Random Forest regressor on 12 structural features.
Requires the scrape pipeline to have been run first (Quick start, step 1). Open with the
`pepxpro` Jupyter kernel.

**Show every default parameter**: `pepxpro config --print`

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
