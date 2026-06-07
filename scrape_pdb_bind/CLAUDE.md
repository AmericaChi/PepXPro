# scrape_pdb_bind

Modular pipeline for downloading, filtering, and preparing peptide–protein systems from PDBbind.

## Environment

```bash
conda activate peppro
```

To create from scratch:
```bash
conda create -n peppro -c conda-forge python=3.12 openmm matplotlib numpy pandas biopython pdbfixer scipy scikit-learn pymol-open-source
```

Key packages: `biopython` (MMCIF2Dict), `openmm` (PDBxFile, Modeller), `pdbfixer`, `scipy` (cKDTree), `numpy`, `matplotlib`.

## Structure

```
scrape_pdb_bind/
├── scripts/
│   ├── 01_download.py       Download assembly CIFs, compute ΔG
│   ├── 02_filter.py         Filter, fix, and extract metadata
│   ├── 03_summary.py        Build summary table
│   └── 04_distributions.py  Box plot distributions of summary columns
├── modules/
│   ├── utils.py             load_tsv, save_tsv
│   ├── download.py          parse_index, download_all
│   ├── metadata.py          extract_sequence, extract_resolution, extract_ss_fractions
│   └── filter.py            process_entry, _check_assembly, build_fixed_cif
├── outputs/
│   ├── peptide_entries.tsv
│   ├── filtered_systems.tsv
│   ├── summary.tsv
│   └── plots/
│       └── summary_distributions.png
├── logs/
│   ├── 02_filter.log
│   └── 03_summary.log
├── inputs/
│   └── INDEX_general_PL.2020    PDBbind flat-file index (manual download, requires registration)
├── pdb_cache/               Raw assembly CIFs (<pdb_id>-<N>.cif)
└── fixed_cif/               Processed CIFs (AMBER14-compatible, missing residues filled)
```

## Running the Pipeline

```bash
cd scrape_pdb_bind/scripts

python -u 01_download.py 2>&1 | tee ../logs/01_download.log
python -u 02_filter.py   2>&1 | tee ../logs/02_filter.log
python -u 03_summary.py
python    04_distributions.py
```

Each script accepts `--clean` to remove its output files (never touches `pdb_cache/` or `fixed_cif/`).

## Outputs

### outputs/peptide_entries.tsv
Columns: `pdb_id`, `n_mer`, `year`, `binding_type`, `dG_kcal`, `dG_approx`, `step1_reject`

### outputs/filtered_systems.tsv
Columns (accepted rows):
`pdb_id`, `n_mer`, `accepted`, `assembly_num`,
`peptide_chain`, `peptide_length`, `peptide_added`, `peptide_seq`, `peptide_ss`,
`receptor_chains`, `receptor_length`, `receptor_added`, `receptor_seq`, `receptor_ss`,
`resolution`, `reject_reason`

- `peptide_seq` — 1-letter amino acid sequence
- `receptor_seq` — per-chain sequences: `"A:MKTAY...|B:ACDEF..."`
- `peptide_ss` / `receptor_ss` — H/E/C fractions from deposited `_struct_conf` / `_struct_sheet_range`: `"H:0.30,E:0.10,C:0.60"` or `"NA"`
- `resolution` — crystallographic resolution in Å, or `"NA"` for NMR/missing

### outputs/summary.tsv
Columns: `pdb_id`, `receptor`, `peptide`, `dG_kcal_mol`, `binding_type`, `resolution`, `year`,
`peptide_added`, `peptide_seq`, `peptide_ss`, `receptor_added`, `receptor_seq`, `receptor_ss`
Sorted strong→weak by ΔG. Full metadata — superset of all non-redundant columns from
`peptide_entries.tsv` and `filtered_systems.tsv`.

## Step Details

**Step 1:** Parses INDEX for `(N-mer)` entries. Keeps only Kd/Ki (drops IC50/NA). Computes ΔG = RT·ln(Kd) at 300 K. Rejects covalent complexes. Downloads all assembly CIFs (assemblies 1, 2, 3, ... until 404; 8 concurrent workers). Extracts deposition year from INDEX column 2.

**Step 2:** Main filtering and structure preparation.
1. Selects best assembly (fewest chains containing a peptide candidate)
2. Identifies peptide chain (shortest protein chain ≤ max(50, n_mer×2) residues)
3. Checks AMBER14 whitelist for peptide and receptor chains
4. Identifies receptor chains within 4 Å of the peptide
5. Rejects entries with receptor gaps > 10 residues or covalent continuity between peptide and receptor
6. Extracts sequences (1-letter), resolution, and secondary structure fractions from CIF
7. Builds fixed CIF: fills missing residues, removes water, retains AMBER14 ions
8. Rejects if backbone CA-CA gap > 5 Å after fixing

**Step 3:** Joins filtered_systems.tsv with peptide_entries.tsv, formats columns, sorts by ΔG.

**Step 4:** Auto-detects numeric columns in summary.tsv, generates box plot grid, saves PNG.

## Planned Improvements

- Secondary structure: replace `_struct_conf`-based fractions with DSSP per-residue assignment
- Metadata: R-factor, experimental method (X-ray/NMR/cryo-EM)
- Sequence redundancy filtering → `gen_sample/` module

## Key Conventions

- All file I/O uses `encoding='ascii'`
- ΔG: `R = 0.001987 kcal/(mol·K)`, `T = 300.0 K`
- Distances in Ångströms throughout
