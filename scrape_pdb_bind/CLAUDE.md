# scrape_pdb_bind

Modular pipeline for downloading, filtering, and preparing peptide–protein systems from PDBbind.

## Environment

```bash
conda activate pepxpro
```

## Structure

```
scrape_pdb_bind/
├── scripts/
│   ├── 01_download.py       Parse INDEX, compute ΔG, download assembly CIFs
│   ├── 02_filter.py         Filter, fix, extract metadata
│   ├── 03_summary.py        Build summary table
│   ├── 04_distributions.py  Box plot distributions
│   └── 05_external.py       Integrate external datasets (e.g. PEPBI.xlsx)
├── modules/
│   ├── utils.py             load_tsv, save_tsv
│   ├── download.py          parse_index, parse_temperature, calc_dG, download_all
│   ├── metadata.py          extract_sequence, extract_resolution, classify_peptide_ss
│   └── filter.py            process_entry, _check_assembly, build_fixed_cif
├── outputs/
├── logs/
├── inputs/
│   ├── INDEX_general_PL.2020         PDBbind v2020 (19,443 entries)
│   ├── INDEX_general_PL.2020R1.lst   PDBbind v2020R1 (19,037 entries, strict subset of v2020)
│   └── external/                     External datasets (PEPBI.xlsx, etc.)
├── pdb_cache/               Raw assembly CIFs (<pdb_id>-<N>.cif)
└── fixed_cif/               Processed CIFs (AMBER14-compatible, missing residues filled)
```

## Running the pipeline

```bash
cd scrape_pdb_bind/scripts

python 01_download.py
python 02_filter.py
python 03_summary.py
python 04_distributions.py
```

Each script accepts `--clean` to remove its output files.

## Key flags

### 01_download.py
- `--index FILE` — PDBbind index file (default: `inputs/INDEX_general_PL.2020`)
  Use `--index ../inputs/INDEX_general_PL.2020R1.lst` for the re-processed version.
- `--temperature K` — fallback temperature in K when not in index notes (default: 300.0).
  Temperature is parsed automatically from index comment field when available.

### 02_filter.py
- `--force` — ignore cache, reprocess all entries
- `--input / --output / --log` — override default paths

### 05_external.py
- `--input FILE` — external dataset (xlsx/tsv/csv)
- `--name NAME` — short identifier for output files
- `--temperature K` — fallback temperature (default: 300.0)
- `--merge-only` — re-merge all existing external summaries without re-running pipeline
- `--no-download / --no-filter / --no-summary / --no-merge` — skip steps

## Outputs

### outputs/peptide_entries.tsv
Columns: `pdb_id`, `n_mer`, `resolution`, `year`, `binding_type`, `dG_kcal`,
         `dG_approx`, `step1_reject`, `temperature`

### outputs/filtered_systems.tsv
Columns (accepted rows):
`pdb_id`, `n_mer`, `accepted`, `assembly_num`,
`peptide_chain`, `peptide_length`, `peptide_added`, `peptide_seq`,
`peptide_ss` (deposited _struct_conf), `peptide_ss_dssp` (mkdssp),
`receptor_chains`, `receptor_length`, `receptor_added`, `receptor_seq`,
`receptor_ss`, `receptor_ss_dssp`,
`resolution`, `receptor_description`, `receptor_func_class`,
`n_contacts` (peptide heavy atoms within 4.5 Å of receptor), `n_hbonds` (N/O/S proxy at 3.5 Å),
`reject_reason`

### outputs/summary.tsv
Columns: `pdb_id`, `receptor`, `peptide`, `dG_kcal_mol`, `binding_type`,
         `resolution`, `year`, `temperature`,
         `peptide_added`, `peptide_seq`, `peptide_ss`, `peptide_ss_dssp`,
         `peptide_pct_helix`, `peptide_pct_beta`, `peptide_pct_coil`, `peptide_ss_class`,
         `receptor_added`, `receptor_seq`, `receptor_ss`, `receptor_ss_dssp`,
         `receptor_pct_helix`, `receptor_pct_beta`, `receptor_pct_coil`,
         `receptor_func_class`, `n_contacts`, `n_hbonds`

Sorted strong→weak by ΔG.

## Step details

**Step 1:** Parses INDEX for `(N-mer)` entries. Keeps only exact Kd/Ki (drops IC50, inequalities).
Parses temperature from comment field (regex); falls back to `--temperature` (default 300 K).
Computes ΔG = RT·ln(Kd/Ki). Rejects covalent complexes. Downloads all assembly CIFs
(assemblies 1, 2, 3 … until 404; 8 concurrent workers).

**Step 2:** Selects best assembly (fewest chains with a peptide candidate). Identifies peptide
(shortest protein chain ≤ max(50, n_mer×2) residues). Checks AMBER14 whitelist. Identifies
receptor chains within 4 Å of peptide. Rejects receptor gaps > 10 residues or covalent
peptide–receptor continuity. Fills missing residues (PDBFixer). Rejects CA-CA gap > 5 Å.
Extracts sequences, deposited SS (_struct_conf), DSSP SS (mkdssp), receptor description.
Computes interface contacts (`n_contacts`) and H-bond proxies (`n_hbonds`) from the fixed CIF.
Caches results — use `--force` to reprocess.
Missing lightweight columns (`receptor_description`, `receptor_func_class`, `n_contacts`, `n_hbonds`)
are patched cheaply from existing `fixed_cif/` files without full reprocessing.

**Step 3:** Joins filtered + entries tables. Auto-detects extra columns from entries TSV.
Computes ss_percentages and ss_class from DSSP string.

## Key constants in modules/filter.py (hardcoded, no flags yet)

- `NEIGHBOR_CUTOFF_A = 4.0` — receptor identification distance (Å)
- `MAX_RECEPTOR_GAP = 10` — max gap in receptor residue numbering
- `PEPTIDE_BOND_CUTOFF = 1.8` — covalent terminus detection (Å)
- `PEPTIDE_COPY_THRESHOLD = 0.7` — sequence similarity to exclude peptide copies
- `CONTACT_CUTOFF_A = 4.5` — heavy-atom contact distance for `n_contacts`
- `HBOND_CUTOFF_A = 3.5` — N/O/S proximity distance for `n_hbonds`
- peptide_max_len = `max(50, n_mer × 2)`

## ΔG formula

ΔG = R·T·ln(conc_M)  where R = 0.001987 kcal/(mol·K)
Temperature: parsed from index comment when available, else `--temperature` (default 300 K).
Per-entry temperature stored in `peptide_entries.tsv` and `summary.tsv`.
