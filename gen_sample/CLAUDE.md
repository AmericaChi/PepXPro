# gen_sample

Modular pipeline for filtering, sampling, and exporting peptide–protein benchmark datasets
from the outputs of `scrape_pdb_bind`.

## Structure

```
gen_sample/
├── scripts/
│   └── gen_sample.py       # entry point (argparse)
└── modules/
    ├── utils.py            # load_tsv, save_tsv, enrich_rows
    ├── filters.py          # all filter functions + apply_all_filters
    ├── sampling.py         # diverse, kennard-stone, kde inverse density
    └── export.py           # write_sample, save_run_config, link_cifs
```

## Inputs

Reads from `scrape_pdb_bind/`:
- `summary.tsv` — main input (pdb_id, receptor, peptide, dG_kcal_mol, binding_type)
- `filtered_systems.tsv` — required only for `--max-added` / `--max-receptor-added`
- `fixed_cif/` — source CIF files for `--link-cifs` / `--copy-cifs`

## Outputs

Each run writes to `gen_sample/outputs/run_XX/`:
- `sample.tsv` — selected entries
- `run_config.json` — all parameters used (filters, sampling method, seed, paths)
- `cif/` — symlinks or copies of selected CIF files (optional)

## Running

```bash
cd gen_sample/scripts

python gen_sample.py --dry-run
python gen_sample.py --sample 20
python gen_sample.py --sample 20 --sampling kennard-stone --seed 42
python gen_sample.py --sample 20 --sampling kde
python gen_sample.py --no-sample
python gen_sample.py --sample 20 --single-chain --max-peptide 15 --link-cifs
python gen_sample.py --sample 20 --min-dg -15 --max-dg -8 --binding-type kd
python gen_sample.py --sample 20 --max-added 3 --max-receptor-added 5
python gen_sample.py --sample 20 --exclude bad_ids.txt --copy-cifs
```

## Filters available

| Flag | Description |
|---|---|
| `--single-chain` | Receptor with exactly one chain |
| `--min/max-receptor N` | Receptor residue count bounds |
| `--min/max-peptide N` | Peptide residue count bounds |
| `--min/max-dg F` | Binding affinity range (kcal/mol) |
| `--binding-type kd\|ki\|both` | Measurement type |
| `--max-added N` | Max residues PDBFixer added to peptide |
| `--max-receptor-added N` | Max residues PDBFixer added to receptor |
| `--exclude FILE` | Blacklist of pdb_ids (one per line) |

## Sampling methods

| Method | Description |
|---|---|
| `diverse` (default) | 6×6 quantile-bin grid on ΔG × peptide length |
| `kennard-stone` | Maximizes minimum pairwise distance in (ΔG, peptide_length) space |
| `kde` | Inverse density sampling — favors underrepresented ΔG regions |

## Planned modules

- `structural_similarity.py` — PyMOL RMSD diversity filter (from scrape_pdb_bind step 4)
- `sequence_similarity.py` — CD-HIT / MMseqs2 redundancy removal
- `peptide_features.py` — charge, hydrophobicity, aromaticity
- `metadata.py` — resolution, experimental method, deposition date
- `secondary_structure.py` — DSSP helix/sheet fraction
