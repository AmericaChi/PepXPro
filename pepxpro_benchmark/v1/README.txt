PepPro Benchmark — v1
==========================================
Generated : 2026-06-08

CONTENTS
--------
sample.tsv                  87 entries, full metadata
cif/                        87 fixed CIFs (AMBER14-compatible)
pdb/                        87 PDB conversions
sequences/
  peptide_sequences.fasta
  receptor_sequences.fasta
plots/
  summary_distributions.png
  run_comparison_dg_ss_class.png
  run_comparison_dg_func_class.png
pool_quality.tsv            Entries after quality filters
pool_redundancy.tsv         Entries after redundancy filters
run_config.json             Full provenance

PIPELINE
--------
  Input:           /Users/amerikachi/Projects/Claude_projects/PREPPRO/PepPro/scrape/outputs/summary_merged.tsv
  Quality pool:    582 entries  (--single-chain)
  Redundancy pool: 89 entries  (seq receptor 40%, peptide 70%, TM-score 0.5)
  Final sample:    87 entries  (kennard-stone-3D)

SELECTION CRITERIA
------------------
Single receptor chain.
Sequence redundancy: receptor 40% CD-HIT, peptide 70% Biopython pairwise.
Structural redundancy: receptor TM-score < 0.5 (different fold).
Sampling: Kennard-Stone in dG x peptide-SS-class x receptor-function space.
Source: PDBbind 2020 + PEPBI (summary_merged.tsv), seed=42.

PARAMETERS
----------
{
  "sample_size": 87,
  "pool_size": 89,
  "sampling_method": "kennard-stone-3D",
  "dims": [
    "dg",
    "ss_class",
    "func_class"
  ],
  "seed": 42,
  "filters": {
    "single_chain": true,
    "seq_id_receptor": "strict",
    "seq_id_peptide": "strict",
    "tm_score_receptor": 0.5
  },
  "sources": {
    "summary_tsv": "/Users/amerikachi/Projects/Claude_projects/PREPPRO/PepPro/scrape/outputs/summary_merged.tsv",
    "fixed_cif_dir": "/Users/amerikachi/Projects/Claude_projects/PREPPRO/PepPro/scrape/fixed_cif"
  },
  "generated": "2026-06-08"
}