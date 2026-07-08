PepXPro Benchmark — v2
==========================================
Generated : 2026-06-29

CONTENTS
--------
sample.tsv                  89 entries, full metadata
cif/                        89 fixed CIFs (AMBER14-compatible)
pdb/                        89 PDB conversions
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
  Input:           /Users/amerikachi/Projects/Claude_projects/PREPPRO/PepXPro/scrape/outputs/summary_merged.tsv
  Quality pool:    740 entries  (--single-chain)
  Redundancy pool: 100 entries  (seq receptor 40%, peptide 70%, TM-score 0.5)
  Final sample:    89 entries  (kennard-stone-3D)

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
  "sample_size": 89,
  "pool_size": 100,
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
    "summary_tsv": "/Users/amerikachi/Projects/Claude_projects/PREPPRO/PepXPro/scrape/outputs/summary_merged.tsv",
    "fixed_cif_dir": "/Users/amerikachi/Projects/Claude_projects/PREPPRO/PepXPro/scrape/fixed_cif"
  },
  "generated": "2026-06-29",
  "pca_features": [
    "dg",
    "pep_len",
    "rec_len",
    "n_rec_chains",
    "pep_helix",
    "pep_beta",
    "pep_coil",
    "rec_helix",
    "rec_beta",
    "rec_coil",
    "resolution",
    "binding_kd"
  ]
}