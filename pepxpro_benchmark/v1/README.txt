PepXPro Benchmark — v1  (profile GC2)
==========================================
Generated : 2026-07-10

CONTENTS
--------
sample.tsv                  75 entries, full metadata (incl. external dataset fields)
cif/                        75 fixed CIFs (AMBER14-compatible)
pdb/                        75 PDB conversions
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
  Quality pool:    737 entries  (--single-chain)
  Redundancy pool: 165 entries  (seq receptor 40% [strict], peptide 70% [strict], TM-score 0.7)
  Final sample:    75 entries  (kennard-stone-4D, elbow N=75)

SELECTION CRITERIA
------------------
Single receptor chain.
Sequence redundancy: receptor 40% CD-HIT (strict), peptide 70% Biopython pairwise (strict).
Structural redundancy: receptor TM-score < 0.7 (different fold).
Sampling: kennard-stone-4D in dg x pep_len x rec_len x n_contacts space, N chosen at the D_composite-vs-N elbow for this grid config (pepxpro_benchmark/optimization/results.tsv).
Source: PDBbind 2020 + external dataset (summary_merged.tsv), seed=42.

PARAMETERS
----------
{
  "profile": "GC2",
  "command": "python 01_make_benchmark.py --profile GC2",
  "sample_size": 75,
  "pool_size": 165,
  "sampling_method": "kennard-stone-4D",
  "dims": [
    "dg",
    "pep_len",
    "rec_len",
    "n_contacts"
  ],
  "seed": 42,
  "filters": {
    "single_chain": true,
    "seq_id_receptor": "strict",
    "seq_id_peptide": "strict",
    "tm_score_receptor": 0.7
  },
  "sources": {
    "summary_tsv": "/Users/amerikachi/Projects/Claude_projects/PREPPRO/PepXPro/scrape/outputs/summary_merged.tsv",
    "fixed_cif_dir": "/Users/amerikachi/Projects/Claude_projects/PREPPRO/PepXPro/scrape/fixed_cif"
  },
  "generated": "2026-07-10",
  "pca_features": [
    "dg",
    "pep_len",
    "rec_len",
    "pep_helix",
    "pep_beta",
    "pep_coil",
    "rec_helix",
    "rec_beta",
    "rec_coil",
    "n_contacts",
    "n_hbonds",
    "sasa_burial"
  ]
}