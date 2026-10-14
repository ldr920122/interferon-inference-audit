# Additional file 9: version 1.7 analyses (supersedes earlier versions)
All analyses use saved project outputs or public GEO metadata; no new experimental data were generated.

| File | Purpose |
|---|---|
| stratified_validation.py → stratified_validation_per_dataset_per_unit.csv | Table 6. Per-dataset, per-unit summary (official and reproduced signatures). v1.6 fix: the v1.5 code counted only sign patterns with sum ≥ observed, i.e. a one-sided P, while the text said two-sided. Both are now reported; the manuscript uses the two-sided value. |
| patient_level_sensitivity.py → .json | GSE72754 exact permutation, detection-limit alternatives and within-patient permutation; GSE11903 TNFA exact sign-flip, exact permutation and patient bootstrap (all two-sided, unchanged from v1.5; post hoc). |
| superseries_relations_66.csv | GEO SuperSeries/SubSeries relations for the 64 GEO candidates, and related series in CytoSig final training (queried 13 Oct 2026). |
| GSE60236_gsm_level_check.csv, geo/ | GSM-level membership of SuperSeries GSE60236 (2,054 GSMs = GSE60234 + GSE60235 + GSE60341), stimulation and donor IDs. |
| eligibility_tiers.py → eligibility_tiers_66.csv | Three-tier grading (same as Additional file 2). Input: bg_candidate_eligibility_final.csv (manual review). v1.6: GSE60236 moved from Tier 2 to Tier 3. |
| locate_missing_column.py → missing_column_search.csv | Locates the 1,307th CytoSig final column (GSE13887 NO T cell 600 uM 24 h). |
| build_centroid_1307.py → centroid_1307_vs_official.csv, .log | Rebuilds the signature from 1,307 columns: 43/43 families identical to the official signature (max abs diff 4.4e-16), with SHA-256 of outputs. |
| harmonize_labels.py | v1.7: one controlled training-identity vocabulary (cytosig_final_training, cytosig_candidate_pool_only, superseries_contains_cytosig_final_training, not_in_cytosig) in the leakage index, Additional files 2 and 3; GSE60236 counts in Additional file 3 corrected (411 IFN-β, 0 IFN-γ); gsm_both_beta_and_gamma column (41 of 66); counts json recomputed. |
| make_additional_file_5.py → centroid_1306_vs_1307_vs_official.csv | Regenerates Additional file 5 (1,306 columns 42/43; 1,307 columns 43/43). |

Run with Python 3.10 (pandas, numpy, scipy, matplotlib). Paths are relative to the project tree, except CytoSig inputs under /workspace/cytosig_repro.

Order: eligibility_tiers.py, then harmonize_labels.py (which adds the controlled labels to its output).
