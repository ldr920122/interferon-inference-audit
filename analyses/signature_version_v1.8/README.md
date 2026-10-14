# Fixed-permutation signature-version verification v1.8

Post hoc comparison of 238 existing contrasts, not new independent validation.

Run `python replay_signature_check.py --out replay_outputs` from any directory after installing the recorded requirements. The replay checks all 24 dataset-by-signature output blocks against the archived values with 1e-10 tolerances. It does not download data, train a model or modify the reference outputs.

`contrast_metadata.csv` records the exact dataset order and contrast order. `*_expression_contrasts.tsv.gz` contain already-constructed stimulated-minus-control log-expression contrasts, aligned to signature genes. `download_manifest.json` records the original GEO and official CytoSig input files with URLs and SHA-256. `build_from_downloads.py` records the original contrast construction and comparison; run it with --input-dir PATH --resource-dir PATH --output-dir PATH after obtaining the files in download_manifest.json. The resource-dir is the unpacked full resource archive; the replay itself requires no downloads.

Matrices: legacy1306 is the earlier reconstructed matrix; resolved1307 adds the missing NO column; official is fetched from data2intelligence/CytoSig commit ce53cc7f. CytoSig-derived matrices are supplied under upstream GNU GPL v3 terms, not MIT. See matrices/UPSTREAM_LICENSE.txt, GPL-3.0-only.txt and the package licence scope. The newly written replay code is MIT; original comparison tables are CC BY 4.0. Public processed expression contrasts retain source-study terms and are not relicensed as original audit tables.

Main analyses remain the previously frozen official-signature outputs. Existing CSV `signature=reproduced` labels in older additional files mean legacy1306; no old result is silently relabelled as resolved1307. The new check uses the same gene permutation stream across matrices, not the legacy run's different permutation stream.
