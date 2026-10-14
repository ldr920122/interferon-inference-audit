# IFN independence resource v1.10

Current local manuscript candidate: v1.10, 14 October 2026 Asia/Shanghai. Not deposited, published or submitted. This package merges the full v1.8 resource with selected v1.9 additions and v1.10 fixed-input replay. Older analysis folders retain frozen outputs and historical labels; current root tables below govern the source query.

## Current materials

- query_tool/: local source lookup, current CSV/data.js and portable refresh script. Training membership is reconstructed, not an official final-column list. huCIRA non-GEO sample identities are not certified.
- eligibility_audit/: 66 candidate records; tiers 2/1/22/41, 42 with both beta and gamma arms; only the two newly admitted sources received full fixed-input replay in this revision. Legacy manual flags are not human approval.
- overlap_matrix/: current 2,080-source table, with GSE237879 counts corrected to4/4 and original28/28 retained in explicit legacy_auto fields. Other-method-role columns from v1.8 are preserved.
- cytosig_reproduction/: legacy1306 and resolved1307 signatures, source columns and original filter records. Resolved1307 matches all43 official centroids. Nothing here asserts a new exhaustive training reconstruction.
- analyses/signature_version_v1.8/: preserved 238-contrast matched-permutation verification, inputs and code; current_replay_v1.10/ contains the newly rerun24 output-block verification.
- analyses/v1.9_posthoc_frozen/: untouched historical v1.9 audit, code and outputs; absolute paths there are provenance records.
- analyses/v1.10_replay/: exact-hash source inputs, current portable replay, all43 outputs, six-seed checks, processed contrasts, clinical arithmetic verification and interpretation boundaries.
- figures/: six manuscript figures. Figure6 expresses two-channel relative order, not a complete-source accuracy estimate.
- Other analysis/search folders: retained frozen material from v1.8, including per-patient tests and failure logs.

## Interpretation limits

The two donor-paired studies are post hoc cases, not a prospectively frozen benchmark. Cell type, dose, exposure time and study pipeline differ together. Normalization sensitivity concerns log-transformed scaled counts with a fixed pseudocount; it cannot identify the cell-type cause or a general error rate across inference methods. Repeated stimulation conditions and seeds do not increase donor n.

## Reproduction and terms

Use analyses/v1.10_replay/README.md for the current replay. Older scripts may have build-machine paths; they are retained rather than silently rewritten. See LICENSE_SCOPE.md for original-code MIT, original-audit-table CC BY4.0 and third-party exclusions. Source-data redistribution terms and final author approval must be checked before external release.
