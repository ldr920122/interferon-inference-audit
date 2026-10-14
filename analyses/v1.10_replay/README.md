# v1.10 fixed-contract replay and numerical checks

This is a post hoc reproducibility check, not new independent validation or a tuned model. All six fixed inference inputs are bundled in inputs/ and identified in replay_input_acquisition.json. Their hashes match the v1.9 records, including the HGNC snapshot. Source data retain original terms.

Run from any working directory:

```sh
python replay_posthoc.py --inputs inputs --out replay_outputs
python recheck_clinical.py
```

Dependencies: numpy, pandas; scipy is additionally required for clinical statistics. Actual versions and alpha=10000, nrand=1000, primary seed20261014 and alternative seeds1–5 are recorded in outputs/replay_contract_and_checks.json. No parameter or threshold was tuned. Paths in the preserved v1.9 snapshot are historical; the commands above are the current portable entry points.

The main manuscript uses only IFN-beta/gamma contrasts from 9 donors (18 correlated donor–stimulus contrasts). Other IFN arms were already scored in v1.9 and are retained for provenance; their presence does not expand the paper’s validation claim.

The 72 source×scale×donor×stimulus primary summaries (including the other previously scored IFN arms) reproduce v1.9 IFN-channel values within 2.6e-13. GSE236156 also has archived all43 values for comparison. GSE237879 all43 values are reconstructed and now supplied. All43 outputs and ridge coefficients for all six seeds are in outputs/.

Current clinical checks recompute tests from frozen per-sample activity outputs, not raw-expression inference. The signature-version replay rechecks the 238 previously processed contrasts in the separate signature_version_v1.8/current_replay_v1.10 folder. Raw RNA-seq QC, the complete original training filter and every source eligibility decision were not rerun.

A bootstrap-definition correction is documented in clinical_statistics/clinical_statistics_verification.json: the archived lupus interval is a patient-cluster bootstrap for the per-sample correlation, not an interval for the correlation of six patient medians. Its saved numeric value is retained.
