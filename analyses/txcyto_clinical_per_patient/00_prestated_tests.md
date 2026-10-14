# Prespecified tests: per-sample vs per-patient reanalysis of TxCyto clinical series (written 2026-10-13 13:30 UTC+8, before any activity was computed)
Series: GSE100093 (AMG 811, DLE skin, IFNG), GSE130588 and GSE59294 (dupilumab, AD skin, IL4 and IL13), GSE11903 (etanercept, psoriasis skin, TNFA and IL17A), GSE80060 (canakinumab, SJIA blood, IL1B), GSE99802 (fezakinumab, AD skin, IL22). GSE42048 is excluded in advance: ACHN cell-line xenografts in nude mice, no patients. GSE72754 and GSE78193 are reported separately.
Inference: official CytoSig signature.centroid, official C ridge extension (data_significance, alpha = 1e4, nrand = 1000); response = post-treatment log2 expression minus the same patient's baseline in the same tissue (lesional skin, or blood).
Tests (two-sided):
- A. Treated arm, target activity change vs 0. Per sample: Wilcoxon signed-rank over all post samples. Per patient: Wilcoxon over patient medians.
- B. Where a comparison arm exists (placebo; for GSE11903 responders vs non-responders), compare arms. Per sample: Mann-Whitney over post samples. Per patient: Mann-Whitney over patient medians.
Expected direction for blockade: target activity decreases in treated patients (more than placebo or non-responders).
A result "depends on the unit" if P < 0.05 per sample and P >= 0.05 per patient. Every series is reported, including failures to parse or run.
