# Tests stated in advance (written 2026-10-13 07:40 Asia/Shanghai, BEFORE any GSE72754 expression data were downloaded or scored)
Original claim (Jiang 2021 Nat Methods, Fig. 4b; Supplementary Table 3):
- "IFN1 activity in blood correlates with antibody titer upon IFNA vaccine in systemic lupus. Each dot represents a blood sample. X = titer of IFNA antibody post immunization; Y = IFN1 differential activity predicted by CytoSig. Spearman correlation with a two-sided t-test p-value (n=36)."
- Supplementary Table 3: GSE72754, IFN-K, group 'IFNA&high' (response Yes) activity -21.52 (random mean -0.20, sd 1.10); 'IFNA&low' (response No) -13.47 (random 0.45, sd 0.59).
- Original source data (Source Data Fig. 4b, Springer MOESM7) list 36 samples from 6 patient IDs over 6 to 48 months; titer 100 appears to be a floor value.

Pre-stated primary test (T1, replication of the original analysis):
- Unit: post-baseline sample (expected n=36). Y = CytoSig IFN1 activity of (post-treatment log expression minus the same patient's baseline). X = anti-IFNA titer.
- Statistic: Spearman rho with a two-sided p-value. Effect: rho with a 95% CI (Fisher z, and patient-cluster bootstrap).
- Replication criterion: rho is negative with p<0.05, using the official signature (alpha=1e4).

Pre-stated secondary tests:
- T2: the same as T1 with the reproduced signature.
- T3 (pseudo-replication check): patient-level analysis. Take per-patient medians of Y and X, then compute Spearman rho (n = number of patients); this is reported descriptively only.
- T4 (IFN1 vs IFNG separation): in post-treatment samples, is the median IFN1 change more negative than the median IFNG change? Use a paired Wilcoxon test across samples, with a patient-level version as well.
- T5: the group comparison from Supplementary Table 3 (high- vs low-titer group mean IFN1 activity), reproducing the -21.52 / -13.47 direction.
Ambiguities (to be resolved from GEO metadata without changing the tests): what counts as baseline; how titers are assigned per sample; placebo arm handling (the original used IFN-K samples only).
