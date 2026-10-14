# GSE78193 case: pre-stated tests (written before any inference, 2026-10-13 Asia/Shanghai)
Context: GSE78193 contributes 2 IFNG columns (ex vivo IFN-γ stimulation of normal donor whole blood, 24 h and 48 h) to CytoSig final training; the SLE AMG 811 (anti-IFN-γ) treatment arms of the same series are listed by TxCyto SI3 as IFNG clinical validation targets.
Question: does removing the series' own training columns (leave-series-out signature) change the inferred IFNG response to AMG 811 blockade?
- T1: per SLE subject, post-dose minus own baseline log2 expression; infer IFNG activity with (a) full reproduced signature, (b) signature rebuilt without GSE78193 columns. Report per-arm medians and the Spearman correlation between (a) and (b).
- T2: AMG 811 vs placebo difference in IFNG change (Mann-Whitney, subject-level), under (a) and (b).
- T3: rank of IFNG (expect most negative) among cytokines for the arm-mean profile, under (a) and (b).
Inference: same closed-form ridge as leakage_pilot/loao_pilot.py (alpha 1e4, z-scored), no permutation z.
