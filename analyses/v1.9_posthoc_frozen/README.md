# v1.9 post hoc re-audit and donor-paired rescoring (14 October 2026, Asia/Shanghai)
Post hoc: performed after the frozen analyses; not part of the prespecified validation.
- gsm_reaudit_17.csv: GSM-level recheck of candidates whose classification rested on study-design text. Plus GSE236156 (trigger). Outcome: GSE236156 and GSE237879 -> Tier 1; GSE188847 reason corrected; others unchanged.
- elig/update_eligibility_v19.py: exact edits applied to eligibility_audit/*.csv.
- GSE236156/: rescore.py (author featureCounts, hawn-lab/AM_MDM_TB_public@0076bde4 MDM_IFN/data_raw; download not redistributed, see input_sha256.json), per-donor CytoSig outputs, all 43 channels, seed check.
- GSE237879/: same contract on GEO supplementary Salmon counts (download from GEO; SHA-256 recorded).
Contract: log2(1+count/libsize*scale), scale 1e6 (prespecified) or 1e4; stimulated minus same-donor control; official CytoSig signature.centroid; ridge alpha 1e4; 1000 gene-row permutations; seed 20261014.
Key counts (correct type favoured; donors): GSE236156 IFN-g 5/5 (CPM), 5/5 (CP10k); IFN-b 1/5 (CPM), 5/5 (CP10k). GSE237879 IFN-b 4/4, 4/4; IFN-g 3/4, 4/4.
