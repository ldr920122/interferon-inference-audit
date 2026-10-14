# Frozen legacy analyses

The saved tables from v1.7 are retained. In these tables `signature=reproduced` means the 1,306-column reconstruction, not the resolved 1,307-column signature. `signature=official` means the frozen official-signature analysis. These outputs used their original permutation streams. Updated 1,307-column outputs with matched permutations are in Additional file 10; they do not overwrite these values.

Historical scripts document the original environment and can contain /workspace paths. They are not advertised as standalone reproductions. The v1.8 standalone replay script is in Additional file 10 and was executed from a fresh process. Full source filtering still requires the original large CytoSig/TCGA/GTEx inputs; it was not rerun in this revision.

The v1.8 stratified CSV restores the exact small probability for 18 technical units (2/2^18), which previously rounded to zero. The method and donor caveats are unchanged.
