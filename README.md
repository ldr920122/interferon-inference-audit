# Interferon inference under audit

Code and resource package accompanying *Interferon inference under audit: training-source overlap and patient-level validation of transcriptomic signatures* by Xiaohua Sun and Dongrui Liu.

This first public release, v1.11, accompanies manuscript v1.11. The scientific files are the complete 211-file v1.10 resource, with release documentation and third-party licence notices added. Numerical results, frozen inputs, tables, figures and historical scripts are retained. The manuscript and internal author-review files are not included in this software deposit.

## Download and cite

- [Fixed release](https://github.com/ldr920122/interferon-inference-audit/releases/tag/v1.11) includes a standalone resource ZIP and SHA-256 checksum.
- The corresponding Zenodo record is created from this tagged GitHub source snapshot. Its version DOI will be recorded in the GitHub release and manuscript after publication.
- Authors: Xiaohua Sun, Department of General Practice, and Dongrui Liu, Department of Pharmacy, The Affiliated Taizhou People's Hospital of Nanjing Medical University, Taizhou, Jiangsu, China.

## Contents

- `query_tool/index.html`: offline lookup of reconstructed training-source roles. Open in a browser; no server is needed.
- `overlap_matrix/`: 2,080-source index covering 43 cytokines, including explicit legacy fields.
- `eligibility_audit/`: 66 IFN-beta/gamma candidate records and eligibility tiers.
- `cytosig_reproduction/`: source-column records and legacy/resolved reconstructed signatures.
- `analyses/signature_version_v1.8/`: archived processed contrasts, three signatures and matched-permutation replay.
- `analyses/v1.10_replay/`: frozen inputs, portable donor-paired replay, all 43 channels, alternative seeds and clinical-statistics check.
- Other `analyses/` folders: frozen primary and historical analyses, including negative results and failure records.
- `figures/`: six manuscript figures and their available plotting scripts.
- `provenance/`: original local v1.10 documentation and original file hashes. Historical statements about an unpublished local candidate describe that earlier snapshot.

## Reproduce the current checks

Python 3.12 is required. The donor and signature replays were verified with numpy 2.3.5 and pandas 2.2.3. The clinical-statistics check additionally uses scipy; its separately recorded verification environment used numpy 2.5.4, pandas 3.0.6 and scipy 1.18.1. These are two recorded environments, not a claim that every historical script runs in one environment.

```sh
python verify_release.py
python analyses/v1.10_replay/replay_posthoc.py --inputs analyses/v1.10_replay/inputs --out replay_outputs
python analyses/signature_version_v1.8/replay_signature_check.py --out signature_replay_outputs
python analyses/v1.10_replay/recheck_clinical.py
```

Use `--seeds 20261014` for only the primary donor replay; the default also runs seeds 1 through 5. Current replay scripts use paths relative to the script or supplied arguments. Older scripts preserve build-machine paths as historical provenance. Full raw-sequence processing and the original complete training filter are outside this release's replay checks. See the analysis READMEs and JSON records for the exact scope.

## Interpretation

The two donor-paired studies are post hoc cases involving nine donors. Eighteen IFN-beta/gamma contrasts are correlated observations from those donors. Two-channel relative family ordering is distinct from all-43-source classification accuracy. Cell type, dose, exposure time and study pipeline differ together. This resource does not certify independent clinical validation, a general error rate or a uniquely identified ligand source.

## Licences and attribution

This is a multi-licence research resource. Independently written verification code uses MIT and original audit tables use CC BY 4.0. CytoSig-derived matrices and the MDM source repository retain GPL v3; HGNC data use CC0. GEO data and metadata retain their source-study terms. Consult [LICENSE_SCOPE.md](LICENSE_SCOPE.md) and [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md) before reuse. No blanket MIT licence applies to the archive.

AI tools assisted code, analysis verification and documentation preparation. The authors are responsible for checking the released content. This resource release is separate from journal submission.
