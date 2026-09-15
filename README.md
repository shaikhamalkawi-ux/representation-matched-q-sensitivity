# Representation-Matched q-Sensitivity — Reproducibility Archive

Public repository: https://github.com/shaikhamalkawi-ux/representation-matched-q-sensitivity

Companion computational archive for the manuscript:

**Representation-Matched q-Sensitivity: Rank Reversal and Decision Robustness in q-Rung Orthopair Fuzzy Models**

Authors: Ghassan Malkawi, Ahmed Abdelaziz Elsayed, Mohammed Alhagyan, Nazihah Ahmad, Haslinda Ibrahim, and Wan Suhana Wan Daud.

## What this archive contains

- `benchmark/` — seven executable published q-rung orthopair fuzzy case reproductions used for the primary representation-matched benchmark.
- `analysis/` — Qin stagewise sensitivity, local sensitivity, radius-certificate replay, directed upper-witness verification, and He structural-control analysis.
- `results/` — concise public result summaries and release metadata.
- `run_all.py` — top-level verification entry point.

## Quick start

```bash
python -m pip install -r analysis/requirements.txt
python run_all.py
```

The default command runs the seven-case benchmark plus the fast Qin/He verification layer. Use `python run_all.py --radius-replay` to replay the retained radius certificate, or `python run_all.py --full` to rebuild the complete radius cover.

The analysis layer compiles an MPFR-backed C++ interval kernel locally. A C++17 compiler plus GMP/MPFR runtime/development support may be required.

## Scientific locks

The archive reproduces the manuscript's bounded conclusions. In particular:

- the primary benchmark is purposive and is not a prevalence sample;
- the matched primary paths do not change rank through the declared benchmark range;
- for the Qin application, the exact global winner-reversal radius remains **unknown / not claimed**;
- the reported bound is
  `0.000995 <= rho_win <= 0.003514796840392034`;
- the upper endpoint is a feasible witness, not a certified global minimum.

## Data and redistribution boundary

Publisher-supplied binaries, publisher PDFs, and respondent-level third-party data are not included. Source studies are identified by bibliographic metadata/DOI in the manuscript and benchmark provenance records. This public archive contains executable reconstructions, derived numeric records, verification outputs, and source code prepared for reproducibility.

No license file is included in this release candidate. Until the authors select a reuse license, copyright remains reserved by the respective rights holders.

## Citation

See `CITATION.cff`. A Zenodo DOI can be added after the GitHub repository is enabled in Zenodo and release `v1.0.0` is archived.
