# Representation-Matched q-Sensitivity - V18 reproducibility archive

Companion code and certificate archive for **Representation-Matched q-Sensitivity: Rank Changes and Winner Robustness in Orthopair Fuzzy Models**.

Authors: Ghassan Malkawi, Ahmed Elsayed, Mohammed Alhagyan, Nazihah Ahmad, Haslinda Ibrahim, and Wan Suhana Wan Daud. Corresponding author: Mohammed Alhagyan. See [AUTHORS.md](AUTHORS.md).

## Current release

The current public code is in [reproducibility/V18](reproducibility/V18/). Start with its README and scope inventory. Earlier Git history and baseline release notes are preserved; old titles, dates and QA summaries refer to their original versions, not to V18.

Download the standalone [QROF V18 reproducibility ZIP](releases/QROF_V18_Reproducibility.zip). Its SHA-256 is `208a14374e757ea264c4e77bead0c4b076c07f1aba8dc6c3c808c175cfdef6fb` (1,775,127 bytes). It contains 116 files including the integrity manifest. The included `PUBLIC_QA.json` records the complete fresh replay, including the slower Qin retained cover; the packaged ZIP was also extracted and its Seikh core recompiled and replayed again.

After cloning this repository, run integrity checking from `reproducibility/V18`:

```sh
python -B verify_manifest.py
```

Follow that directory's README to install the required dependencies separately and run numerical replay in a new output directory. The released scripts do not install dependencies or access the network.

The reserved Zenodo identifier is **10.5281/zenodo.22793242**. It must not be treated as a published, accessible deposit until publication is verified. GitHub access does not by itself establish Zenodo publication.

## Scientific scope

The curated release includes the recoverable seven-case benchmark, Qin and He checks, Qin interval covers, Zhang's score-table control, and Seikh source-coordinate certificates and rational witnesses. The inventory distinguishes fresh replay from inherited evidence. It is not a complete archive of every historical analysis.

- Seikh uses the same 40-coordinate normalized-baseline-raw L-infinity domain, center, and metric at the compared rungs.
- The certified bounds imply `rho_src(4) <= 0.02321968 < 0.0235 <= rho_src(16)`.
- The all-finite-real-`q >= 4` uniform bracket is `0.019 <= rho_src <= 0.02321968`.
- These are bounds, not exact radii or a proof of monotonicity for all `q`. The `q=8` interval overlaps the others.
- The Qin model H, q=3, canonical L2 bound `0.000995 <= rho_win <= 0.003514796840392034` is separate; it is not an exact global minimum or a bound for every rung/model.
- The seven-case benchmark is purposive, not a prevalence sample. The historical 30-record candidate-state ledger remains unavailable.

## Public-release boundary

Publisher PDFs and binaries, respondent-level third-party data, internal development notes, credentials, and superseded internal manuscript packages are not redistributed. Source publications remain separately cited. Deterministic numerical replay is not independent human peer review or proof-assistant verification.

OpenAI ChatGPT/Codex assisted with drafting, literature review, mathematical development, and code development and checking. The authors reviewed the outputs and take responsibility for the work.

## License and citation

The existing [CC0 1.0 license](LICENSE) applies to original material only; third-party rights are not waived. See [CITATION.cff](CITATION.cff). Repository publication is not journal submission, acceptance, or publication of the manuscript.
