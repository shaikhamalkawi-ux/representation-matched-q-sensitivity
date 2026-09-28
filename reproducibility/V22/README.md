# V22 reproducibility code and certificates

Companion computational material for **Winner-Path Complexity and Score Margins
in a Two-Criterion Entropy-Weighted Orthopair Model** (manuscript V22).

Authors: Ghassan Malkawi, Ahmed Elsayed, Mohammed Alhagyan, Nazihah Ahmad,
Haslinda Ibrahim, and Wan Suhana Wan Daud. The manuscript corresponding author
is Mohammed Alhagyan. This software release is not journal submission,
acceptance, peer review, or publication of the manuscript.

## Start here

1. Run `python -B verify_manifest.py` from this directory.
2. Read `empirical/README.md` for the eleven frozen encodings of six official-data
   domains, the WA/WG/cosine-TOPSIS grid comparison, and the V22 housing analysis.
3. Read `theory/README.md` for the constructed full-box witness, the refined
   source-radius certificates and exact upper witnesses, and any separately
   identified continuum controls.
4. Read the module claim maps and QA records before interpreting a replay.

Dependencies are installed separately. No dependency binaries, publisher PDFs,
private manuscript packages, credentials, internal correspondence or raw
respondent-level data are included. Execution should use a fresh output
directory; frozen input and certificate files must not be overwritten.

## What this release does and does not establish

This is a curated computational companion, not every historical research file.
Analytical proofs of the joint capacity/resource law and the restricted Frank
corollary belong to the manuscript. Numerical replay is not a substitute for
their mathematical proof or an independent human specialist review.

The central model is restricted: two benefit criteria, common raw
nonmembership, and the stated entropy-derived weights. The application grids
are analyst-defined encodings of official aggregate data, not native fuzzy
elicitation, a random prevalence sample, or a decision-quality ground truth.
Eleven encodings are not eleven independent empirical trials. Grid stability
is not a continuum certificate; explicitly certified controls are identified
separately. The housing mechanism analysis is post hoc and is not a new
held-out example. No superiority, novelty-priority or revolutionary-discovery
claim follows from publishing this archive.

Older V18 material remains separately available in the repository's
`reproducibility/V18` directory and its unchanged ZIP. Its historical bounds,
title, dates and QA records describe V18, not this release. In particular,
V22 uses the refined source-radius bounds described in its own module.

## Rights and attribution

The existing `LICENSE` is preserved. CC0-1.0 applies only to original material;
it does **not** relicense the third-party data extracts or waive their source
conditions. Read `THIRD_PARTY_NOTICES.md` and the empirical source/provenance
map. Cite the data producers and retain their attribution and reuse conditions
when redistributing derived inputs. Source papers remain separately cited and
are not distributed here.

## AI assistance and review status

OpenAI ChatGPT/Codex assisted with literature review, mathematical development,
drafting, programming, and verification. Automated replays and independent
agent checks are identified as such; they do not establish completed human
author review or peer review of the new V22 results. The authors remain
responsible for the work and its final scientific approval.

## Citation and access

See `CITATION.cff` for the software authors and version. The repository is
https://github.com/shaikhamalkawi-ux/representation-matched-q-sensitivity .
The associated Zenodo identifier is `10.5281/zenodo.22793242`; its live
publication/access status is reported on the repository landing page.
A reserved identifier alone is not evidence of public accessibility.
