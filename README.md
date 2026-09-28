# Winner-Path Complexity and Score Margins — V22 code

Computational companion for **Winner-Path Complexity and Score Margins in a
Two-Criterion Entropy-Weighted Orthopair Model** (manuscript V22).

Authors: Ghassan Malkawi, Ahmed Elsayed, Mohammed Alhagyan, Nazihah Ahmad,
Haslinda Ibrahim, and Wan Suhana Wan Daud. Manuscript corresponding author:
Mohammed Alhagyan. See [AUTHORS.md](AUTHORS.md).

## Current public code

Start at [reproducibility/V22](reproducibility/V22/). Download the standalone
[V22 reproducibility ZIP](releases/QROF_V22_Reproducibility.zip).
SHA-256: `e09ac8027edd50e33379b56be1224917c9206511eca6b887067ac1514d00a87f`; 5,159,988 bytes;
138 members including its exact SHA-256 manifest.

The module READMEs, claim maps and QA records state exactly which computations
were rerun and which results remain analytic manuscript proofs. The artifact
contains frozen derived aggregate inputs, grid comparisons, a post-hoc housing
analysis and explicitly identified numerical certificates. It is not a
complete archive of every historical calculation or independent human peer
review. Install the documented dependencies separately and run the supplied
commands with fresh external output directories. No runtime or publisher PDF
is bundled. Start integrity checking with:

```sh
cd reproducibility/V22
python -B verify_manifest.py
```

## Zenodo access status

The identifier **10.5281/zenodo.22793242** is currently reserved in a draft;
publication is not yet verified. Do not cite it as an accessible public archive
until this status has been updated after anonymous record/file verification.
GitHub V22 access is independent of Zenodo publication.

## Preserved history

[V18 code](reproducibility/V18/) and its unchanged
[V18 ZIP](releases/QROF_V18_Reproducibility.zip) remain available.
V18 ZIP SHA-256 is `208a14374e757ea264c4e77bead0c4b076c07f1aba8dc6c3c808c175cfdef6fb`.
Older titles, source-radius bounds and QA records describe their historical
versions, not the current refined V22 evidence. Prior commits are preserved;
V19–V21 private manuscript packages are not redistributed here.

## Scope, rights and review

The main law is restricted to the stated two-criterion, common-nonmembership,
entropy-weighted model. Sampled grid stability is not a continuum certificate;
constructed examples do not establish prevalence or real-world superiority.
Read the [V22 scope and claim maps](reproducibility/V22/README.md).

The existing [CC0 licence](LICENSE) applies to original material only.
Derived third-party data retain their source conditions; see the
[third-party notices](reproducibility/V22/THIRD_PARTY_NOTICES.md), including
CAA's no-onward-sale restriction and source attribution requirements.

OpenAI ChatGPT/Codex assisted literature review, mathematical development,
drafting, programming and verification. Automated replays and agent reviews
do not establish completed human scientific approval of the new V22 results.
The authors remain responsible for their final work. Repository publication
is not journal submission, acceptance or publication of the manuscript.
See [CITATION.cff](CITATION.cff) for the versioned software citation.
