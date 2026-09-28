# V22 empirical reproducibility module

This portable, curated module reproduces the existing six official-data domains
and eleven cohort/encoding analyses. It is not a new analysis, raw-source
re-extraction, validation sample or public copy of the private working archive.
The unchanged mathematical kernels are extracted from the frozen R15 code;
new wrappers replace workspace-specific input/admission/output routing.
`PROVENANCE.json` pins original kernels, exact ordered Fraction tables and
archived numerical references. No original input or definition was repaired.

## Run offline

Use Python 3.11 or later (tested with Python 3.12.14) with NumPy already available. Nothing installs packages,
downloads files, accesses accounts or modifies the module. From **any** working
directory, supply an absolute script path and a new output directory outside it:

```
python /path/to/empirical/run_empirical.py --out /new/results/grid
python /path/to/empirical/verify_decimal.py --out /new/results/decimal
python /path/to/empirical/verify_housing.py --out /new/results/housing
```

On Windows replace those paths with quoted absolute Windows paths. Existing output
directories are refused. Full-grid scores are generated at run time, not embedded
as a large duplicate payload. The grid job normally takes minutes; the ordinary
90-digit seven-node independent replay may take several minutes. The 60-digit
housing replay covers all 961 nodes. A failure is recorded, never converted into
a success; all three core commands plus the printed-example command below must finish successfully for this module's checks.
Python optimized mode (`-O`/`-OO`) is refused because assertions are part of the
integrity and numerical checks.

`run_empirical.py` evaluates 961 nodes for all three operators on all eleven tables,
compares every frozen summary and seven-node reference, and checks historical WA
pair diagnostics separately. `verify_decimal.py` uses separate Decimal arithmetic,
not the primary NumPy kernel, for seven nodes, all rows, three methods and three
weight arms. `verify_housing.py` independently recomputes all housing weights,
scores, ranks, plotted fields and crossing brackets. Its first argument is no
longer a private workspace. All checks verify `MANIFEST.json`.

This distribution contains compact reference values, not the old large full-grid
arrays. During release preparation the extracted module was also compared against
every original full-grid array; the release's external replay receipt states the
actual outcome. Merely reading this README does not establish that verification.

## Read before interpretation

Read `SCOPE.md`, `SOURCE_CITATIONS.md`, `THIRD_PARTY_DATA_NOTICE.md`, and
`CLAIM_MAP.json`. Exact canonical inputs are `data/*.json`; all selection and
exclusion decisions are in `exclusions/`. `expected/` contains preserved numerical
checks, not newly acquired observations. No participant-level or patient-level
records are included. Display names refer to public aggregate reporting units.

City continuum reproduction requires the separate directed-arithmetic module
if present in the final release; consult its actual coverage record. No city
continuum claim follows from this module's grid alone. This module is only one
component of V22 reproducibility.

Original program code is released under CC0 as specified by the repository.
Third-party-derived data are **not** blanket CC0: consult the rights notice,
particularly CAA's no-onward-sale restriction and required attribution.

## Optional figure regeneration and tested environment

The numerical scripts require only NumPy beyond the Python standard library.
`requirements.txt` pins the tested NumPy version; `requirements-renderer.txt`
adds the tested Matplotlib version for the optional figure renderer only.
Users may prepare a separate environment using their normal dependency workflow;
the replay programs themselves never install or access the network.

```
python /path/to/empirical/render_housing.py --out /new/results/figure
```

This reads frozen plotted values and executes the unchanged V22 R2 plotting
block; only public input/output routing and integrity gates are new. The figure
is an explicitly post-hoc illustration. The manuscript defines the winning gap
as North Dakota's score minus the maximum other-state score. No mathematical
calculation, new dataset or scientific claim is added by regeneration.

## Published operator example audit

```
python /path/to/empirical/verify_published_examples.py --out /new/results/published
```

This reruns the unchanged standard-library Decimal90 arithmetic on literal
already-transcribed source-example numeric tables and weights, then compares
every original output field except the timestamp. It preserves **29 comparisons,
28 compatible, one incompatible**, including the unresolved printed WG first
score and resulting top-order discrepancy. No source text/PDF is shipped or
reauthenticated. Liu/Wang formulas/examples were read in author-posted full text;
Khan's cosine example was read in the published PDF. See source citations and
rounding/access qualifications in `SCOPE.md` and the expected output itself.
An earlier mpmath-based checker failed before arithmetic because that library
was unavailable. The successful standard-library Decimal implementation was
separate, not a successful execution of the failed mpmath script.

## Wrapper revision history

Early public-wrapper Decimal runs were intentionally stopped after identifying
a repeated reference-NPZ decompression inside scalar comparisons. R4 caches
that unchanged array once. This changes only I/O performance, not arithmetic,
data, tolerances or reference values. Early incomplete runs are not counted
as successful numerical replays; final actual-archive replay is separately recorded.
