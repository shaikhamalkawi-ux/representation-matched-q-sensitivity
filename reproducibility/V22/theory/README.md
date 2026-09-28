# V22 theory and numerical-certificate support

This module supports **specified numerical claims of manuscript V22**. It is
not the manuscript, a proof assistant, a replay of every historical calculation,
or an assertion that scientific author/specialist review is complete.
See `CLAIM_MAP.md` before interpreting a successful run.

The original mathematical programs and exact JSON inputs are byte-preserved in
`original/`. The dated relative directory hierarchy is retained so that their
existing imports and SHA-256 guards work without editing scientific code. In
particular, an upstream file bearing V17/V18 in its name is a reused kernel,
not a statement that the older V18 radius bounds are the current V22 bounds.
`THEORY_MANIFEST.json` lists every included file and hash. No reference PDFs,
private manuscript PDFs, credentials, dependency binaries or private ZIPs belong
to this module.

## Run from any directory

Use Python 3.12 with the packages in `requirements.txt` already available.
The runner never installs or downloads anything. Pure Decimal and integer
arithmetic are in Python; MPFR is supplied by `gmpy2`, and the fixed-q gradient
diagnostic also uses `mpmath`. These third-party libraries are dependencies, not
bundled code or binaries.

```text
python /path/to/theory/run_theory.py --profile quick --out /new/external/quick
python /path/to/theory/run_theory.py --profile full --out /new/external/full
```

Both output directories must be new and outside the package. Two jobs run at
most concurrently; default wall limit is 600 seconds per job, configurable only
within 1..600 seconds. A timeout/failure is retained and never silently retried
or accepted as a negative scientific result. The full profile runs eighteen jobs:

| Job family | Quick | Full | What is actually checked |
|---|---:|---:|---|
| Six-decimal 16-row box, directed MPFR256 | yes | yes | Literal transported and canonical two-product scores at four integer nodes; all 64 coordinates independently variable |
| Same box, independent outward Decimal100 | yes | yes | Separate model implementation; 256 contained score enclosures and positive all-rival margins |
| Decimal100 under Python `-OO` | no | yes | Identical result without assertions/docstrings; not an extra theorem |
| Refined Seikh fixed-q integer224 | yes | yes | Three covers, independent ordinary 100-digit derivative diagnostic, malformed-record rejection |
| Refined Seikh fixed-q MPFR256 | yes | yes | Three separately recomputed directed covers (11, 3, 3 nodes) |
| Refined Seikh uniform integer224 | no | yes | Complete 811-node cover and malformed-record tests |
| Refined Seikh uniform MPFR256 | no | yes | Independent directed evaluation of the same complete 811-node cover |
| Exact Seikh point witnesses | yes | yes | Three rational feasible upper witnesses, unique Y1, and 11 malformed-record tests |
| QOL Decimal80 and MPFR256 fixed-source covers | no | yes | Unchanged independent evaluators; complete 32/61-leaf covers, all 82 rivals |
| S4 fragile seven-row example | no | yes | Directed MPFR256, independent Decimal100 and exact rational domain/auxiliary exclusion checks; seven fixed nodes and all 28 source coordinates |
| S4 rounding protocol | no | yes | All 6–15-place trials plus 16-place baseline; failed 6–9-place tests retained, not suppressed |
| S6 cached ordinary diagnostic | no | yes | 252 literal Decimal-tuple probes and all 18 frozen cases; unchanged 180-bisection body, not the original whole-script replay |
| Source2022 baseline and matched screen | no | yes | Natural/other log conventions and printed weights, failed nearest-rounding tests, post-hoc truncation diagnostic, 74 finite matched nodes |

The quick profile therefore runs seven jobs and **does not recheck the
all-finite-q lower bound, QOL continuum cover, S4/S6 or 2022 source diagnostics**.
Use the full profile for those claims. Result counts
are checks, nodes or arithmetic comparisons, not independent experiments.

The MPFR output adapter changes only where a fresh result is saved. It calls
the original `verify(record)` function; it does not change formulas, precision,
acceptance tests, derivatives, input boxes or certificate trees. The other
certificate programs receive explicit new output paths through their original CLI.

The S6 adapter calls the already audited cached diagnostic's original functions,
with a 120-second equivalence budget and a 420-second complete-comparison budget.
It preserves the original 180-bisection body and every retained numerical value;
its clean receipt omits host paths. This is not a successful rerun of the original
uncached whole script, whose historical timeouts remain a limitation.

The source2022 adapter changes only input/output/provenance paths. The public
input removes one operational `access` metadata field; all other keys and nested
values, including the 50 literal raw coordinates, remain unchanged. The clean
operator note replaces internal provenance prose, not mathematics. Exact hashes,
deleted fields and original-to-public mappings are explicit in the manifest and
expected scientific projections. Baseline parity excludes only the six listed
time/name/hash/provenance fields; screen parity excludes only `utc` and
`input_pins`. Every remaining nested scientific value must equal the frozen
reference. The source PDF was not authenticated; both the failed nearest-rounding
checks and the explicitly post-hoc truncation compatibility are retained.

On Windows, the runner and adapters use extended-path prefixes for their own
resolved paths. Original programs remain unchanged. A short extraction root is
recommended for archive tools that do not support long paths. On other systems
ordinary absolute paths are used; no machine-specific path is required.

The QOL adapter calls the unchanged Decimal80 `prepare/evaluate` and independent
MPFR256 `Model.eval` functions on the frozen complete leaf cover. It replaces
only the legacy CLIs' workspace-bound paths and dependency-folder bookkeeping,
checks exact input/code hashes, rejects missing/duplicate/reversed leaf covers,
and writes to a new output. It does not rerun spreadsheet extraction or assert
binary identity of the reader's MPFR installation. Both exact input extractions
are included; the original spreadsheet hash is provenance, not a claim that the
spreadsheet was reread. Zürich is certified only within the fixed analyst model
and q interval [1,16], not under survey/source uncertainty or for every finite q.

The pinned Decimal primitive originally lives in a broader Pinar checker. The
floor16 program imports only its interval class and primitive self-tests. The
unrelated Pinar CLI and its data are not part of this module's supported jobs.
Similarly, `model_inputs.json` retains a historical unused Qin entry because
changing it would invalidate the existing exact-input guards. No Qin replay is
claimed here.

## Interpretation

At the six-decimal source, the winners at q = 34, 142, 238, 721 are C8, C7, C6,
C5 for every point in the halfwidth 1/2,000,000 raw box. Continuity then forces
at least three changes somewhere between those observations. This does not
locate crossings or establish their exact count, universal decimal sufficiency,
prevalence, or measurement calibration.

The Seikh enclosures concern the specific printed normalized four-row,
five-criterion 2023 centre and declared natural-log algebraic WA realization.
They do not silently resolve source normalization ambiguity or recover the
authors' software. They are distinct from the five-row 2022 Frank source and
from the main sharp capacity theorem.

All universal capacity/resource/Frank results remain analytic claims whose
proofs are in the manuscript and supplement; finite computation is not their
proof. No public-access, DOI, submission or journal-readiness claim is made by
this module.

## Source attribution

The four-row Seikh centre is derived numeric data from M. R. Seikh and U. Mandal,
*q-Rung Orthopair Fuzzy Archimedean Aggregation Operators: Application in the Site
Selection for Software Operating Units*, Symmetry 15(9) (2023), 1680,
https://doi.org/10.3390/sym15091680. Its formula/normalization qualifications above
are part of the audited model, not a claim to recover source-author software.

The distinct five-row source is M. R. Seikh and U. Mandal, *q-Rung orthopair fuzzy
Frank aggregation operators and its application in multiple attribute
decision-making with unknown attribute weights*, Granular Computing 7 (2022),
709–730, https://doi.org/10.1007/s41066-021-00290-2. See
`SOURCE2022_OPERATOR_SCOPE.md`. Only derived numerical transcription and original
reproduction code are included, not publisher prose, PDF or table images.
