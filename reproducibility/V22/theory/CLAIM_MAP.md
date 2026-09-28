# V22 claim-to-evidence map

Target: V22 reviewed R4. The following source identifiers locate claims in the
manuscript supplied separately to authors/reviewers. This package does not
publish the private manuscript or substitute code for its proofs.

Core source SHA-256: `055a2e53dd116e990cbd35cceadc7a8397bd11c52c9a4e09ce543b3213b778e5`.
Supplement source SHA-256: `763cf2f8cf144d5905a64686c88253a3c23e5f8547b53567e6be2e58349da93a`.
These identify the reviewed mathematical text; they are not proof certificates.

## Analytic results: proof lives in the manuscript

| Claim/location | Exact scope | What this public module does not claim |
|---|---|---|
| Theorem 1 / `thm:joint-capacity`; Supplement S1–S3 | Algebraic WA, two benefit criteria, positive strictly r-admissible fixed table, globally common raw nonmembership, recomputed same-table natural-log weights; r >= 1, m >= 4, 0 < gamma <= B/2, B = 77 ln(2)/122880. Sharp order Theta(min(m², sqrt(m/gamma), 1/gamma)) for chronological changes between all-rival-margin-qualified observations | No code certifies the universal quantifiers; no exact leading constant, arbitrary-criterion sharp law, empirical prevalence or priority clearance |
| `cor:integer-rational`; S3/S10 | Rational distinct-row integer-observation existence; quantitative rung/precision control only in the stated balanced central-scale construction | General existence is not an efficient algorithm or universal decimal-resolution guarantee |
| `cor:unequal-margin`; S7/S10 | Two criteria with arbitrary nonmembership; matching central/larger-margin order only for B/m³ <= gamma <= B/2 | No full low-margin arbitrary-nonmembership law or externally prescribed profile attainment |
| `prop:resource-qualified`, `cor:budget-capacity`, `cor:analytic-resources`; S8/S10 | Algebraic two-product scores; any fixed d for the cap/window upper bounds. Qualified count limited by all-grade raw-log cap and logarithmic rung span; matching two-criterion budget law only on its stated N/m/gamma range | These are analytic estimates, not numerically tested universal bounds. Membership-only necessity needs common nonmembership; ordinary rational encoding incurs exponential auxiliary bit cost |
| `cor:frank-capacity`; S11 | Frank WA parameter exactly 2, same-table natural-log weights, two criteria, globally common nonmembership; m >= 4, r >= 1, 0 < gamma <= 77/706560, real finite observations | No numerical certificate here proves this analytic extension. No all-parameter, WG, unequal-nonmembership, integer/rational or resource refinement is asserted for Frank |

The argument credits classical exponential-zero, concavity/tangent and
generator methods. Numerical replay does not establish that the model-specific
same-table construction is novel relative to all prior literature. That remains
a substantive specialist assessment, not a successful-test flag.

## Executable numerical evidence supplied here

| V22 location | Exact claim | Program/data |
|---|---|---|
| Main Table 1 / S9 (`tab:compact-path`, `tab:s-floor16-nodes`) | Four strict winners throughout the 64-coordinate half-unit raw box, floor 0.1249995, baseline r=4 | `check_floor16_mpfr.py`, `check_floor16_decimal.py`, pinned `FLOOR16_SOURCE_AND_INTEGER_NODES_R1.json`, MPFR reference record; runner jobs `floor16_*` |
| S5 `tab:s-radii`, q=4/8/16 lower bounds | rho_src >= 0.0232195 / 0.0250515 / 0.0258973 respectively; exact endpoints; complete source boxes | Three `Q*_H_*.json` records; unchanged integer solver/replayer and independent MPFR verifier |
| S5 uniform lower | inf over every finite real q >= 4 of rho_src(q) >= 0.0232195; parameter t=1/q includes zero only as continuous endpoint | `UNIFORM_JOINT_H_0_0232195.json`; full-profile integer/MPFR jobs, 811 nodes and 407 positive leaves |
| S5 witness uppers | rho_src(4) <= 0.02321968; rho_src(8) <= 0.02505164; rho_src(16) <= 0.02589746; uniform upper follows q=4 | `FIXED_RUNG_WITNESSES.json` and `verify_fixed_rung_witnesses.py`: exact feasible rational sources, unique Y1; point witnesses are not lower proofs |
| Main observed-data section / S12 QOL continuum control | Zürich remains strict WA model winner for every q in [1,16], both two/ten-criterion sets, 83 cities | Unchanged Decimal80/MPFR256 evaluators, their two exact input tables and frozen 32/61-leaf cover; full-profile `qol_*_cover` jobs, no new spreadsheet extraction |
| S4 fragile seven-row construction | Seven fixed strict winners throughout the all-28-coordinate halfwidth 10^-11 box, hence at least six changes; auxiliary rows globally excluded | `COMPACT_SOURCE_AND_NODES.json`, unchanged directed MPFR256/Decimal100 and exact rational domain checker; no optimum-radius or calibrated-error claim |
| S4 precision failure qualification | Simultaneously rounding every raw grade at fixed nodes fails at 6–9 places; all 6–15 tests and 16-place baseline retained | `check_decimal_rounding.py`; directed Decimal100 specific-protocol results, not a universal precision lower bound |
| S6 diagnostic crosswalk | 18 particular ordinary Decimal90 instances match their frozen diagnostic values | Original script and JSON plus unchanged cached checker; 252 literal probes and same 180-bisection body. Not a continuum certificate or original uncached whole-script completion |
| Separate 2022 Frank source control | Five-row baseline with cost swaps, declared natural-log weights, nearest-rounding mismatches, post-hoc truncation agreement; no observed winner change at 74 finite matched nodes | Original Decimal100 scripts, clean 50-coordinate input/provenance, exact scientific-result comparison. HTML source only; no all-q or source-author convention claim |

The six-decimal minimum all-rival gap lower bounds printed in Table 1 are
0.000223681874187378, 0.000014183414599287, 0.000014821884633083 and
0.000014674290057407 in increasing q order. They are rounded downward from the
independent outward Decimal endpoints. Lower covers and upper witnesses jointly
separate the three S5 radii, but do not give exact optima or monotonicity at all q.

## Explicit limits of this module's coverage

- It is the **theory/certificate module**, not the complete combined V22
  reproducibility release. The separate empirical module and top-level inventory
  must state their own coverage and results.
- The full profile replays the listed current numerical claims, not every
  historical discovery scan, optimizer run, failed construction search, source
  extraction or manuscript proof. Expected records never substitute for fresh
  numerical checks in a reported successful job.
- No proof-assistant certification, source-author code recovery, calibrated
  uncertainty, reviewer access or completed human scientific approval follows.
