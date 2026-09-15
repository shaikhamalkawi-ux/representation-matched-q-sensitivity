# Public release QA — 2026-09-16

Status: PASS for the public quick workflow.

Verified in a clean working copy:

- Seven-case executable benchmark: PASS.
- Benchmark release verifier: PASS.
- Qin stagewise verification: PASS (54/54 checks in the independent verification layer).
- Structural controls: PASS (21/21).
- Negative checker fixtures: PASS (4/4 rejected as intended).
- Directed upper witness: PASS.
- Directed A1-A3 margin interval: [-3.3874873674535827e-16, -3.3874873674535820e-16].
- Directed L2 witness-distance interval: [0.003514796840391997, 0.0035147968403920335].
- The retained radius-cover replay metadata are included in GitHub. The larger compressed cover trace is distributed with the Zenodo archive; a full rebuild remains available with `python run_all.py --full`.

Scientific boundary remains unchanged: exact global winner-reversal radius is unknown / not claimed.
