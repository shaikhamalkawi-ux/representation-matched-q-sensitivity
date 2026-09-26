# Case 3 provenance — Liu & Wang 2018

Primary source: Peide Liu and Peng Wang, “Some q-Rung Orthopair Fuzzy Aggregation Operators and their Applications to Multiple-Attribute Decision Making”, *International Journal of Intelligent Systems* 33(2), 259–280. DOI 10.1002/int.21927.

The executable endpoint uses Example 4, Table III, the published attribute weights, the standard q-ROFWA operator, Table IV at q=3 and Table V sensitivity results.

A 2026 correction (DOI 10.1155/int/9864340) corrects Definition 1's indeterminacy degree. It does not replace the q-ROFWA operator, Example 4 matrix, or Table V.

## Bound q-ROFWG audit

The source also reports q-ROFWG results in Table VI. Applying the printed q-ROFWG operator to the published Table III data reproduces X2–X5, but not X1. Therefore q-ROFWG is not used as a clean benchmark endpoint; the discrepancy is preserved in `results/qrofwg_mismatch_audit.json`.
