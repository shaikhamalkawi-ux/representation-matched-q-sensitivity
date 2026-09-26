# QROF V18 reproducibility artifact

Companion code and numerical evidence for **Representation-Matched q-Sensitivity:
Rank Changes and Winner Robustness in Orthopair Fuzzy Models**.

This is a curated scientific reproducibility release, not the private manuscript
delivery archive. It contains author-written code, small numerical model inputs
transcribed from cited published examples, derived outputs, exact certificate
records and verification programs. Publisher PDFs, runtime binaries, private
history, old manuscript drafts and nested private archives are not included.

## Fresh verification

Provide Python 3.10+ with NumPy, SciPy and mpmath, and a preinstalled C++17 compiler
with MPFR/GMP development headers and libraries. The programs perform no network
access or dependency installation. Do not run Python with `-O`.

From the extracted release directory:

```sh
python -B verify_manifest.py
python -B reproduce.py --output ../qrof-v18-fresh
```

For a compiler prefix not on PATH, append `--toolchain /path/to/compiler-prefix`.
The prefix must contain `bin`; on Windows use a compatible MinGW-w64 toolchain
with MPFR/GMP DLLs in that directory. Optionally use `--python-deps /path/to/deps`
for an existing directory containing mpmath. No machine-specific path is built
into this release. The output directory must not already exist and must be
outside the release. All calculations run in a new working copy, preserving the
released source and certificate files.

The slower historical Qin q=3 radius-cover verification is explicit:

```sh
python -B reproduce.py --output ../qrof-v18-full-check --radius-replay
```

This validates the retained 18,515-node cover; it does not rerun its optimizer or
cover search. To review only the new Seikh module, with mpmath and MPFR/GMP:

```sh
python -B seikh/reproduce_v17.py --output ../qrof-seikh-fresh
```

The `v17` names inside that module intentionally identify the unchanged scientific
records incorporated in V18. V18 changes presentation and public availability,
not those model formulas or certificate values.

## Evidence map

| Directory | Included evidence | Verification scope |
|---|---|---|
| `seikh/` | Four integer certificates, three rational failure witnesses, model input pins and independent MPFR implementation | Regenerates all substantive certificate fields; 40 positive leaves, three witnesses, corruption tests |
| `benchmark/` | Seven published-example implementations, attribution register and derived comparison ledger | Executes all seven deterministic examples and theory/ledger checks |
| `analysis/` | Qin source-form and He structural controls, stagewise checks, q=3 lower-radius cover and directed upper witness | Quick checks by default; 18,515-node retained cover via `--radius-replay` |
| `qin/unbounded/` | H/A all-finite-real-q>=3 winner covers on canonical cube radius 1/4096; compactified kernels and independent source checks | Fresh compilation and 256-bit replay of retained covers, including endpoint and malformed-record tests |
| `qin/fullorder/` | H/A complete-order covers for real q in [3,10], canonical cube radius 1/65536 | Standalone topology/domain checker and freshly compiled kernel |
| `zhang/` | Inverse-score and ordered-weight source-control implementation and derived tables | Regenerates 15-score control and source-pipeline comparisons |

See `SCOPE_AND_PROVENANCE.md` for model conventions, what is deliberately missing,
the trust boundary and interpretation of historical records. `PUBLIC_QA.json`
records the verification actually performed for this release; included code is
not by itself evidence that every historical search was rerun.

## Main Seikh result retained in this release

For the same admissible 40-coordinate normalized baseline-raw L-infinity domain,
center, transport and recomputed entropy-weight model:

`rho_src(4) <= 0.02321968 < 0.0235 <= rho_src(16)`.

The uniform finite-real-q>=4 radius is bracketed by
`0.019 <= rho_src <= 0.02321968`.

These are strict ordering and certified bounds, not exact optimal radii or a
claim of monotonicity for every q. The q=8 bracket overlaps the other brackets.
Finite tests do not establish formula identities; the proof arguments and the
stated perturbation domains are in the manuscript and supplement.

## Rights

`LICENSE` reproduces the repository's existing CC0 notice byte-for-byte and applies
only to original material released by the authors. Third-party source articles,
data rights, libraries and dependencies remain subject to their own terms. No
MPFR, GMP, NumPy, SciPy, mpmath or compiler source/binary is redistributed here.
Published-model numerical transcriptions are attributed in code and source
registers; the CC0 notice does not purport to relicense those sources.
