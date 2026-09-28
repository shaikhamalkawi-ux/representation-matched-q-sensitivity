# Public provenance and model convention: Seikh and Mandal (2022)

M. R. Seikh and U. Mandal, *q-Rung orthopair fuzzy Frank aggregation operators
and its application in multiple attribute decision-making with unknown
attribute weights*, Granular Computing 7 (2022), 709–730.
DOI: https://doi.org/10.1007/s41066-021-00290-2.

The exact 5-by-5 numeric input was transcribed from the article's authorized
full-text HTML, not visually authenticated against a downloaded source PDF.
This public note restates the audited mathematical conventions; it replaces
only internal provenance prose, not a formula or numeric input. The original
provenance-note SHA-256 is recorded in the curated baseline projection's
original-record provenance. No subscription article/PDF is redistributed.

For normalized pairs (mu_j,nu_j), normalized weights w_j, and tau=2, define

    F_tau(x,w) = log_tau(1 + product_j (tau^(x_j)-1)^(w_j)).

The source Eq. 7 WA powered components are

    mu_out^q = 1 - F_tau(1-mu_j^q,w)
    nu_out^q = F_tau(nu_j^q,w).

The source Eq. 8 WG powered components are

    mu_out^q = F_tau(mu_j^q,w)
    nu_out^q = 1 - F_tau(1-nu_j^q,w).

The source score is (1+mu_out^q-nu_out^q)/2. Cost columns 1, 4, 5 swap the
two grades. Raw-grade entropy uses e_j=1+mean_i(mu_ij log(mu_ij)+nu_ij
log(nu_ij)), normalized over columns. The inspected source does not explicitly
specify the logarithm base: natural log is the declared audited model and is
numerically compatible with the printed Eq. 10 weights, not an attributed
source-author declaration. Other log bases and the printed weights remain
separate baseline diagnostics.

The matched screen fixes canonical powers at r=4 and recomputes these raw-grade
weights along (mu,nu)->(mu^(4/q),nu^(4/q)). It does not reproduce the article's
different fixed-raw/tau=4 sensitivity procedure. All 74 finite nodes and the
equal-weight limit evaluation are ordinary Decimal100, not directed or all-q
certificates. Nearest-four-decimal score mismatches and the post-hoc truncation
compatibility are both retained; truncation is not asserted to be the authors'
rounding convention.
