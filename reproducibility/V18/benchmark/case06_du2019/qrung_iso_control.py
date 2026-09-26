"""Core q-rung representation-control utilities.

Canonical powered coordinates:
    Phi_q(mu,nu) = (mu**q, nu**q)

Representation-preserving transport:
    T_{r->q}(mu,nu) = (mu**(r/q), nu**(r/q))
"""
from math import prod

def admissible(x, q, tol=1e-15):
    mu, nu = x
    return mu**q + nu**q <= 1.0 + tol

def minimum_integer_rung(values, qmax=100):
    for q in range(1, qmax + 1):
        if all(admissible(x, q) for x in values):
            return q
    raise ValueError("No admissible integer rung found in search range.")

def transport_pair(x, r, q):
    mu, nu = x
    return mu**(r/q), nu**(r/q)

def transport_row(row, r, q):
    return [transport_pair(x, r, q) for x in row]

def qrof_wpm(row, weights, q, power=2.0):
    """Du (2019) q-ROF weighted power mean, for positive power exponent."""
    mem_prod = prod((1.0 - mu**(power*q))**w
                   for (mu, _), w in zip(row, weights))
    non_prod = prod((1.0 - (1.0 - nu**q)**power)**w
                   for (_, nu), w in zip(row, weights))
    mu_out = (1.0 - mem_prod)**(1.0/(power*q))
    inner = 1.0 - (1.0 - non_prod)**(1.0/power)
    nu_out = inner**(1.0/q)
    return mu_out, nu_out

def q_score(x, q):
    mu, nu = x
    return mu**q - nu**q

def euclidean_closeness(x, p=2.0):
    """Raw-coordinate closeness endpoint used as a residual-dependence probe."""
    mu, nu = x
    numerator = (1.0-mu)**p + nu**p
    denominator = mu**p + (1.0-nu)**p
    return 1.0 / (1.0 + (numerator/denominator)**(1.0/p))

def rank_desc(scores):
    return tuple(sorted(scores, key=lambda k: (-scores[k], k)))
