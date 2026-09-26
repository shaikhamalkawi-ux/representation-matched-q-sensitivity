# Isomorphism-controlled q-sensitivity — formal core

Let

\[
D_q=\{(\mu,\nu)\in[0,1]^2:\mu^q+\nu^q\le 1\},
\qquad
\Delta=\{(u,v)\in[0,1]^2:u+v\le1\}.
\]

Define the powered-coordinate chart

\[
\Phi_q:D_q\to\Delta,\qquad
\Phi_q(\mu,\nu)=(\mu^q,\nu^q),
\]

with inverse

\[
\Phi_q^{-1}(u,v)=(u^{1/q},v^{1/q}).
\]

For r,q >= 1, define the representation-preserving transport

\[
T_{r\to q}=\Phi_q^{-1}\circ\Phi_r,
\qquad
T_{r\to q}(\mu,\nu)=(\mu^{r/q},\nu^{r/q}).
\]

Then

\[
T_{q\to s}\circ T_{r\to q}=T_{r\to s},
\qquad
T_{r\to q}^{-1}=T_{q\to r}.
\]

The construction is a specialization/use of known lattice-isomorphism theory.
The novelty claim of this project is **not** the existence of the isomorphism.
The prospective contribution is its use as a counterfactual control for q-sensitivity.

## q-natural stage theorem

Suppose an n-ary stage has the factorization

\[
F_q=\Phi_q^{-1}\circ G\circ\Phi_q^{\otimes n},
\]

where G is independent of q. Then

\[
F_q\circ T_{r\to q}^{\otimes n}
=
T_{r\to q}\circ F_r.
\]

Thus F is transport-equivariant.

If a scalar functional has the form

\[
s_q=h\circ\Phi_q,
\]

then

\[
s_q(T_{r\to q}x)=s_r(x).
\]

Therefore, a composition of q-natural aggregation stages followed by
q-natural scalar/ranking stages is invariant under the controlled change of q.

## TOPSIS corollary

For the Alkan-Kahraman pipeline:

- q-ROFWG aggregation is q-natural.
- q-ROF multiplication is q-natural in powered coordinates.
- cost normalization C(mu,nu)=(nu,mu) commutes with transport.
- the score S_q=mu^q-nu^q and accuracy H_q=mu^q+nu^q are invariant.
- Euclidean distances in Eqs. (29)-(30) and (35)-(36) use powered
  coordinate differences, so they are invariant.
- entropy in Eqs. (32)-(34) is a function only of powered coordinates;
  hence the entropy-derived crisp criterion weights are invariant.

Consequently both complete controlled TOPSIS pipelines are expected to
be invariant, apart from floating-point roundoff.

## Non-additive sensitivity decomposition

Do not write “raw sensitivity = representation effect + residual method
effect” as a numerical identity. For any metric d on rankings define

\[
D_{\rm raw}=d(R_{\rm raw},R_0),\quad
D_{\rm repr}=d(R_{\rm raw},R_{\rm ctrl}),\quad
D_{\rm resid}=d(R_{\rm ctrl},R_0).
\]

Then

\[
D_{\rm raw}\le D_{\rm repr}+D_{\rm resid}
\]

by the triangle inequality. For a fully q-natural pipeline,

\[
D_{\rm resid}=0.
\]
