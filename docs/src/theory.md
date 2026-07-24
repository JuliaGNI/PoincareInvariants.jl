# Theory

Poincaré integral invariants are built from two ingredients: *differential forms*, which are the objects that can be integrated over families of trajectories, and *integral invariants*, which are the differential forms whose integrals are preserved by a dynamical flow. This page condenses the minimum theory required to state the first and second invariants computed by this package. A more detailed treatment of differential forms can be found in [Abraham & Marsden](https://doi.org/10.1007/978-1-4612-1029-0).

## Differential forms

Let the phase space be ``P = \mathbb{R}^d`` with cartesian coordinates ``z^i``, ``i = 1, \dots, d``. The coordinate differentials ``dz^i`` are treated as ``d`` independent, anti-commuting quantities, so that ``dz^i \wedge dz^j = - \, dz^j \wedge dz^i``, where ``\wedge`` denotes the (wedge) product. A homogeneous polynomial of degree ``k`` in these differentials is an *algebraic ``k``-form*, and a *differential ``k``-form* (or just ``k``-form) is an assignment of an algebraic ``k``-form to each point ``\bm{z} \in P``,
```math
\alpha (\bm{z}) = \sum_{i_1 < \dots < i_k} \alpha_{i_1 \dots i_k} (\bm{z}) \, dz^{i_1} \wedge \dots \wedge dz^{i_k} ,
```
with smooth real-valued coefficients ``\alpha_{i_1 \dots i_k} (\bm{z})``. By anti-commutativity the only ``k``-form with ``k > d`` is ``0``.

A ``k``-form can be integrated over a parametrized ``k``-dimensional submanifold ``S \subset P``. Suppose ``S`` is parametrized by smooth functions ``\bm{s} (\bm{u}) = (s^1 (\bm{u}), \dots, s^d (\bm{u}))`` with ``\bm{u} = (u^1, \dots, u^k)`` ranging over a parameter set ``U \subset \mathbb{R}^k``. Applying the chain rule to ``z^i = s^i (\bm{u})`` pulls the coordinate differentials back to the parameter differentials,
```math
dz^i = \sum_{j=1}^{k} \frac{\partial s^i}{\partial u^j} \, du^j ,
```
and each order-``k`` monomial collapses to a Jacobian determinant,
```math
dz^{i_1} \wedge \dots \wedge dz^{i_k} = \frac{\partial (s^{i_1}, \dots, s^{i_k})}{\partial (u^1, \dots, u^k)} \, du^1 \wedge \dots \wedge du^k .
```
The integral of ``\alpha`` over ``S`` is then the ordinary integral over ``U``,
```math
\int_S \alpha = \int_U \left( \sum_{i_1 < \dots < i_k} \alpha_{i_1 \dots i_k} (\bm{s} (\bm{u})) \, \frac{\partial (s^{i_1}, \dots, s^{i_k})}{\partial (u^1, \dots, u^k)} \right) du^1 \dots du^k .
```
Although this definition uses a particular parametrization, the change-of-variables rule for Jacobian determinants guarantees that the value of the integral is independent of the parametrization chosen for ``S``.

## Integral invariants

Consider a general non-autonomous ordinary differential equation on ``P``,
```math
\dot{\bm{z}} = X_t (\bm{z}) ,
```
where ``X_t`` is a time-dependent vector field. An **absolute integral invariant** is a time-dependent ``k``-form ``\alpha_t`` such that the integral
```math
I_t = \int_{S_t} \alpha_t
```
is independent of time whenever ``S_t`` is a compact ``k``-dimensional submanifold advected along solutions of the ODE. A **relative integral invariant** is one for which ``I_t`` is constant only when ``S_t`` is in addition *closed* (compact and without boundary).

The classic example is *Kelvin's circulation theorem* for the barotropic Euler equations: the ``1``-form ``\bm{u}_t (\bm{x}) \cdot d\bm{x}`` is a relative integral invariant, while its exterior derivative, the vorticity ``2``-form built from ``\bm{\omega}_t = \nabla \times \bm{u}_t``, is an absolute integral invariant. This illustrates a general fact used repeatedly below: the exterior derivative of a relative integral invariant is always an absolute integral invariant.

### Canonical Hamiltonian systems

Let ``P = \mathbb{R}^n \times \mathbb{R}^n \ni (\bm{q}, \bm{p})`` with coordinates ``q^i, p_i``. A canonical Hamiltonian system arises as the Euler–Lagrange equations of the action
```math
S[(\bm{q}, \bm{p})] = \int_{t_1}^{t_2} \left( \sum_{i=1}^{n} p_i \, \dot{q}^i - H_t (\bm{q}, \bm{p}) \right) dt ,
```
that is, Hamilton's equations
```math
\dot{q}^i = \frac{\partial H_t}{\partial p_i} , \qquad \dot{p}_i = - \frac{\partial H_t}{\partial q^i} .
```
The Hamiltonian ``H_t`` is a conserved scalar only when it is time-independent, but canonical Hamiltonian systems always admit integral invariants, regardless of any time-dependence of ``H_t``. The fundamental one is the order-one relative invariant, the Lagrangian one-form
```math
\vartheta = \sum_{i=1}^{n} p_i \, dq^i ,
```
whose exterior derivative is the order-two absolute invariant, the symplectic two-form
```math
\omega = \sum_{i=1}^{n} dq^i \wedge dp_i .
```
Higher-order invariants are obtained by taking wedge powers: at each odd order a relative invariant is given by ``\vartheta \wedge \omega``, ``\vartheta \wedge \omega \wedge \omega``, ``\dots``, and at each even order an absolute invariant by ``\omega``, ``\omega \wedge \omega``, ``\dots``. In fact, canonical Hamiltonian systems are exactly those ODEs that admit ``\vartheta`` as a relative integral invariant.

### Noncanonical Hamiltonian systems

A broader class arises as the Euler–Lagrange equations of the action
```math
S[\bm{z}] = \int_{t_1}^{t_2} \left( \sum_{i=1}^{d} \vartheta_{ti} (\bm{z}) \, \dot{z}^i - H_t (\bm{z}) \right) dt ,
```
for a time-dependent one-form ``\vartheta_t = \sum_i \vartheta_{ti} \, dz^i`` and function ``H_t``. As in the canonical case, the order-one relative integral invariant is ``\vartheta_t`` and the order-two absolute integral invariant is its exterior derivative
```math
\omega_t = \sum_{i < j} \left( \frac{\partial \vartheta_{tj}}{\partial z^i} - \frac{\partial \vartheta_{ti}}{\partial z^j} \right) dz^i \wedge dz^j ,
```
with higher-order invariants again constructed from wedge products of these.

The first and second Poincaré invariants ``I_1`` and ``I_2`` computed by this package are the integrals of ``\vartheta`` and ``\omega`` over an advected loop and surface, respectively; see the [Home](index.md) page and the guides for their numerical evaluation.
