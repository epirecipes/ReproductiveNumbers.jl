# Mathematical background

This page summarises the construction the package implements. Notation follows
Diekmann, Heesterbeek and Roberts (2010) (DHR below); the equivalent notation of van den
Driessche and Watmough (2002) (vdDW) is given where it differs.

## The infected subsystem and its linearisation

Write the model as ``\dot{x} = f(x, y)``, ``\dot{y} = g(x, y)`` where ``x`` collects the
*infected* compartments (those an individual can occupy while carrying the infection, such
as exposed, infectious, or treated) and ``y`` the uninfected ones. Let ``(0, y^*)`` be the
infection-free steady state. Near it the infected subsystem is linear,

```math
\dot{x} = J\,x, \qquad J = \left.\frac{\partial f}{\partial x}\right|_{(0, y^*)},
```

and ``J`` is split as ``J = T + \Sigma`` where

  - ``T_{ij}`` is the rate at which an individual in infected state ``j`` produces *new*
    infections in state ``i`` (transmissions; vdDW's ``F``), and
  - ``\Sigma_{ij}`` collects every other change of state: progression, recovery, death
    (transitions; vdDW's ``-V``).

``T \ge 0`` entrywise, the off-diagonal entries of ``\Sigma`` are non-negative and its
diagonal entries negative, and ``-\Sigma^{-1} \ge 0``: its entry ``(i, j)`` is the expected
time an individual now in state ``j`` will spend in state ``i`` over its remaining
infected life. [`validate_decomposition`](@ref) checks these signs numerically.

## Next-generation matrices

The **next-generation matrix with large domain** is

```math
K_L = -T\,\Sigma^{-1},
```

whose entry ``(i, j)`` is the expected number of new infections in state ``i`` produced by
an individual currently in state ``j``. Only some states can be entered by infection: the
**states-at-infection** are the non-zero rows of ``T``. Let ``E`` be the matrix whose
columns are the corresponding unit vectors. The **next-generation matrix** is the
restriction

```math
K = E^{\mathsf T} K_L E = -E^{\mathsf T} T \Sigma^{-1} E,
```

indexed by states-at-infection, whose entries have the classical interpretation: the
expected number of new infections of type ``i`` produced over its whole infected life by
one individual who entered the infected population in state ``j``. ``K`` and ``K_L`` have
the same non-zero eigenvalues (because ``T = E E^{\mathsf T} T``, so ``K_L = B A`` and
``K = A B`` with ``A = E^{\mathsf T} K_L`` and ``B = E``), hence

```math
R_0 = \rho(K) = \rho(K_L).
```

When the states-at-infection are always entered in fixed proportions, ``T = C R`` has rank
one and DHR reduce further to the **small-domain** matrix ``K_S = -R \Sigma^{-1} C``,
which is a scalar; see [`small_domain_matrix`](@ref).

## The threshold property

``R_0 > 1`` if and only if the infection-free steady state is unstable, that is the
spectral bound ``s(T + \Sigma)`` is positive, and ``R_0 = 1`` exactly when
``\det(T + \Sigma) = 0``, since ``T + \Sigma = (K_L - I)(-\Sigma)``. The regression tests
check this at random parameter values and against simulated epidemics.

## What counts as a transmission is a choice

Different splittings of ``J`` into ``T`` and ``\Sigma`` (all satisfying the sign
conditions) give different next-generation matrices and different values of ``R_0``, but
the same threshold. DHR's vector–host example is the standard illustration: counting a
mosquito infected by a human as a new infection gives ``R_0 = \sqrt{K_{12} K_{21}}``,
while counting only human infections gives ``K_{12} K_{21}``. The package offers a
default and lets you override it, see
[Choosing what counts as a transmission](@ref).

## How the package computes ``R_0`` symbolically

 1. [`irreducible_blocks`](@ref) permutes ``K`` to block triangular form using the strongly
    connected components of its non-zero pattern; ``\rho(K)`` is the maximum of the
    spectral radii of the diagonal blocks.

 2. Each irreducible block is handled in closed form if it is ``1 \times 1``, has rank one
    (``\rho = \operatorname{trace}``), or is ``2 \times 2``:
    
    ```math
    \rho\begin{pmatrix} a & b \\ c & d \end{pmatrix}
      = \frac{a + d + \sqrt{(a - d)^2 + 4bc}}{2},
    ```
    
    which reduces to ``\sqrt{bc}`` when the diagonal vanishes.
 3. Otherwise a [`NoClosedFormError`](@ref) is thrown and `basic_reproduction_number(ngm, p)`
    evaluates ``\rho(K)`` numerically.

Rank is tested through the ``2 \times 2`` minors: symbolically when they are rational
functions of the parameters, and by evaluation at random points as a fallback (see
[`ReproductiveNumbers.symbolic_iszero`](@ref)). Each of these steps is proved correct in
Lean, see [Formal proofs](@ref).

## Type reproduction numbers

For a set ``S`` of states-at-infection with projection ``P``, the type reproduction
number of Roberts and Heesterbeek (2003) is

```math
T_S = \rho\bigl(P K (I - (I - P) K)^{-1}\bigr),
```

the expected number of infections of types in ``S`` produced by an individual of a type in
``S``, counting chains through the other types. Reducing transmission from ``S`` by a
factor ``1 - 1/T_S`` eliminates the infection, and ``T_S > 1 \iff R_0 > 1`` provided
``\rho((I - P)K) < 1``. See [`type_reproduction_number`](@ref).

## References

  - Diekmann O, Heesterbeek JAP, Roberts MG (2010). The construction of next-generation
    matrices for compartmental epidemic models. *J R Soc Interface* 7:873–885.
  - van den Driessche P, Watmough J (2002). Reproduction numbers and sub-threshold endemic
    equilibria for compartmental models of disease transmission. *Math Biosci* 180:29–48.
  - Roberts MG, Heesterbeek JAP (2003). A new method for estimating the effort required to
    control an infectious disease. *Proc R Soc Lond B* 270:1359–1364.
  - Diekmann O, Heesterbeek JAP, Metz JAJ (1990). On the definition and the computation of
    the basic reproduction ratio ``R_0`` in models for infectious diseases in heterogeneous
    populations. *J Math Biol* 28:365–382.
