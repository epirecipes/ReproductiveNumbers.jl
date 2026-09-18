# How the symbolic algorithms are routed

Every symbolic step of the pipeline can be done in more than one way, with different costs
and different chances of success. The package chooses an order of methods from the size of
the problem and falls through to the next method when one fails. Set

```julia
ENV["JULIA_DEBUG"] = "ReproductiveNumbers"
```

to log the route taken (or use `Logging.with_logger` with a debug-level logger).

| Step                                            | Methods, in the order tried                                                                                    | What decides the order                                                                                                              | Force a method                                                               |
|:----------------------------------------------- |:-------------------------------------------------------------------------------------------------------------- |:----------------------------------------------------------------------------------------------------------------------------------- |:---------------------------------------------------------------------------- |
| Inverting `Σ` (and the other symbolic inverses) | Laplace expansion (`inv`, compact but factorial cost) and LU (`Symbolics.sym_lu`, cubic cost)                  | Laplace first up to `LAPLACE_INVERSE_LIMIT[]` (4) compartments, LU first above; Laplace is never tried above 8                      | [`ReproductiveNumbers.symbolic_inverse`](@ref)`(A; method = :laplace / :lu)` |
| Tidying entries                                 | sign normalisation, then `simplify_fractions`                                                                  | `simplify_fractions` is skipped for entries larger than `TIDY_SIZE_LIMIT[]` (2000) expression nodes                                 | `tidy(x; fractions = false)`                                                 |
| Deciding whether an expression is zero          | numeric probes at fixed pseudo-random points; exact rational-function simplification; probabilistic acceptance | probes first because a non-zero value is cheap and conclusive; exact simplification only for rational expressions                   | `symbolic_iszero(x; numeric = false)`                                        |
| Closed-form `R₀`                                | irreducible blocks of `K` (`1 × 1`, rank one, `2 × 2`); the small-domain matrix `K_S`                          | small domain first when `K` is larger than `2 × 2` and `T` has low rank, blocks first otherwise; a `NoClosedFormError` if both fail | `basic_reproduction_number(ngm; method = :blocks / :small_domain)`           |
| Infection-free steady state                     | linear solve; polynomial solve (needs Nemo); error asking for `equilibrium`                                    | linear systems are solved exactly; polynomial systems need Nemo and a unique admissible root                                        | `equilibrium = Dict(...)`                                                    |
| Numbers                                         | everything above with parameter values substituted, then dense numeric linear algebra                          | always available; no closed form is needed                                                                                          | `basic_reproduction_number(ngm, p)`                                          |

The limits are `Ref`s that can be changed at run time, for instance
`ReproductiveNumbers.TIDY_SIZE_LIMIT[] = 10_000` to insist on exact simplification of large
entries, or `ReproductiveNumbers.LAPLACE_INVERSE_LIMIT[] = 0` to always use LU.

## Rationale

  - `inv(::Matrix{Num})` expands by cofactors. For the `12 × 12` transition matrix of a
    six-group SEIR model it would not finish; LU takes a few seconds. For a `2 × 2` matrix
    the cofactor formula gives the more readable result.
  - `Symbolics.simplify_fractions` computes polynomial GCDs of numerator and denominator. On
    the entries of small models it is instantaneous and produces the compact forms shown in
    the vignettes; on the closed form of a `2 × 2` block with square roots, or on entries
    with thousands of nodes, it can take minutes without changing the value.
  - Testing all `2 × 2` minors of a block symbolically is exact but only cheap for rational
    entries; the numeric probe rejects most non-rank-one blocks in microseconds.
  - The small-domain matrix has the same non-zero spectrum as `K` and is smaller whenever `T`
    has low rank (section 3.3 of Diekmann et al. 2010). For the vertical-transmission model of
    vignette 8, `K` is `4 × 4` with a `3 × 3` irreducible block and no closed form, while `K_S`
    is `2 × 2`.

## Retries

[`ReproductiveNumbers.first_success`](@ref) runs a ladder of `name => thunk` pairs, logs
each failure at debug level and returns the first result; the last exception is rethrown
when every method fails. New routes can be added by extending a ladder.
