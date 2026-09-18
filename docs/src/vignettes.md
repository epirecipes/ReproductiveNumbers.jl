# Vignettes

The vignettes are Quarto documents (Julia engine) in the `vignettes/` directory of the
repository, each in its own folder. They are rendered by the documentation workflow and
bundled into this site; the links below open the rendered pages. To render them locally:

```sh
cd vignettes
julia --project=. -e 'using Pkg; Pkg.develop(path=".."); Pkg.instantiate()'
quarto render
```

| # | Vignette                                                                                           | What it shows                                                                                                                     |
|:- |:-------------------------------------------------------------------------------------------------- |:--------------------------------------------------------------------------------------------------------------------------------- |
| 1 | [The SIR model](vignettes/01-sir/sir.html)                                                         | `T`, `Σ`, `K` and `R₀ = β/γ`; numeric evaluation; the threshold by simulation; a reproduction number at any state                 |
| 2 | [SEIR with demography](vignettes/02-seir/seir.html)                                                | Large versus restricted domain, states-at-infection, automatic infection-free steady state, characteristic polynomial             |
| 3 | [Two latent stages](vignettes/03-two-latent-stages/two-latent-stages.html)                         | Example 2.1 of Diekmann et al.: a rank-one `K` and the small-domain matrix                                                        |
| 4 | [Vector-borne transmission](vignettes/04-vector-borne/vector-borne.html)                           | Ross–Macdonald: `R₀ = √(…)` and what changes when only host-to-host generations are counted                                       |
| 5 | [Reaction networks with Catalyst](vignettes/05-catalyst/catalyst.html)                             | Transmission reactions identified from stoichiometry; vertical transmission; overriding by predicate, flags or indices            |
| 6 | [Heterogeneous mixing](vignettes/06-heterogeneous-mixing/heterogeneous-mixing.html)                | Two host types (`2 × 2` closed form), separable mixing, multi-strain block structure, an age-structured model that needs numerics |
| 7 | [Type reproduction numbers](vignettes/07-type-reproduction-numbers/type-reproduction-numbers.html) | Targeted control of a core group; type reproduction numbers in the vector-borne model                                             |
