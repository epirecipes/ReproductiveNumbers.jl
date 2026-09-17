# Vignettes

Each subdirectory holds one Quarto document (Julia engine) that can be rendered on its own:

```sh
cd vignettes
julia --project=. -e 'using Pkg; Pkg.develop(path=".."); Pkg.instantiate()'
quarto render 01-sir/sir.qmd          # one vignette
quarto render                         # all of them
```

The vignettes go from simple to complex:

| # | Vignette                       | What it shows                                                                                                                 |
|:- |:------------------------------ |:----------------------------------------------------------------------------------------------------------------------------- |
| 1 | `01-sir`                       | The SIR model: `T`, `Σ`, `K` and `R₀ = β/γ`; numeric evaluation                                                               |
| 2 | `02-seir`                      | SEIR with demography: `K_L` versus `K`, states-at-infection, automatic infection-free steady state                            |
| 3 | `03-two-latent-stages`         | Example 2.1 of Diekmann et al.: a rank-one `K` and the small-domain matrix                                                    |
| 4 | `04-vector-borne`              | Ross–Macdonald malaria: `R₀ = √(…)` and what changes when only host-to-host generations are counted                           |
| 5 | `05-catalyst`                  | Reaction networks with Catalyst: transmission reactions identified from stoichiometry, vertical transmission                  |
| 6 | `06-heterogeneous-mixing`      | Two host types (2×2 closed form), separable mixing, multi-strain block structure, an age-structured model that needs numerics |
| 7 | `07-type-reproduction-numbers` | Type reproduction numbers and targeted control                                                                                |
