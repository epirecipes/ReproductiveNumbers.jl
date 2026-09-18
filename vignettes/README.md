# Vignettes

Each subdirectory holds one Quarto document (Julia engine) that can be rendered on its own:

```sh
cd vignettes
julia --project=. -e 'using Pkg; Pkg.develop(path=".."); Pkg.instantiate()'
quarto render 01-sir/sir.qmd          # one vignette
quarto render                         # all of them
```

The vignettes go from simple to complex:

| #  | Vignette                           | What it shows                                                                                                                 |
|:-- |:---------------------------------- |:----------------------------------------------------------------------------------------------------------------------------- |
| 1  | `01-sir`                           | The SIR model: `T`, `Σ`, `K` and `R₀ = β/γ`; numeric evaluation                                                               |
| 2  | `02-seir`                          | SEIR with demography: `K_L` versus `K`, states-at-infection, automatic infection-free steady state                            |
| 3  | `03-two-latent-stages`             | Section 2.1 of Diekmann et al.: a rank-one `K` and the small-domain matrix                                                    |
| 4  | `04-vector-borne`                  | Ross–Macdonald malaria: `R₀ = √(…)` and what changes when only host-to-host generations are counted                           |
| 5  | `05-catalyst`                      | Reaction networks with Catalyst: transmission reactions identified from stoichiometry, vertical transmission                  |
| 6  | `06-heterogeneous-mixing`          | Two host types (2×2 closed form), separable mixing, multi-strain block structure, an age-structured model that needs numerics |
| 7  | `07-type-reproduction-numbers`     | Type reproduction numbers and targeted control                                                                                |
| 8  | `08-vertical-transmission-sti`     | Section 4.2 of Diekmann et al.: `K` has no closed form but the small-domain matrix `K_S` does                                 |
| 9  | `09-bovine-viral-diarrhoea`        | Section 4.3 of Diekmann et al.: persistently infected animals as a second state-at-infection                                  |
| 10 | `10-staged-progression`            | Staged progression: `K` stays `1×1` while `K_L` grows; elasticities by stage                                                  |
| 11 | `11-multi-strain`                  | Cross-immunity, invasion reproduction numbers, mutation and co-infection                                                      |
| 12 | `12-waning-immunity`               | SEIRS with reinfection: `R₀` checked against steady states from SteadyStateDiffEq and NonlinearSolve                          |
| 13 | `13-within-host`                   | Viral dynamics within a host: the same two-generation choice as for vector-borne infections                                   |
| 14 | `14-contact-matrix`                | Six age groups with a contact matrix: numeric `R₀`, Perron vectors, elasticities, targeted vaccination                        |
| 15 | `15-plain-functions`               | Models as plain Julia functions via ForwardDiff; leaky versus all-or-nothing vaccines                                         |
| 16 | `16-effective-reproduction-number` | `R_t` in closed form, at a state, along a solution, and as an observed variable of the system                                 |
