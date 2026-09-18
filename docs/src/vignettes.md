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
|:-- |:-------------------------------------------------------------------------------------------------- |:--------------------------------------------------------------------------------------------------------------------------------- |
| 1 | [The SIR model](vignettes/01-sir/sir.html)                                                         | `T`, `Σ`, `K` and `R₀ = β/γ`; numeric evaluation; the threshold by simulation; a reproduction number at any state                 |
| 2 | [SEIR with demography](vignettes/02-seir/seir.html)                                                | Large versus restricted domain, states-at-infection, automatic infection-free steady state, characteristic polynomial             |
| 3 | [Two latent stages](vignettes/03-two-latent-stages/two-latent-stages.html)                         | section 2.1 of Diekmann et al.: a rank-one `K` and the small-domain matrix                                                        |
| 4 | [Vector-borne transmission](vignettes/04-vector-borne/vector-borne.html)                           | Ross–Macdonald: `R₀ = √(…)` and what changes when only host-to-host generations are counted                                       |
| 5 | [Reaction networks with Catalyst](vignettes/05-catalyst/catalyst.html)                             | Transmission reactions identified from stoichiometry; vertical transmission; overriding by predicate, flags or indices            |
| 6 | [Heterogeneous mixing](vignettes/06-heterogeneous-mixing/heterogeneous-mixing.html)                | Two host types (`2 × 2` closed form), separable mixing, multi-strain block structure, an age-structured model that needs numerics |
| 7 | [Type reproduction numbers](vignettes/07-type-reproduction-numbers/type-reproduction-numbers.html) | Targeted control of a core group; type reproduction numbers in the vector-borne model                                             |
| 8 | [Vertical and sexual transmission](vignettes/08-vertical-transmission-sti/vertical-transmission-sti.html) | Section 4.2 of Diekmann et al.: `K` has no closed form but the small-domain matrix `K_S` does |
| 9 | [Bovine viral diarrhoea](vignettes/09-bovine-viral-diarrhoea/bovine-viral-diarrhoea.html) | Section 4.3 of Diekmann et al.: persistently infected animals as a second state-at-infection |
| 10 | [Staged progression](vignettes/10-staged-progression/staged-progression.html) | `K` stays `1 × 1` while `K_L` grows; elasticities by stage |
| 11 | [Multiple strains](vignettes/11-multi-strain/multi-strain.html) | Cross-immunity, invasion reproduction numbers, mutation and co-infection |
| 12 | [Waning immunity](vignettes/12-waning-immunity/waning-immunity.html) | SEIRS with reinfection: `R₀` checked against steady states from SteadyStateDiffEq and NonlinearSolve |
| 13 | [Within-host dynamics](vignettes/13-within-host/within-host.html) | Viral dynamics within a host: the same two-generation choice as for vector-borne infections |
| 14 | [Contact matrix](vignettes/14-contact-matrix/contact-matrix.html) | Six age groups: numeric `R₀`, Perron vectors, elasticities, targeted vaccination |
| 15 | [Plain functions and vaccines](vignettes/15-plain-functions/plain-functions.html) | Models as plain Julia functions via ForwardDiff; leaky versus all-or-nothing vaccines |
| 16 | [The effective reproduction number](vignettes/16-effective-reproduction-number/effective-reproduction-number.html) | `R_t` in closed form, at a state, along a solution, and as an observed variable of the system |
| 17 | [Seasonal forcing](vignettes/17-seasonal-forcing/seasonal-forcing.html) | Seasonally forced SIR with births and deaths: `R_t` oscillating about one on the endemic attractor, with an intervention |
