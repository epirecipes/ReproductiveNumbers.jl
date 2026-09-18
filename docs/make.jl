using Documenter
using ReproductiveNumbers
using ModelingToolkit, Catalyst, Symbolics, LinearAlgebra

# Bundle the rendered Quarto vignettes and the Lean proof report into the site by copying
# them into `docs/src` (Documenter copies every file there to the build). The vignettes are
# rendered separately (`quarto render` in `vignettes/`, see the Documentation workflow); if
# they are absent the links on the vignettes page are reported as invalid and left dead.
const SRC = joinpath(@__DIR__, "src")
const VIGNETTES = normpath(joinpath(@__DIR__, "..", "vignettes"))
let copied = 0
    for dir in filter(d -> occursin(r"^\d\d-", d), readdir(VIGNETTES))
        for f in filter(f -> endswith(f, ".html"), readdir(joinpath(VIGNETTES, dir)))
            dest = joinpath(SRC, "vignettes", dir)
            mkpath(dest)
            cp(joinpath(VIGNETTES, dir, f), joinpath(dest, f); force = true)
            copied += 1
        end
    end
    copied == 0 &&
        @warn "no rendered vignettes found in $(VIGNETTES); run `quarto render` there first"
    @info "bundled $(copied) rendered vignette(s)"
end
let proofs = normpath(joinpath(
        @__DIR__, "..", "proofs", "docs", "ReproductiveNumbersProofs.html"))
    if isfile(proofs)
        mkpath(joinpath(SRC, "proofs"))
        cp(proofs, joinpath(SRC, "proofs", "ReproductiveNumbersProofs.html"); force = true)
    else
        @warn "rendered Lean report not found at $(proofs)"
    end
end

DocMeta.setdocmeta!(ReproductiveNumbers, :DocTestSetup,
    :(using ReproductiveNumbers); recursive = true)

makedocs(;
    sitename = "ReproductiveNumbers.jl",
    modules = [ReproductiveNumbers],
    authors = "Simon D. W. Frost and contributors",
    repo = Remotes.GitHub("epirecipes", "ReproductiveNumbers.jl"),
    format = Documenter.HTML(;
        prettyurls = get(ENV, "CI", nothing) == "true",
        edit_link = "main",
        canonical = "https://epirecip.es/ReproductiveNumbers.jl",
        assets = String[],
        size_threshold = 400 * 1024,
        size_threshold_warn = 200 * 1024),
    pages = [
        "Home" => "index.md",
        "Mathematical background" => "background.md",
        "Tutorial" => "tutorial.md",
        "Choosing what counts as a transmission" => "transmissions.md",
        "Vignettes" => "vignettes.md",
        "Formal proofs" => "proofs.md",
        "API reference" => "api.md"
    ],
    checkdocs = :exports,
    warnonly = [:cross_references],
    linkcheck = false)

deploydocs(;
    repo = "github.com/epirecipes/ReproductiveNumbers.jl",
    devbranch = "main",
    push_preview = true)
