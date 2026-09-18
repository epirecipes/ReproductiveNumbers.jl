module ReproductiveNumbersLatexifyExt

using ReproductiveNumbers
using ReproductiveNumbers: NextGenerationMatrix, isnumeric
using Latexify

"""
    latexify(ngm::NextGenerationMatrix)

LaTeX rendering of a next-generation matrix decomposition: the infected states, the
transmission matrix `T`, the transition matrix `Σ`, the next-generation matrix `K` with
its states-at-infection, and `R₀` when a closed form exists.
"""
function Latexify.latexify(ngm::NextGenerationMatrix; kwargs...)
    mat(A) = Latexify.latexraw(A; kwargs...)
    parts = String[]
    push!(parts, "x = " * mat(ngm.infected))
    push!(parts, "T = " * mat(ngm.T))
    push!(parts, "\\Sigma = " * mat(ngm.Σ))
    push!(parts,
        "K = " * mat(ngm.K) * " \\quad \\text{indexed by } " * mat(ngm.states_at_infection))
    R0 = try
        basic_reproduction_number(ngm)
    catch err
        err isa NoClosedFormError ? nothing : rethrow()
    end
    R0 === nothing || push!(parts, "R_0 = " * Latexify.latexraw(R0; kwargs...))
    # a single line: renderers that treat the string as Markdown (Documenter) would
    # otherwise read lines beginning with `-` as list items
    body = replace(join(parts, " \\\\ "), "\n" => " ")
    return Latexify.LaTeXString("\\begin{align*} " * body * " \\end{align*}")
end

function Base.show(io::IO, ::MIME"text/latex", ngm::NextGenerationMatrix)
    print(io, Latexify.latexify(ngm))
end

end
