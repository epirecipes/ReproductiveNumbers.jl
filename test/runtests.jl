using SafeTestsets, Test

@testset "ReproductiveNumbers.jl" begin
    @safetestset "Quality assurance (Aqua)" include("qa.jl")
    @safetestset "Unit: symbolic utilities" include("unit/test_utils.jl")
    @safetestset "Unit: spectral radius" include("unit/test_spectral.jl")
    @safetestset "Unit: decomposition" include("unit/test_decomposition.jl")
    @safetestset "Unit: analysis" include("unit/test_analysis.jl")
    @safetestset "Unit: routing" include("unit/test_routing.jl")
    @safetestset "Integration: ModelingToolkit models" include("integration/test_modelingtoolkit.jl")
    @safetestset "Integration: Catalyst models" include("integration/test_catalyst.jl")
    @safetestset "Integration: numeric evaluation" include("integration/test_numeric.jl")
    @safetestset "Integration: features" include("integration/test_features.jl")
    @safetestset "Regression: literature formulae" include("regression/test_literature.jl")
    @safetestset "Regression: threshold and simulation" include("regression/test_threshold.jl")
    @safetestset "Regression: steady states" include("regression/test_steady_state.jl")
    @safetestset "Regression: paper examples" include("regression/test_paper_examples.jl")
end
