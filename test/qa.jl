using Aqua, ReproductiveNumbers, Test

Aqua.test_all(ReproductiveNumbers;
    ambiguities = false,       # ambiguities in upstream symbolic packages
    piracies = false,          # no piracy in this package; upstream noise
    stale_deps = (ignore = [:DocStringExtensions],),
    deps_compat = true,
    persistent_tasks = false)
Aqua.test_ambiguities(ReproductiveNumbers)
