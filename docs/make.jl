using Documenter
using HomogeneousManifoldTesting

ENV["GKSwstype"] = "100"

makedocs(
    sitename = "HomogeneousManifoldTesting",
    format = Documenter.HTML(),
    modules = [HomogeneousManifoldTesting],
    pages=[
        "Home" => "index.md",
        "Sample generation" => "samples_generation.md",
        "Signature computation generation" => "signatures_computation.md"
    ]
)

deploydocs(;
    repo="github.com/clairebrecheteau/HomogeneousManifoldTesting",
)
