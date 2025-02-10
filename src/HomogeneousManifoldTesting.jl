module HomogeneousManifoldTesting

using DocStringExtensions

using LinearAlgebra
using Plots
using Random
using Statistics
using NearestNeighbors

using RCall

import Distributions: Normal, MvNormal, Multinomial, Chisq, quantile, mean, std
import Manifolds: Grassmann, rand, distance

export RandomUniform
export plot_data

include("abstractdata.jl")
include("abstracttest.jl")

include("bolzadata.jl")
include("circledata.jl")
include("spheredata.jl")
include("torusdata.jl")
include("grassmanndata.jl")
include("abstractdatalist.jl")

include("test_homogeneity.jl")
include("test_homogeneity_param_selection.jl")
include("test_homogeneity_aggregated.jl")

include("compute_true_signatures.jl")
include("compute_mean_vect.jl")
include("distance_signature_test_.jl")
include("debiased_dtm_signature.jl")

include("plot_barycenters_true_signatures.jl")

include("compare_power.jl")

include("parameter_selection.jl")

end
