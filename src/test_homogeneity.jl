export Test_Homogeneity

"""
$(TYPEDEF)

Test_Homogeneity to test H_0"The samples of type AbstractData are samples from the uniform distribution on the corresponding homogeneous space".
    - np : number of sample points
    - nb_nearest_neighbours : vector of number of nearest neighbours : parameters for the dtm signature. 
    - monte_carlo_size : number of signatures of iid uniform samples to compute the barycenter and to estimate quantiles of the statistic under hypothesis H_0.
    - barycenter_signatures : Wasserstein barycenters of the dtm-signatures with parameters in nb_nearest_neighbours.
    - stats_null_hypothesis :
"""


struct Test_Homogeneity <: AbstractTest

    np::Integer
    data_type::UnionAll
    proportion_nearest_neighbours::Vector{<:Real}
    monte_carlo_size::Integer
    barycenter_signatures::Vector{Matrix{<:Real}}
    true_signatures::Vector{<:Real}
    samples_signatures::Vector{<:Vector{<:Matrix{<:Real}}}
    stats_null_hypothesis_homogeneity::Vector{<:Vector{<:Real}}
    stats_null_hypothesis_iidness::Vector{<:Vector{<:Real}}

end

function Test_Homogeneity(
    np::Integer,
    proportion_nearest_neighbours::Vector{<:Real},
    monte_carlo_size::Integer,
    data_type::UnionAll
)

    @assert np>0 "The sample size argument np should be non negative."
    @assert data_type in AbstractDataList "The data_type should belong to [CircleData,SphereData,TorusData,BolzaData,GrassmannData23,GrassmannData24]."
    @assert all((proportion_nearest_neighbours.>0) .&& (proportion_nearest_neighbours.<=1)) "proportion_nearest_neighbours is a vector with elements in (0,1]."
    @assert monte_carlo_size>0 "monte_carlo_size argument should be non negative."

    barycenter_signatures, samples_signatures = compute_mean_vect(data_type,[np],proportion_nearest_neighbours,monte_carlo_size)

    barycenter_signatures = barycenter_signatures[1]
    true_signatures = compute_true_signatures(data_type,proportion_nearest_neighbours)

    n_params = size(proportion_nearest_neighbours)[1]

    stats_null_hypothesis_homogeneity = [sort(distance_signatures_test_homogeneity(samples_signatures[1][j], true_signatures[j],np)) for j in 1:n_params]
    stats_null_hypothesis_iidness = [sort(distance_signatures_test_iidness(samples_signatures[1][j], barycenter_signatures[j],np)) for j in 1:n_params]

    
    return Test_Homogeneity(np,data_type,proportion_nearest_neighbours,monte_carlo_size,barycenter_signatures,true_signatures,samples_signatures,stats_null_hypothesis_homogeneity,stats_null_hypothesis_iidness)

end




export apply_test

for datatype in [:CircleData, :SphereData, :TorusData, :BolzaData, :GrassmannData23, :GrassmannData24]

    @eval function apply_test(
        t::Test_Homogeneity,
        data::$datatype,
        method::String
    ) 

        @assert method in ["homogeneity","iidness"] "The method is either 'homogeneity' or 'iidness'."
        @assert t.data_type == $datatype "The data type should be the same as the data type for the test."
        @assert t.np == data.np "The sample data should be of the same size as the test method."

        signatures = debiased_dtm_signature(squared_distance_matrix(data),t.proportion_nearest_neighbours)
        p_value = [0]

        n_params = size(t.proportion_nearest_neighbours)[1]
        statistic = zeros(n_params)
        p_value = zeros(n_params)

        if method == "homogeneity"
            for j in 1:n_params
                statistic[j] = distance_signatures_test_homogeneity(signatures[j,:], t.true_signatures[j],t.np)
                p_value[j] = mean([stat >= statistic[j] for stat in t.stats_null_hypothesis_homogeneity[j]])
            end
        elseif method == "iidness"
            for j in 1:n_params
                statistic[j] = distance_signatures_test_iidness(sort(signatures[j,:]), t.barycenter_signatures[j],t.np)
                p_value[j] = mean([stat >= statistic[j] for stat in t.stats_null_hypothesis_iidness[j]])
            end
        end

        return [p_value, statistic]

    end

end