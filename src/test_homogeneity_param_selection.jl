export Test_Homogeneity_Param_Selection

"""
$(TYPEDEF)

Test_Homogeneity_Param_Selection to test H_0"The samples of type AbstractData are samples from the uniform distribution on the corresponding homogeneous space".
    - np : number of sample points
    - nb_nearest_neighbours : vector of number of nearest neighbours : parameters for the dtm signature. 
    - monte_carlo_size : number of signatures of iid uniform samples to compute the barycenter and to estimate quantiles of the statistic under hypothesis H_0.
    - barycenter_signatures : Wasserstein barycenters of the dtm-signatures with parameters in nb_nearest_neighbours.
    - stats_null_hypothesis :
"""


struct Test_Homogeneity_Param_Selection <: AbstractTest

    np::Integer
    data_type::UnionAll
    proportion_nearest_neighbours::Vector{<:Real}
    monte_carlo_size::Integer
    barycenter_signatures::Vector{Matrix{<:Real}}
    true_signatures::Vector{<:Real}
    samples_signatures_homogeneity::Vector{<:Vector{<:Matrix{<:Real}}}
    samples_signatures_iidness::Vector{<:Vector{<:Matrix{<:Real}}}
    stats_null_hypothesis_homogeneity::Vector{<:Vector{<:Real}}
    stats_null_hypothesis_iidness::Vector{<:Vector{<:Real}}
    test_0::Test_Homogeneity

end

function Test_Homogeneity_Param_Selection(
    np::Integer,
    proportion_nearest_neighbours::Vector{<:Real},
    monte_carlo_size::Integer,
    data_type::UnionAll
)

    @assert np>0 "The sample size argument np should be non negative."
    @assert data_type in AbstractDataList "The data_type should belong to [CircleData,SphereData,TorusData,BolzaData]."
    @assert all((proportion_nearest_neighbours.>0) .&& (proportion_nearest_neighbours.<=1)) "proportion_nearest_neighbours is a vector with elements in (0,1]."
    @assert monte_carlo_size>0 "monte_carlo_size argument should be non negative."

    test_0 = Test_Homogeneity(np,proportion_nearest_neighbours,monte_carlo_size,data_type)

    # fonction compute_mean_vect à récrire avec l'argument test_0 à la fin
    # barycenter_signatures vector of n_params vectors (barycenter signatures). The test will compare to the correct barycenter.
    # samples_signatures is a vector of n_params vectors of monte_carlo_size signatures. 
    barycenter_signatures, samples_signatures_homogeneity, samples_signatures_iidness = compute_mean_vect(data_type,[np],proportion_nearest_neighbours,monte_carlo_size,test_0)

    barycenter_signatures = barycenter_signatures[1]
    true_signatures = compute_true_signatures(data_type,proportion_nearest_neighbours)

    n_params = size(proportion_nearest_neighbours)[1]

    stats_null_hypothesis_homogeneity = [sort(distance_signatures_test_homogeneity(samples_signatures_homogeneity[1][j], true_signatures[j],np)) for j in 1:n_params]
    stats_null_hypothesis_iidness = [sort(distance_signatures_test_iidness(samples_signatures_iidness[1][j], barycenter_signatures[j],np)) for j in 1:n_params]

    
    return Test_Homogeneity_Param_Selection(np,data_type,proportion_nearest_neighbours,monte_carlo_size,barycenter_signatures,true_signatures,samples_signatures_homogeneity,samples_signatures_iidness,stats_null_hypothesis_homogeneity,stats_null_hypothesis_iidness,test_0)

end




export apply_test

for datatype in [:CircleData, :SphereData, :TorusData, :BolzaData, :GrassmannData23, :GrassmannData24]

    @eval function apply_test(
        t::Test_Homogeneity_Param_Selection,
        data::$datatype,
        method::String
    ) 

        @assert method in ["homogeneity","iidness"] "The method is either 'homogeneity' or 'iidness'."
        @assert t.data_type == $datatype "The data type should be the same as the data type for the test."
        @assert t.np == data.np "The sample data should be of the same size as the test method."

        if method == "homogeneity"
            index = argmin(apply_test(t.test_0, data, "homogeneity")[1])
        elseif method == "iidness"
            index = argmin(apply_test(t.test_0, data,"iidness")[1])
        end

        signatures = debiased_dtm_signature(squared_distance_matrix(data),[t.proportion_nearest_neighbours[index]])

        #n_params = size(t.proportion_nearest_neighbours)[1]
        statistic = 0
        p_value = 0

        if method == "homogeneity"
            statistic = distance_signatures_test_homogeneity(signatures[1,:], t.true_signatures[index],t.np)
            p_value = mean([stat >= statistic for stat in t.stats_null_hypothesis_homogeneity[index]])
        elseif method == "iidness"
            statistic = distance_signatures_test_iidness(sort(signatures[1,:]), t.barycenter_signatures[index],t.np)
            p_value = mean([stat >= statistic for stat in t.stats_null_hypothesis_iidness[index]])
        end

        return [p_value, statistic, index, t.proportion_nearest_neighbours[index]]

    end

end