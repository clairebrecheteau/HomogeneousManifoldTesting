
export Test_Homogeneity_Aggregated

"""
$(TYPEDEF)

Test_Homogeneity to test H_0"The samples of type AbstractData are samples from the uniform distribution on the corresponding homogeneous space".
"""


struct Test_Homogeneity_Aggregated <: AbstractTest

    np::Integer
    data_type::UnionAll
    proportion_nearest_neighbours::Vector{<:Real}
    weights::Vector{<:Real}
    monte_carlo_size::Integer
    monte_carlo_size_2::Integer
    vect_alpha::Vector{<:Real}
    vect_quantiles_homogeneity::Vector{Vector{<:Real}}
    vect_quantiles_iidness::Vector{Vector{<:Real}}
    test_0::Test_Homogeneity

end


function Test_Homogeneity_Aggregated(
    np::Integer,
    proportion_nearest_neighbours::Vector{<:Real},
    weights::Vector{<:Real},
    monte_carlo_size::Integer,
    monte_carlo_size_2::Integer,
    data_type::UnionAll,
    vect_alpha::Vector{<:Real},
    error_approx_u = 0.001
)

    @assert np>0 "The sample size argument np should be non negative."
    @assert data_type in AbstractDataList "The data_type should belong to [CircleData,SphereData,TorusData,BolzaData,GrassmannData23,GrassmannData24]."
    @assert all((proportion_nearest_neighbours.>0) .&& (proportion_nearest_neighbours.<=1)) "proportion_nearest_neighbours is a vector with elements in (0,1]."
    @assert monte_carlo_size>0 "monte_carlo_size argument should be non negative."
    @assert size(weights)[1] == size(proportion_nearest_neighbours)[1] "The vector weights should have the same length as the vector proportion_nearest_neighbours."

    test_0 = Test_Homogeneity(np,proportion_nearest_neighbours,monte_carlo_size,data_type)

    n_params = size(proportion_nearest_neighbours)[1]

        # ATTENTION !!!
        # Besoin d'un nouvel échantillon pas trié
        # statistics_2, taille monte_carlo_size_2 

        

    statistics_2_homogeneity = [Vector{Float64}(undef,monte_carlo_size_2) for _ in 1:n_params]
    statistics_2_iidness = [Vector{Float64}(undef,monte_carlo_size_2) for _ in 1:n_params]

    for n_ in 1:monte_carlo_size_2
        data = data_type(np,"uniform")
        signatures = debiased_dtm_signature(squared_distance_matrix(data),proportion_nearest_neighbours)

        statistic = zeros(n_params)

        for j in 1:n_params
            statistics_2_homogeneity[j][n_] = distance_signatures_test_homogeneity(signatures[j,:], test_0.true_signatures[j],np)
        end

        for j in 1:n_params
            statistics_2_iidness[j][n_] = distance_signatures_test_iidness(sort(signatures[j,:]), test_0.barycenter_signatures[j],np)
        end

    end

    statistics_2 = [statistics_2_homogeneity, statistics_2_iidness]

    # Compute the value of u_α

    quant_type = eltype(test_0.stats_null_hypothesis_homogeneity[1])

    vect_quantiles_homogeneity = Vector{Vector{quant_type}}(undef,size(vect_alpha)[1])
    vect_quantiles_iidness = Vector{Vector{quant_type}}(undef,size(vect_alpha)[1])

    for (idx_alpha, alpha) in enumerate(vect_alpha)
        quantiles = Vector{Vector{quant_type}}(undef,2)

        for (idx_,statistics) in enumerate([test_0.stats_null_hypothesis_homogeneity, test_0.stats_null_hypothesis_iidness])

            Δ = 1
            u_min = 0
            u_max = 1
            u = mean([u_min,u_max])
            while Δ > error_approx_u
                centered_stats = hcat([ statistics_2[idx_][idx] .- stat[Int(ceil(monte_carlo_size*(1 - u*weights[idx])))] for (idx,stat) in enumerate(statistics)]...) 
                T = mean(maximum(centered_stats, dims = 2) .> 0)
                if T > alpha
                    u_max = u
                else
                    u_min = u
                end
                Δ = u_max - u_min
                u = mean([u_min,u_max])
            end

            quantiles[idx_] = [stat[Int(ceil(monte_carlo_size*(1 - u*weights[idx])))] for (idx, stat) in enumerate(statistics)]
        end
        
        vect_quantiles_homogeneity[idx_alpha] = quantiles[1]
        vect_quantiles_iidness[idx_alpha] = quantiles[2]
    end

    return Test_Homogeneity_Aggregated(np,data_type,proportion_nearest_neighbours,weights,monte_carlo_size,monte_carlo_size_2,vect_alpha,vect_quantiles_homogeneity,vect_quantiles_iidness,test_0)

end


export apply_test

for datatype in [:CircleData, :SphereData, :TorusData, :BolzaData, :GrassmannData23, :GrassmannData24]

    @eval function apply_test(
        t::Test_Homogeneity_Aggregated,
        data::$datatype,
        method::String
    ) 

        @assert method in ["homogeneity","iidness"] "The method is either 'homogeneity' or 'iidness'."
        @assert t.data_type == $datatype "The data type should be the same as the data type for the test."
        @assert t.np == data.np "The sample data should be of the same size as the test method."
   
        signatures = debiased_dtm_signature(squared_distance_matrix(data),t.proportion_nearest_neighbours)

        n_params = size(t.proportion_nearest_neighbours)[1]

        statistic = zeros(n_params)

        vect_stat = Vector{eltype(t.vect_quantiles_homogeneity[1])}(undef,size(t.vect_alpha)[1])
        vect_param_max = Vector{eltype(t.proportion_nearest_neighbours)}(undef,size(t.vect_alpha)[1])

        if method == "homogeneity"
            for j in 1:n_params
                statistic[j] = distance_signatures_test_homogeneity(signatures[j,:], t.test_0.true_signatures[j],t.test_0.np)
            end
            for (idx_alpha, quant_hom) in enumerate(t.vect_quantiles_homogeneity)
                diff_quantile = [st - quant_hom[idx] for (idx, st) in enumerate(statistic)]
            #diff_quantile = [st - t.quantiles_homogeneity[idx] for (idx, st) in enumerate(statistic)]
                stat, index_max = findmax(diff_quantile)
                param_max = t.proportion_nearest_neighbours[index_max]
                vect_stat[idx_alpha] = stat
                vect_param_max[idx_alpha] = param_max
            end
        elseif method == "iidness"
            for j in 1:n_params
                statistic[j] = distance_signatures_test_iidness(sort(signatures[j,:]), t.test_0.barycenter_signatures[j],t.test_0.np)
            end
            for (idx_alpha, quant_iid) in enumerate(t.vect_quantiles_iidness)
                diff_quantile = [st - quant_iid[idx] for (idx, st) in enumerate(statistic)]
                stat, index_max = findmax(diff_quantile)
                param_max = t.proportion_nearest_neighbours[index_max]
                vect_stat[idx_alpha] = stat
                vect_param_max[idx_alpha] = param_max
            end
        end

        return [vect_stat, vect_param_max]

    end

end
