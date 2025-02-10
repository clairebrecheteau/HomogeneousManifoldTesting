function compute_mean_vect(
    data_type ::UnionAll,
    sample_size_vect::Vector{<:Integer},
    proportion_nearest_neighbours::Vector{<:Real},
    monte_carlo_size::Integer
)

    mean_vect = Vector{Vector{Matrix{Float64}}}(undef, size(sample_size_vect)[1])
    samples_signatures = Vector{Vector{Matrix{Float64}}}(undef, size(sample_size_vect)[1])
    nthreads = Threads.nthreads()
    n_per_threads = Int.(floor.(ones(nthreads).*monte_carlo_size/nthreads))
    it0 = 1

    while sum(n_per_threads)<monte_carlo_size
        n_per_threads[it0] +=1
        it0 +=1
    end
    result =  Vector{Vector{Matrix{Float64}}}(undef, nthreads)

    
    for i in 1:(size(sample_size_vect)[1])
        print(i)
        sample_size = sample_size_vect[i]
        Threads.@threads for it in 1:nthreads
            result[it] = [debiased_dtm_signature(squared_distance_matrix(data_type(sample_size,"uniform")),proportion_nearest_neighbours) for i in 1:n_per_threads[it]]
        end
        dtm_sigs = vcat(result...)
        v   = [hcat([sort(dtm_sig[j,:]) for dtm_sig in dtm_sigs]...) for j in 1:size(proportion_nearest_neighbours)[1]] # All signatures are sorted
        mean_dtm_sig = [mean(v[j],dims = 2) for j in 1:size(proportion_nearest_neighbours)[1]]
        mean_vect[i] = mean_dtm_sig
        samples_signatures[i] = v
    end

    return [mean_vect, samples_signatures]

end


function compute_mean_vect(
    data_type ::UnionAll,
    sample_size_vect::Vector{<:Integer},
    proportion_nearest_neighbours::Vector{<:Real},
    monte_carlo_size::Integer,
    test_0::Test_Homogeneity
)

    n_params = size(proportion_nearest_neighbours)[1]
    monte_carlo_size = monte_carlo_size*n_params
    # We expect enough signatures (around initial monte_carlo_size) per parameter (prop. nn).
    # Contrary to classical compute_mean_vect function, a sample is used for only one parameter, the one for which the pvalue of test test_0 is the smallest.

    mean_vect = Vector{Vector{Matrix{Float64}}}(undef, size(sample_size_vect)[1])
    samples_signatures_homogeneity = Vector{Vector{Matrix{Float64}}}(undef, size(sample_size_vect)[1])
    samples_signatures_iidness = Vector{Vector{Matrix{Float64}}}(undef, size(sample_size_vect)[1])

    nthreads = Threads.nthreads()
    n_per_threads = Int.(floor.(ones(nthreads).*monte_carlo_size/nthreads))
    it0 = 1

    while sum(n_per_threads)<monte_carlo_size
        n_per_threads[it0] +=1
        it0 +=1
    end
    result =  Vector{Vector{Matrix{Float64}}}(undef, nthreads)
    indices_min_pval_homogeneity = Vector{Vector{Integer}}(undef, nthreads)
    indices_min_pval_iidness = Vector{Vector{Integer}}(undef, nthreads)
    
    for i in 1:(size(sample_size_vect)[1])
        sample_size = sample_size_vect[i]
        Threads.@threads for it in 1:nthreads
            result[it] = Vector{Matrix{Float64}}(undef,n_per_threads[it])
            indices_min_pval_homogeneity[it] = Vector{Integer}(undef,n_per_threads[it])
            indices_min_pval_iidness[it] = Vector{Integer}(undef,n_per_threads[it])
            for j in 1:n_per_threads[it]
                sample = data_type(sample_size,"uniform")
                index_min_pval_homogeneity = argmin(apply_test(test_0,sample,"homogeneity")[1])
                index_min_pval_iidness = argmin(apply_test(test_0,sample,"iidness")[1])
                indices_min_pval_homogeneity[it][j] = index_min_pval_homogeneity
                indices_min_pval_iidness[it][j] = index_min_pval_iidness
                result[it][j] = debiased_dtm_signature(squared_distance_matrix(sample),proportion_nearest_neighbours)
            end
        end

        indices_homogeneity = vcat(indices_min_pval_homogeneity...)
        indices_iidness = vcat(indices_min_pval_iidness...)
        dtm_sigs = vcat(result...)

        # Check that there are more than one point in each group ????

        v   = [hcat([sort(dtm_sigs[i][j,:]) for i in 1:size(dtm_sigs)[1] if indices_iidness[i] == j]...) for j in 1:size(proportion_nearest_neighbours)[1]] # All signatures are sorted
        mean_dtm_sig = [mean(v[j],dims = 2) for j in 1:size(proportion_nearest_neighbours)[1]]
        mean_vect[i] = mean_dtm_sig
        samples_signatures_iidness[i] = v

        w   = [hcat([sort(dtm_sigs[i][j,:]) for i in 1:size(dtm_sigs)[1] if indices_homogeneity[i] == j]...) for j in 1:size(proportion_nearest_neighbours)[1]] 
        samples_signatures_homogeneity[i] = w

    end

    return [mean_vect, samples_signatures_homogeneity, samples_signatures_iidness]

end