
export compare_power


for datatype in [:CircleData, :SphereData, :TorusData, :BolzaData]

    @eval function compare_power(
        samples::Vector{<:$datatype},
        t::AbstractTest,
        α::Vector{<:Real}
    )

        data_type = t.data_type
        sample_size = t.np

        @assert data_type == $datatype "The test and the sample should have the same data type."
        @assert sample_size == samples[1].np "The test and the data should share the same parameter np (sample size)."

        prop_near_neigh = t.proportion_nearest_neighbours
        n_params = size(prop_near_neigh)[1]

        tests_homogeneity = [apply_test(t,a,"homogeneity") for a in samples];
        tests_iidness = [apply_test(t,a,"iidness") for a in samples];

        if typeof(t) == Test_Homogeneity
            pval_homogeneity = [ [t[1][n] for t in tests_homogeneity] for n in 1:n_params]
            pval_iidness = [ [t[1][n] for t in tests_iidness] for n in 1:n_params]
        else
            pval_homogeneity = [t_[1] for t_ in tests_homogeneity]
            pval_iidness = [t_[1] for t_ in tests_iidness]
        end

        samples_ = [a.points for a in samples];

        p = plot_data(samples[1])

        if typeof(t) == Test_Homogeneity
            power_homogeneity = [[mean(pval_homogeneity[n] .< alph) for alph in α] for n in 1:n_params]
            power_iidness = [[mean(pval_iidness[n] .< alph) for alph in α] for n in 1:n_params]
        else
            power_homogeneity = [mean(pval_homogeneity .< alph) for alph in α]
            power_iidness = [mean(pval_iidness .< alph) for alph in α]
        end
        
    return p, samples_, pval_homogeneity, pval_iidness, power_homogeneity, power_iidness
    end

end

function compare_power(
    t::AbstractTest,
    n_rep::Integer,
    α::Vector{<:Real},
    sample_method::String,
    params...
)

    data_type = t.data_type
    sample_size = t.np
    prop_near_neigh = t.proportion_nearest_neighbours
    n_params = size(prop_near_neigh)[1]

    samples = [data_type(sample_size,sample_method,params...) for n in 1:n_rep];

    return compare_power(samples,t,α)

    #tests_homogeneity = [apply_test(t,a,"homogeneity") for a in samples];
    #tests_iidness = [apply_test(t,a,"iidness") for a in samples];

    #pval_homogeneity = [ [t[1][n] for t in tests_homogeneity] for n in 1:n_params]
    #pval_iidness = [ [t[1][n] for t in tests_iidness] for n in 1:n_params]

    #samples_ = [a.points for a in samples];

    #p = plot_data(samples[1])

    #power_homogeneity = [[mean(pval_homogeneity[n] .< alph) for alph in α] for n in 1:n_params]
    #power_iidness = [[mean(pval_iidness[n] .< alph) for alph in α] for n in 1:n_params]

    #return p, samples_, pval_homogeneity, pval_iidness, power_homogeneity, power_iidness

end