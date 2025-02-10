export compute_power

function compute_power(
    data_type::UnionAll,
    method::String,
    vect_params,
    sample_size_vect::Vector{<:Integer},
    proportion_nearest_neighbours::Vector{<:Real},
    monte_carlo_size::Integer,
    n_rep::Integer
)

    n_params = size(vect_params)[1]
    n_sample_sizes = size(sample_size_vect)[1]
    n_prop_nn = size(proportion_nearest_neighbours)[1]

    samples = [[[data_type(np,method,params...) for n in 1:n_rep] for params in vect_params] for np in sample_size_vect]

    ans_hom = 0
    ans_iid = 0

    rejection_hom = Vector{Vector{Vector{Real}}}(undef,n_sample_sizes)
    rejection_iid = Vector{Vector{Vector{Real}}}(undef,n_sample_sizes)

    for i in 1:n_sample_sizes
        np = sample_size_vect[i]
        test = Test_Homogeneity(np,proportion_nearest_neighbours,monte_carlo_size,data_type)
        rejection_r_hom = Vector{Vector{Real}}(undef,n_params)
        rejection_r_iid = Vector{Vector{Real}}(undef,n_params)
        for j in 1:n_params
            print(j)
            n_rejected_hom = zeros(n_prop_nn)
            n_rejected_iid = zeros(n_prop_nn)
            for k in 1:n_rep
                ans_hom = apply_test(test,samples[i][j][k],"homogeneity")
                ans_iid = apply_test(test,samples[i][j][k],"iidness")
                n_rejected_hom .+= (ans_hom[1].<= alpha)
                n_rejected_iid .+= (ans_iid[1].<= alpha)
            end
            rejection_r_hom[j] = n_rejected_hom./n_rep
            rejection_r_iid[j] = n_rejected_iid./n_rep
        end
        rejection_hom[i] = rejection_r_hom
        rejection_iid[i] = rejection_r_iid
    end

    return rejection_hom, rejection_iid

end

export plot_power

function plot_power(
    rejection_hom::Vector{<:Real},
    rejection_iid::Vector{<:Real},
    vect_params::Vector{<:Real},
    data_type::UnionAll,


)

plot(vect_params,[vcat(rejection_hom[1]...),vcat(rejection_iid[1]...)],label= ["test homogeneity" "test iidness"], title = "$(data_type), α = 0.05, sample size = $(sample_size_vect[1])",xlab = "r",ylab = "power",ylims = [0,0.3])

end