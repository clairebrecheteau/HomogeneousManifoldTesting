export plot_barycenters_true_signatures

function plot_barycenters_true_signatures(
    data_type ::UnionAll,
    sample_size_vect::Vector{<:Integer},
    proportion_nearest_neighbours::Vector{<:Real},
    monte_carlo_size::Integer,
    title_::String,
    ylims_ = [0,2],
    palett = :rose,
    legend_ = true
)

    tests = [Test_Homogeneity(np,proportion_nearest_neighbours,monte_carlo_size,data_type) for np in sample_size_vect]

    true_signatures = tests[1].true_signatures
    
    h = proportion_nearest_neighbours[1]

    n_sample_size_vect = size(sample_size_vect)[1]
    n_proportion_nearest_neighbours = size(proportion_nearest_neighbours)[1]

    type_line = repeat([:dot,:dashdotdot,:dashdot,:dash],Int(ceil(n_sample_size_vect/4)))

    if legend_
        plot((1:sample_size_vect[1])/sample_size_vect[1],tests[1].barycenter_signatures[1],linewidth=2, ls = type_line[1], alpha = 1/n_sample_size_vect, label = "n = $(sample_size_vect[1])", line_z = h,
color = palett,clims=(0,0.5), ylims = ylims_, title = title_,colorbar = false)
        for j in 2:n_sample_size_vect
            plot!((1:sample_size_vect[j])/sample_size_vect[j],tests[j].barycenter_signatures[1],linewidth=2, ls = type_line[j], line_z = h, alpha = j/n_sample_size_vect, label = "n = $(sample_size_vect[j])",color = palett, clims=(0,0.5))
        end
        plot!([0,1],[true_signatures[1],true_signatures[1]],linewidth=3, line_z = h, label = "true signature",color = palett, clims=(0,0.5))
    else
        plot((1:sample_size_vect[1])/sample_size_vect[1],tests[1].barycenter_signatures[1],linewidth=2, ls = type_line[1], alpha = 1/n_sample_size_vect, ylims = ylims_, title = title_, label="",line_z = h, color = palett, clims=(0,0.5),colorbar = false)
        for j in 2:n_sample_size_vect
            plot!((1:sample_size_vect[j])/sample_size_vect[j],tests[j].barycenter_signatures[1],linewidth=2, ls = type_line[j], alpha = j/n_sample_size_vect, label="",line_z = h, color = palett, clims=(0,0.5))
        end
        plot!([0,1],[true_signatures[1],true_signatures[1]],linewidth=3, line_z = h, label="",color = palett, clims=(0,0.5))
    end

    for i in 1:n_proportion_nearest_neighbours
        i = n_proportion_nearest_neighbours - i + 1
        h = proportion_nearest_neighbours[i]
        for j in 1:n_sample_size_vect
            plot!((1:sample_size_vect[j])/sample_size_vect[j],tests[j].barycenter_signatures[i],linewidth=2, ls = type_line[j], alpha = j/n_sample_size_vect, label = "", line_z = h, color = palett, clims=(0,0.5))
        end
        plot!([0,1],[true_signatures[i],true_signatures[i]],linewidth=3,  label = "", line_z = h, color = palett, clims=(0,0.5))    
    end

    if legend_
        plot!([0,1],[0,0],linewidth=0, label = " ", color= "white")
        for i in 1:n_proportion_nearest_neighbours
            i = n_proportion_nearest_neighbours - i + 1
            h = proportion_nearest_neighbours[i]
            plot!([0,1],[true_signatures[i],true_signatures[i]],linewidth=2, line_z = h, label = "h = $(h)",color = palett, clims=(0,0.5))
        end
        plot!(legend=:topleft, legendcolumns=1)
    end
    plot!([0,1],[0,0],linewidth=0, label = "", color= "white")
end


export plot_mean_vect_samples

function plot_mean_vect_samples(
    test::Test_Homogeneity,
    n_sample_signatures::Integer,
    title_:: String,
    ylims_ = [0,2]::Vector{<:Real},
    ls_ = :dot,
    palett = :rose
)

    @assert n_sample_signatures <= test.monte_carlo_size "We cannot plot more signatures, n_sample_signatures, than monte_carlo_size"

    n_proportion_nearest_neighbours = size(test.proportion_nearest_neighbours)[1]
    
    h = test.proportion_nearest_neighbours[n_proportion_nearest_neighbours]
    plot((1:test.np)/test.np,test.barycenter_signatures[n_proportion_nearest_neighbours],linewidth=3, ls = ls_, alpha = 1, label = "barycenter, h = $(h)",title = title_, ylims = ylims_,line_z = test.proportion_nearest_neighbours[n_proportion_nearest_neighbours], color = palett, clims=(0,0.5),colorbar = false)

    for i in 2:n_proportion_nearest_neighbours
        i = n_proportion_nearest_neighbours - i + 1
        h = test.proportion_nearest_neighbours[i]
        plot!((1:test.np)/test.np,test.barycenter_signatures[i],linewidth=3, ls = ls_, alpha = 1, label = "barycenter, h = $(h)",title = title_, ylims = ylims_,line_z = test.proportion_nearest_neighbours[i], color = palett, clims=(0,0.5),colorbar = false)
    end
    print("ok")

    for i in 1:n_proportion_nearest_neighbours
        h = test.proportion_nearest_neighbours[i]
        for iter in 1:n_sample_signatures
            plot!((1:test.np)/test.np,test.samples_signatures[1][i][:,iter],linewidth=1, alpha = 0.2,line_z = h, color = palett, clims=(0,0.5), label = "")
        end
    end
    plot!([0,0],[0,1],linewidth=0, label = "", color= "white")
    plot!(legend=:topleft)
end

"""

# Plot the quantile functions - black and white
function plot_mean_vect(mean_vect,true_signatures,sample_size_vect, title_, ylims_ = [0,2], legend_ = true)
    h = h_vect[5]

    if legend_
        plot((1:sample_size_vect[1])/sample_size_vect[1],mean_vect[1][5],linewidth=2, ls = :dot, alpha = 1, label = "n = dolard(sample_size_vect[1])", 
        lc="black", ylims = ylims_, title = title_)
        plot!((1:sample_size_vect[2])/sample_size_vect[2],mean_vect[2][5],linewidth=2, ls = :dashdotdot, alpha = 1, label = "n = dolard(sample_size_vect[2])",lc = "black")
        plot!((1:sample_size_vect[3])/sample_size_vect[3],mean_vect[3][5],linewidth=2, ls = :dashdot,  alpha = 1, label = "n = dolard(sample_size_vect[3])",lc = "black")
        plot!((1:sample_size_vect[4])/sample_size_vect[4],mean_vect[4][5],linewidth=2, ls = :dash, alpha = 1, label = "n = dolard(sample_size_vect[4])",lc = "black")
        plot!([0,1],[true_signatures[5],true_signatures[5]],linewidth=3,  alpha = 1, label = "true signature",lc = "black")
    else
        plot((1:sample_size_vect[1])/sample_size_vect[1],mean_vect[1][5],linewidth=2, ls = :dot, alpha = 1, ylims = ylims_, title = title_, label="",lc = "black")
        plot!((1:sample_size_vect[2])/sample_size_vect[2],mean_vect[2][5],linewidth=2, ls = :dashdotdot, alpha = 1, label="",lc = "black")
        plot!((1:sample_size_vect[3])/sample_size_vect[3],mean_vect[3][5],linewidth=2, ls = :dashdot, alpha = 1, label="",lc = "black")
        plot!((1:sample_size_vect[4])/sample_size_vect[4],mean_vect[4][5],linewidth=2, ls = :dash, alpha = 1, label="",lc = "black")
        plot!([0,1],[true_signatures[5],true_signatures[5]],linewidth=3, alpha = 1, label="",lc = "black")
    end

    for i in 2:5
        i = 5 - i + 1
        h = h_vect[i]
        plot!((1:sample_size_vect[1])/sample_size_vect[1],mean_vect[1][i],linewidth=2, ls = :dot, alpha = 0.2*i+0, label = "",lc = "black")
        plot!((1:sample_size_vect[2])/sample_size_vect[2],mean_vect[2][i],linewidth=2, ls = :dashdotdot,  alpha = 0.2*i+0, label = "", lc = "black")
        plot!((1:sample_size_vect[3])/sample_size_vect[3],mean_vect[3][i],linewidth=2, ls = :dashdot, alpha = 0.2*i+0, label = "", lc = "black")
        plot!((1:sample_size_vect[4])/sample_size_vect[4],mean_vect[4][i],linewidth=2, ls = :dash,  alpha = 0.2*i+0, label = "", lc = "black")
        plot!([0,1],[true_signatures[i],true_signatures[i]],linewidth=3,  label = "",alpha = 0.2*i+0, lc = "black")    
    end

    if legend_
        plot!([0,0],[0,0],linewidth=0, label = " ",alpha = 0, lc = "white")
        h = h_vect[5]
        plot!([0,0],100 .+[true_signatures[5],true_signatures[5]],linewidth=2, label = "h = 0.5",alpha = 0.2*5+0,lc = "black")
        h = h_vect[4]
        plot!([0,0],100 .+[true_signatures[4],true_signatures[4]],linewidth=2, label = "h = 0.4",alpha = 0.2*4+0,lc = "black")
        h = h_vect[3]
        plot!([0,0],100 .+[true_signatures[3],true_signatures[3]],linewidth=2, label = "h = 0.3",alpha = 0.2*3+0,lc = "black")
        h = h_vect[2]
        plot!([0,0],100 .+[true_signatures[2],true_signatures[2]],linewidth=2, label = "h = 0.2",alpha = 0.2*2+0,lc = "black")
        h = h_vect[1]
        plot!([0,0],100 .+[true_signatures[1],true_signatures[1]],linewidth=2, label = "h = 0.1",alpha = 0.2*1+0,lc = "black")
        plot!(legend=:topleft, legendcolumns=1)
    end
    plot!([0,1],[0,0],linewidth=0, label = "", lc = "white")
end

# black and white
function plot_mean_vect_samples_(mean_vect,samples_signatures,sample_size_vect,h_vect,nb_sample_signatures,index_sample_size,title_, ylims_ = [0,2], ls_ = :dot)
    i = index_sample_size
    
    plot((1:sample_size_vect[i])/sample_size_vect[i],mean_vect[i][5],linewidth=3, ls = ls_, alpha = 1, label = "barycenter, h = 0.5",title = title_, ylims = ylims_, lc= "black")
    plot!((1:sample_size_vect[i])/sample_size_vect[i],mean_vect[i][4],linewidth=3, ls = ls_, alpha = 0.8, label = "barycenter, h = 0.4",title = title_, ylims = ylims_, lc= "black")
    plot!((1:sample_size_vect[i])/sample_size_vect[i],mean_vect[i][3],linewidth=3, ls = ls_, alpha = 0.6, label = "barycenter, h = 0.3",title = title_, ylims = ylims_, lc= "black")
    plot!((1:sample_size_vect[i])/sample_size_vect[i],mean_vect[i][2],linewidth=3, ls = ls_, alpha = 0.4, label = "barycenter, h = 0.2",title = title_, ylims = ylims_,lc= "black")
    plot!((1:sample_size_vect[i])/sample_size_vect[i],mean_vect[i][1],linewidth=3, ls = ls_, alpha = 0.2, label = "barycenter, h = 0.1",title = title_, ylims = ylims_, lc= "black")
    
    for j in 1:size(h_vect)[1]
        h = h_vect[j]
        for iter in 1:nb_sample_signatures
            plot!((1:sample_size_vect[i])/sample_size_vect[i],samples_signatures[i][j][:,iter],linewidth=1, alpha = 0.1, label = "", lc= "black")
        end
    end
    plot!([0,0],[0,1],linewidth=0, label = "", lc= "white")
    plot!(legend=:topleft)
end

function compute_mean_vect(random_uniform_function,sq_distance_matrix_function,sample_size_vect,h_vect,mc_size)
    mean_vect = Vector{Vector{Matrix{Float64}}}(undef, size(sample_size_vect)[1])
    samples_signatures = Vector{Vector{Matrix{Float64}}}(undef, size(sample_size_vect)[1])
    nthreads = Threads.nthreads()
    n_per_threads = Int.(floor.(ones(nthreads).*mc_size/nthreads))
    it0 = 1
    while sum(n_per_threads)<mc_size
        n_per_threads[it0] +=1
        it0 +=1
    end
    result =  Vector{Vector{Matrix{Float64}}}(undef, nthreads)
    for i in 1:(size(sample_size_vect)[1])
        sample_size = sample_size_vect[i]
        Nb_neig = [Int(ceil(h*sample_size)) for h in h_vect]
        Threads.@threads for it in 1:nthreads
            result[it] = [dtm_signature(sq_distance_matrix_function(random_uniform_function(sample_size)),Nb_neig) for i in 1:n_per_threads[it]]
        end
        dtm_sigs = vcat(result...)
        v   = [hcat([sort(dtm_sig[j,:]) for dtm_sig in dtm_sigs]...) for j in 1:size(Nb_neig)[1]] # All signatures are sorted
        mean_dtm_sig = [mean(v[j],dims = 2) for j in 1:size(Nb_neig)[1]]
        mean_vect[i] = mean_dtm_sig
        samples_signatures[i] = v
    end
    return [mean_vect, samples_signatures]
end
"""