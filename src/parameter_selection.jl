export Bonferroni_procedure

function Bonferroni_procedure(
    proportion_nearest_neighbours::Vector{<:Real},
    pvalues::Vector{<:Vector{<:Real}},
    alpha::Real
)
    I = size(proportion_nearest_neighbours)[1]
    n_rep = size(pvalues)[1]

    @assert (size(pvalues[1])[1] == I) "There should be the same number of pvalues than nearest neighbours parameters."
    @assert (alpha > 0) && (alpha < 1) "Parameter alpha should be between 0 and 1."
    
    threshold = alpha/I

    rejection_indices = [findall(pval .<= threshold) for pval in pvalues]
    selected_parameters = [proportion_nearest_neighbours[ind] for ind in rejection_indices]

    return [rejection_indices, selected_parameters]
end


export Benjamini_Hochberg_procedure

function Benjamini_Hochberg_procedure(
    proportion_nearest_neighbours::Vector{<:Real},
    pvalues::Vector{<:Vector{<:Real}},
    alpha::Real
)

    I = size(proportion_nearest_neighbours)[1]
    n_rep = size(pvalues)[1]

    @assert (size(pvalues[1])[1] == I) "There should be the same number of pvalues than nearest neighbours parameters."
    @assert (alpha > 0) && (alpha < 1) "Parameter alpha should be between 0 and 1."

    thresholds = [i*alpha/I for i in 1:I]

    indices_pval = [sortperm(pval) for pval in pvalues]
    sorted_pvalues = [pvalues[index][indices_p] for (index,indices_p) in enumerate(indices_pval)]
    indices_max = [findlast(pval .<= thresholds) for pval in sorted_pvalues]

    rejection_indices = Vector{Vector{Int64}}(undef,n_rep)

    type_ = eltype(indices_pval[1])

    for n in 1:n_rep

        if indices_max[n]==nothing
            # pas de rejet
            rej_ind = zeros(type_,0)
        else
            ind_max = indices_max[n]
            ind_0 = 1
            while indices_pval[n][ind_0] !=  ind_max
                ind_0 += 1
            end
            rej_ind = indices_pval[n][1:ind_0]
        end

        rejection_indices[n] = rej_ind

    end

    selected_parameters = [proportion_nearest_neighbours[ind] for ind in rejection_indices]

    return [rejection_indices, selected_parameters]
end