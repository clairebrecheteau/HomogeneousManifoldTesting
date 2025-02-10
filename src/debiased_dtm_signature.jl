export debiased_dtm_signature

function debiased_dtm_signature(
    squared_distances_matrix::Vector{<:Vector{<:Real}},
    proportion_nearest_neighbours::Vector{<:Real}
)

    np = size(squared_distances_matrix)[1]

    indices = 1:np
    indicator_vector = [(indices .<= floor( h*(np-1) + 1)) + ((h*(np-1))%1)* (indices .== ceil( h*(np-1) + 1)) for h in proportion_nearest_neighbours]
    dtm_sig = vcat([hcat([sqrt.(1/(proportion_nearest_neighbours[j]*(np-1))*sum( (sort(squared_distances_matrix[i])).*indicator_vector[j] )) for i in 1:np]...) for j in 1:size(proportion_nearest_neighbours)[1]]...)

    return dtm_sig

end