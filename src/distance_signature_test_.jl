function distance_signatures_test_homogeneity(
    signatures::Matrix{<:Real},
    true_signature::Real,
    np::Integer
)

    distance = 1/np*sum(abs.(signatures .- true_signature), dims = 1)

    return  distance[1,:]

end

function distance_signatures_test_homogeneity(
    signatures::Vector{<:Real},
    true_signature::Real,
    np::Integer
)

    return 1/np*sum(abs.(signatures .- true_signature))

end

function distance_signatures_test_iidness(
    signatures::Matrix{<:Real},
    barycenter_signatures::Matrix{<:Real},
    np::Integer
)

    distance = sqrt.(1/np.*sum((signatures .- barycenter_signatures).^2, dims = 1))
    return distance[1,:]

end

function distance_signatures_test_iidness(
    signatures::Vector{<:Real},
    barycenter_signatures::Matrix{<:Real},
    np::Integer
)

    return sqrt(1/np*sum((signatures.- barycenter_signatures).^2))

end