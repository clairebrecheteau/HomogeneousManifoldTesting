

function compute_true_signatures(
    data_type::UnionAll,
    proportion_nearest_neighbours::Vector{<:Real}
)

    @assert data_type in AbstractDataList "The data_type should belong to [CircleData,SphereData,TorusData,BolzaData]."

    if data_type == CircleData
        true_signatures = get_true_signatures_circle(proportion_nearest_neighbours)
    elseif data_type == SphereData
        true_signatures = get_true_signatures_sphere(proportion_nearest_neighbours)
    elseif data_type == TorusData
        true_signatures = get_true_signatures_flat_torus(proportion_nearest_neighbours)
    elseif data_type == BolzaData
        true_signatures = get_true_signatures_Bolza(proportion_nearest_neighbours)
    elseif data_type == GrassmannData23
        true_signatures = get_true_signatures_Grassmann23(proportion_nearest_neighbours)
    elseif data_type == GrassmannData24
        true_signatures = get_true_signatures_Grassmann24(proportion_nearest_neighbours)
    end


return true_signatures

end


function get_true_signatures_circle(
    proportion_nearest_neighbours::Vector{<:Real}
)

    [1/sqrt(3)*π*h for h in proportion_nearest_neighbours]

end

function get_true_signatures_sphere(
    proportion_nearest_neighbours::Vector{<:Real}
)

    nb_parameters = size(proportion_nearest_neighbours)[1]
    true_signatures = zeros(nb_parameters)
    for i in 1:nb_parameters
        h = proportion_nearest_neighbours[i]
        if h >= 0.5
            rh = acos(1-2*h)
        else
            rh = pi - acos(-1 + 2*h)
        end
        true_signatures[i] = sqrt(1/(2*h)*(-rh^2*cos(rh) + 2*rh*sin(rh) + 2*cos(rh) - 2))
    end

    return true_signatures

end

function get_true_signatures_flat_torus(
    proportion_nearest_neighbours::Vector{<:Real}
)

    [sqrt(h/(2*pi)) for h in proportion_nearest_neighbours]

end

function get_true_signatures_Bolza(
    proportion_nearest_neighbours::Vector{<:Real}
)

    nb_parameters = size(proportion_nearest_neighbours)[1]
    true_signatures = zeros(nb_parameters)
    for i in 1:nb_parameters
        h = proportion_nearest_neighbours[i]
        V = 12.5664
        true_signatures[i] = sqrt(2*((asinh(sqrt((h*V)/(4*pi))))^2 + ((asinh(sqrt((h*V)/(4*pi))))*sqrt((4*pi)/(h*V) + 1) - 1)^2))
    end

    return true_signatures

end

function get_true_signatures_Grassmann23(
    proportion_nearest_neighbours::Vector{<:Real}
)
    nb_parameters = size(proportion_nearest_neighbours)[1]
    true_signatures = zeros(nb_parameters)

    #todo

    return true_signatures
end

function get_true_signatures_Grassmann24(
    proportion_nearest_neighbours::Vector{<:Real}
)
    nb_parameters = size(proportion_nearest_neighbours)[1]
    true_signatures = zeros(nb_parameters)
    
    #todo

    return true_signatures
end