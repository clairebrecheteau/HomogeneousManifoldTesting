#--------------------------------------------------------------------------
# Struct TorusData
#--------------------------------------------------------------------------

export TorusData

"""
$(TYPEDEF)

TorusData to store data on the torus
- np: number of points
- points: vector containing the coordinates of the np points (size = (2,np))
"""
struct TorusData{T} <: AbstractData

    np::Int
    points::Matrix{T}

end

#--------------------------------------------------------------------------
# Constructor for TorusData
#--------------------------------------------------------------------------

"""
$(SIGNATURES)

- np: sample size

Returns a sample of size np, on the flat torus.
"""

function TorusData(
    np::Integer,
    method::String,
    params...
)

    @assert method in ["uniform","normal","mixture normal", "large circle incremental", "large circle uniform", "large circle equidistant", "large ellipse incremental", "large ellipse uniform", "large circle line incremental", "large circle line uniform", "dynamical system", "dynamical system2", "dynamical system3", "Faure", "grid"] "Parameter method should be in ['uniform','normal','mixture normal', 'large circle incremental', 'large circle uniform',  'large circle equidistant', 'large ellipse incremental', 'large ellipse uniform', 'large circle line incremental', 'large circle line uniform', 'dynamical system1','dynamical system2', 'dynamical system3', 'Faure', 'grid']"
    if method == "uniform"
        points = random_uniform_torus_data(np)
    elseif method == "normal"
        points = random_normal_torus_data(np, params...)
    elseif method == "mixture normal"
        points = random_mixture_normal_torus_data(np,params...)
    elseif method == "large circle incremental"
        points = random_large_circle_incremental_angle_torus_data(np,params...)
    elseif method == "large circle uniform"
        points = random_large_circle_uniform_angle_torus_data(np, params...)
    elseif method == "large circle equidistant"
        points = random_large_circle_equidistant_angle_torus_data(np,params...)
    elseif method == "large ellipse incremental"
        points = random_large_ellipse_incremental_angle_torus_data(np,params...)
    elseif method == "large ellipse uniform"
        points = random_large_ellipse_uniform_angle_torus_data(np,params...)
    elseif method == "large circle line incremental"
        points = random_large_circle_line_incremental_angle_torus_data(np,params...)
    elseif method == "large circle line uniform"
        points = random_large_circle_line_uniform_angle_torus_data(np,params...)
    elseif method == "dynamical system"
        points = random_dynamical_system_torus_data(np,params...)
    elseif method == "dynamical system2"
        points = random_dynamical_system_two_torus_data(np,params...)
    elseif method == "dynamical system3"
        points = random_dynamical_system_three_torus_data(np,params...)
    elseif method == "Faure"
        points = random_Faure_torus(np,params...)
    elseif method == "grid"
        points = random_grid_normal_torus_data(np, params...)
    end

    return TorusData(size(points)[2],points)

end

function random_uniform_torus_data(
    np::Integer
)

    points = collect(transpose(hcat(rand(np),rand(np))))

    return points
end 

function random_normal_torus_data(
    np::Integer,
    μ::Vector{<:Real},
    Σ::Matrix{<:Real}
)

    mv_normal = MvNormal(μ,Σ)
    normal_points = rand(mv_normal,np)
    points = (normal_points.%1 .+1) .%1

    return points
end

function random_mixture_normal_torus_data(
    np::Integer,
    vect_μ::Vector{<:Vector{<:Real}},
    vect_Σ::Vector{<:Matrix{<:Real}},
    proba::Vector{<:Real}
)

    @assert size(vect_μ)==size(vect_Σ) "The number of centers and the number of standard deviations parameter should be equal."
    @assert size(vect_μ)==size(proba) "The number of centers and the number of probabilities should coincide."
    @assert sum(proba.<0)==0 "The probabilities of each component should be non negative."
    @assert sum(proba)==1 "The sum of the probabilitites should be 1."
    n_components = size(vect_μ)[1]
    n_per_comp = rand(Multinomial(np, proba))
    points = hcat([random_normal_torus_data(n_per_comp[i], vect_μ[i], vect_Σ[i]) for i in 1:n_components]...)

    return points
end

function random_large_circle_incremental_angle_torus_data(
    np::Integer,
    r::Real,
    p0::Vector{<:Real}
)

    U0 = rand(1)[1]
    Theta = [2*pi*i*U0 for i in 1:np]
    points = simulate_large_circle_torus(np,r,p0,Theta)

    return points
end

function random_large_circle_uniform_angle_torus_data(
    np::Integer,
    r::Real,
    p0::Vector{<:Real}
)

    Theta = (2*pi).*rand(np)
    points = simulate_large_circle_torus(np,r,p0,Theta)

    return points
end

function random_large_circle_equidistant_angle_torus_data(
    np::Integer,
    r::Real,
    p0::Vector{<:Real}
)

    Theta = Vector((1:np)*2*pi/np)
    points = simulate_large_circle_torus(np,r,p0,Theta)

    return points
end

function simulate_large_circle_torus(
    np::Integer,
    r::Real,
    p0::Vector{<:Real},
    Theta::Vector{<:Real}
)

    points = hcat([[((p0[1] + r*cos(Theta[i])) - Int(floor(p0[1] + r*cos(Theta[i])))) % 1, (p0[2] + r*sin(Theta[i]) - Int(floor(p0[2] + r*sin(Theta[i])))) % 1] for i in 1:np]...)

    return points
end

function random_large_ellipse_incremental_angle_torus_data(
    np::Integer,
    p0::Vector{<:Real},
    r1::Real,
    r2::Real,
    angle::Real
)

    U0 = rand(1)[1]
    Theta = [2*pi*i*U0 for i in 1:np]
    points = random_large_ellipse_torus_data(np,p0,r1,r2,angle,Theta)

    return points

end

function random_large_ellipse_uniform_angle_torus_data(
    np::Integer,
    p0::Vector{<:Real},
    r1::Real,
    r2::Real,
    angle::Real
)

    Theta = (2*pi).*rand(np)
    points = random_large_ellipse_torus_data(np,p0,r1,r2,angle,Theta)

    return points

end

function random_large_circle_line_incremental_angle_torus_data(
    np::Integer,
    p0::Vector{<:Real},
    r::Real,
    angle::Real,
    angle1::Real
)

    U0 = rand(1)[1]
    Theta = (([2*pi*i*U0 for i in 1:np]).%(2*pi) .+(2*pi)).%(2*pi)
    points = random_large_circle_line_torus_data(np,p0,r,angle,angle1,Theta)

    return points

end

function random_large_circle_line_uniform_angle_torus_data(
    np::Integer,
    p0::Vector{<:Real},
    r::Real,
    angle::Real,
    angle1::Real
)

    Theta = (2*pi).*rand(np)
    points = random_large_circle_line_torus_data(np,p0,r,angle,angle1,Theta)

    return points

end

function random_large_ellipse_torus_data(
    np::Integer,
    p0::Vector{<:Real},
    r1::Real,
    r2::Real,
    angle::Real,
    Theta::Vector{<:Real}
)

    # angle indicates the rotation of the ellipse.
    u = [cos(angle), sin(angle)]
    v = [cos(angle + pi/2), sin(angle + pi/2)]
    points = hcat([[((p0[1] + r1*cos(Theta[i])*u[1] + r2*sin(Theta[i])*v[1]) - Int(floor(p0[1] + r1*cos(Theta[i])*u[1] + r2*sin(Theta[i])*v[1]))) % 1, (p0[2] + r1*cos(Theta[i])*u[2] + r2*sin(Theta[i])*v[2] - Int(floor(p0[2] + r1*cos(Theta[i])*u[2] + r2*sin(Theta[i])*v[2]))) % 1] for i in 1:np]...)

    return points
end

function random_large_circle_line_torus_data(
    np::Integer,
    p0::Vector{<:Real},
    r::Real,
    angle::Real,
    angle1::Real,
    Theta::Vector{<:Real}
)

    # angle indicates the rotation of the shape.
    # angle1 indicates the angle for which the circle is transformed into a line.
    u = [cos(angle), sin(angle)]
    v = [cos(angle + pi/2), sin(angle + pi/2)]
    angle1 = max(angle1 % pi,2*pi/np) # To have a segment at most on half the circle.
    angles_segment = Theta[(Theta .<= angle1) .|| (Theta .>= 2*pi - angle1)]
    angles_circle = Theta[(Theta .> angle1) .&& (Theta .< 2*pi - angle1)]
    points_segment = hcat([[((p0[1] + r*cos(angle1)*u[1] + r*cos(angle1)*sin(ang)/cos(ang)*v[1]) - Int(floor(p0[1] + r*cos(angle1)*u[1] + r*cos(angle1)*sin(ang)/cos(ang)*v[1]))) % 1, (p0[2] + r*cos(angle1)*u[2] + r*cos(angle1)*sin(ang)/cos(ang)*v[2] - Int(floor(p0[2] + r*cos(angle1)*u[2] + r*cos(angle1)*sin(ang)/cos(ang)*v[2]))) % 1] for ang in angles_segment]...)
    points_circle = hcat([[((p0[1] + r*cos(ang)*u[1] + r*sin(ang)*v[1]) - Int(floor(p0[1] + r*cos(ang)*u[1] + r*sin(ang)*v[1]))) % 1, (p0[2] + r*cos(ang)*u[2] + r*sin(ang)*v[2] - Int(floor(p0[2] + r*cos(ang)*u[2] + r*sin(ang)*v[2]))) % 1] for ang in angles_circle]...)
    points = hcat(points_segment,points_circle)
    
    return points
end



function random_dynamical_system_torus_data(
    np::Integer,
    alp::Real,
    bet::Real,
    Nstart::Integer = 0
)

    points = hcat([[n*alp % 1, n*bet % 1] for n in Nstart:(np+Nstart-1)]...)

    return points
end

function random_dynamical_system_two_torus_data(
    np::Integer,
    p0::Vector{<:Real}
)

    A = [2 1; 1 1]
    points = Matrix{typeof(p0[1])}(undef, (2,np))
    points[:,1] = p0
    for i in 2:np
        points[:,i] = (A*points[:,i-1]).% 1
    end

    return points
end

function random_dynamical_system_three_torus_data(
    np::Integer,
    p0::Vector{<:Real}
)

    points = Matrix{typeof(p0[1])}(undef, (2,np))
    points[:,1] = p0
    for i in 2:np
        points[:,i] = (points[:,i-1].*points[:,i-1]).% 1
    end

    return points
end

function random_Faure_torus(
    np::Integer,
    p0::Vector{<:Real},
    #p1::Vector{<:Real},
    r::Integer
)
 # On part des points p0, p0 + p1, ... p0 + (np-1)p1
 # On leur applique la r fois la transformation cat map (2 1 ; 1 1)

    points = Matrix{typeof(p0[1])}(undef, (2,np))

    points[:,1] = (p0).% 1
    for i in 2:np
        points[:,i] = (p0).% 1 #(p0 .+ (i .* p1) ).% 1
    end

    transfo = [3 2 ; 1 1]

    catmap = [2 1; 1 1]

    #for j in 1:r
    #    for i in 1:np
    #        points[:,i] = (catmap * points[:,i]) .% 1
    #    end
    #end

    for i in 1:np
        for j in 1:i
            points[:,i] = ((transfo) * points[:,i]) .% 1 #((transfo^i) * points[:,i]) .% 1
        end
    end

    for i in 1:np
        for j in 1:r
            points[:,i] = ((catmap) * points[:,i]) .% 1
        end
    end

    #for i in 1:np
    #    points[:,i] = ((catmap^r) * points[:,i]) .% 1
    #end

    return points
end

function random_grid_normal_torus_data(
    np::Integer,
    σ::Real,
    grid_size::Integer
)
    normal_distrib = Normal(0,1)

    dat = (((rand(1:grid_size,np*2) ./ grid_size .+ σ .* rand(normal_distrib,np*2)) .% 1) .+ 1) .%1;
    points = collect(transpose(hcat(dat[1:np],dat[(np + 1):(2*np)])))

    return points
end

#--------------------------------------------------------------------------
# Function to Plot data
#--------------------------------------------------------------------------

"""
$(SIGNATURES)

- data_torus: data on the torus
- alpha : transparency parameter for points

Plot the data on the torus.
"""

function plot_data(
    data_torus::TorusData,
    alpha = 0.2
)

    plot([0, 1, 1, 0, 0], [0, 0, 1, 1, 0],color = "black",aspect_ratio=:equal, legend = false, xlims = [-0.6,1.7], ylims = [-0.2,1.2],showaxis = false,grid = false)
    plot!(data_torus.points[1,:],data_torus.points[2,:],seriestype=:scatter, alpha = alpha, color = 2)

end


#----------------------------------------------------------
# True Signature
#----------------------------------------------------------

function squared_distance_matrix(
    data::TorusData
)

    @assert sum(data.points.<0)==0 "The points should have non-negative coordinates."
    @assert sum(data.points.>1)==0 "The points should have coordinates not superior to 1."

    return [sq_distance_flat_torus(view(data.points,:,i),data.points) for i in 1:data.np]

end

function sq_distance_flat_torus(
    x::SubArray{T, 1, Matrix{T}, Tuple{Base.Slice{Base.OneTo{S}}, S}, true},#Vector{<:Real},
    vect_y::Matrix{T}
) where {T<:Real} where {S<:Real}

    # The distance between x and y is the min between |x-y+1|, |x-y-1| and |y-x|
    return vec(sum(min.(mod.(abs.(vect_y .- x),1),1 .- mod.(abs.(vect_y .- x),1)).^2, dims = 1))
end

function sq_distance_flat_torus(x::Vector{Float64},y::Vector{Float64})
    # The distance between x and y is the min between |x-y+1|, |x-y-1| and |y-x|
    return sum(min.(mod.(abs.(x-y),1),1 .- mod.(abs.(x-y),1)).^2)
end


"""
function squared_distance_matrix(
    data::TorusData
)

    @assert sum(data.points.<0)==0 "The points should have non-negative coordinates."
    @assert sum(data.points.>1)==0 "The points should have coordinates not superior to 1."

    return [[sq_distance_flat_torus(data.points[:,i],data.points[:,j]) for i in 1:data.np] for j in 1:data.np]

end
"""