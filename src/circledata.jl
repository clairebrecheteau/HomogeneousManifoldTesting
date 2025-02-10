#--------------------------------------------------------------------------
# Struct CircleData
#--------------------------------------------------------------------------

export CircleData

"""
$(TYPEDEF)

CircleData to store data on a circle
- np: number of points
- points: vector containing the coordinates of the np points (size = np)
"""
struct CircleData{T} <: AbstractData

    np::Int
    points::Vector{T}

end

#--------------------------------------------------------------------------
# Constructor for CircleData
#--------------------------------------------------------------------------

"""
$(SIGNATURES)

- np: sample size

Returns a sample of size np, uniformly distributed on the circle.
Values are in the interval [0,2π].
"""
    

function CircleData(
    np::Integer,
    method::String,
    params...
)

    @assert method in ["uniform","normal","mixture normal", "dynamical system", "dynamical system 2", "FvM", "MFvM"] "Parameter method should be in ['uniform','normal','mixture normal','dynamical system', 'dynamical system 2', 'FvM', 'MFvM']"
    if method == "uniform"
        points = random_uniform_circle_data(np)
    elseif method == "FvM"
        points = random_von_Mises_Fischer_circle_data(np, params...)
    elseif method == "MFvM"
        points = random_mixture_von_Mises_Fischer_circle_data(np,params...)
    elseif method == "normal"
        points = random_normal_circle_data(np, params...)
    elseif method == "mixture normal"
        points = random_mixture_normal_circle_data(np,params...)
    elseif method == "dynamical system"
        points = random_dynamical_system_circle_data(np,params...)
    elseif method == "dynamical system 2"
        points = random_dynamical_system_2_circle_data(np,params...)
    end

    return CircleData(size(points)[1],points)

end


function random_von_Mises_Fischer_circle_data(
    np::Integer,
    κ::Real,
    center_ = [1.0,0]
)

    R"""
    library(sphunif);
    """;

    @rput(np)
    kappa = κ
    @rput(kappa)
    @rput(center_)

    R"""
    sample = r_alt(n = np, p = 2, alt = "vMF",kappa = kappa, mu = center_)[, , 1]
    """
    points = collect(transpose(@rget(sample)))
    angles = acos.(points[1,:]) + 2*(pi .- acos.(points[1,:])).* (points[2,:] .<0)

    return angles
end


function random_mixture_von_Mises_Fischer_circle_data(
    np::Integer,
    κ::Real
)

    R"""
    library(sphunif);
    #library(foreach);
    #library(doParallel);
    #nCores <- 10;
    #registerDoParallel(nCores);
    """;

    @rput(np)
    kappa = κ
    @rput(kappa)
    @rput(center)

    R"""
    sample = r_alt(n = np, p = 2, alt = "MvMF",kappa = kappa)[, , 1]
    """
    points = collect(transpose(@rget(sample)))
    angles = acos.(points[1,:]) + 2*(pi .- acos.(points[1,:])).* (points[2,:] .<0)

    return angles
end


# A FAIRE (passer en dim 1)


function random_uniform_circle_data(
    np::Integer
)

    points = rand(np)*2*π

    return points
end


function random_dynamical_system_circle_data(
    np::Integer,
    increment::Real,
    index_start=0::Real
)

    points = [(n*increment) % (2*π) for n in index_start:(np+index_start-1)]

    return points
end


function random_normal_circle_data(
    np::Integer,
    mean_::Real,
    std_::Real
)

    norm_dist = Normal(mean_, std_)
    norm_sample = rand(norm_dist,np)
    points = norm_sample.%(2*pi) .+ (norm_sample.<0).*(2*pi)

    return points

end


function random_mixture_normal_circle_data(
    np::Integer,
    mean_::Vector{<:Real},
    std_::Vector{<:Real},
    proba::Vector{<:Real}
)

    @assert size(mean_)==size(std_) "The number of centers and the number of standard deviations parameter should be equal."
    @assert size(mean_)==size(proba) "The number of centers and the number of probabilities should coincide."
    @assert sum(proba.<0)==0 "The probabilities of each component should be non negative."
    @assert sum(proba)==1 "The sum of the probabilitites should be 1."
    n_components = size(mean_)[1]
    n_per_comp = rand(Multinomial(np, proba))
    points = vcat([random_normal_circle_data(n_per_comp[i], mean_[i], std_[i]) for i in 1:n_components]...)
    return points

end

function random_dynamical_system_2_circle_data(
    np::Integer,
    r::Real
)

    U = 2*π*rand(Float64,1)[1]
    points = [(exp(r*log(2))*(n + U)) % (2*π) for n in 1:np]

    return points
end

function random_dynamical_system_Faure_1_circle_data(
    np::Integer,
    r::Real
)

    U = 2*π*rand(Float64,1)[1]
    points = [(n + U) % (2*π) for n in 1:np]
    for i in 1:ceil(r)
        points2 = copy(points)
        points = (2 .* points2) .% (2*π)
    end

    return points
end


#--------------------------------------------------------------------------
# Function to Plot data
#--------------------------------------------------------------------------

"""
$(SIGNATURES)

- data_circle: data on the circle
- alpha : transparency parameter for points

Plot the data on the circle.
"""

function plot_data(
    data_circle::CircleData,
    alpha = 0.2
)

    x = [2*π/10000*k for k in 1:10000]
    plot(cos.(x),sin.(x),color = "black",aspect_ratio=:equal, legend = false, xlims = [-1.7,1.7], ylims = [-1.2,1.2],showaxis = false, grid = false)
    plot!(cos.(data_circle.points),sin.(data_circle.points),seriestype=:scatter, alpha = alpha, color = 2)

end


#----------------------------------------------------------
# True Signature
#----------------------------------------------------------

export squared_distance_matrix

function squared_distance_matrix(
    data::CircleData
)

    @assert sum(data.points.<0)==0 "The points should have non-negative coordinates."
    @assert sum(data.points.>2*π)==0 "The points should have coordinates not superior to 2*π."

    return [[sq_distance_circle(x,y) for x in data.points] for y in data.points]

end

function squared_distance_matrix(
    data1::CircleData,
    data2::CircleData
)

    @assert sum(data1.points.<0)==0 "The points should have non-negative coordinates."
    @assert sum(data1.points.>2*π)==0 "The points should have coordinates not superior to 2*π."
    @assert sum(data2.points.<0)==0 "The points should have non-negative coordinates."
    @assert sum(data2.points.>2*π)==0 "The points should have coordinates not superior to 2*π."

    return [[sq_distance_circle(x,y) for x in data1.points] for y in data2.points]

end

function sq_distance_circle(
    x::Float64,
    y::Float64
)
    # The distance between x and y is the min between |x-y-2\pi radius| and |y-x|
    #return min(abs(x-y),abs(x-y-2*π*radius),abs(x-y+2*π*radius))
    return min(mod.(abs.(x-y),2*π),2*π .- mod.(abs.(x-y),2*π))^2 #???

end