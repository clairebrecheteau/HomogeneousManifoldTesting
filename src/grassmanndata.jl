#--------------------------------------------------------------------------
# Struct GrassmannData
#--------------------------------------------------------------------------

Grassmanienne23 = Grassmann(3,2)
Grassmanienne24 = Grassmann(4,2)

export GrassmannData23
export GrassmannData24

"""
$(TYPEDEF)

GrassmannData to store data on a circle
- np: number of points
- points: vector containing the coordinates of the np points (size = np)
"""
struct GrassmannData23{T} <: AbstractData

    np::Int
    points::Vector{T}

end

struct GrassmannData24{T} <: AbstractData

    np::Int
    points::Vector{T}

end


#--------------------------------------------------------------------------
# Constructor for GrassmannData23 et GrassmannData24
#--------------------------------------------------------------------------

"""
$(SIGNATURES)

- np: sample size

Returns a sample of size np, uniformly distributed on the circle.
Values are in the interval [0,2π].
"""

function GrassmannData23(
    np::Integer,
    method::String,
    params...
)

    @assert method in ["uniform","mixture_normal"] "Parameter method should be in ['uniform','mixture_normal']"
    if method == "uniform"
        points = random_uniform_Grassmann23_data(np)
    elseif method == "mixture_normal"
        points = random_mixture_Grassmann23_data(np,params...)
    end

    return GrassmannData23(size(points)[1],points)

end

function GrassmannData24(
    np::Integer,
    method::String,
    params...
)

    @assert method in ["uniform", "mixture_normal"] "Parameter method should be in ['uniform','mixture_normal']"
    if method == "uniform"
        points = random_uniform_Grassmann24_data(np)
    elseif method == "mixture_normal"
        points = random_mixture_Grassmann24_data(np,params...)
    end

    return GrassmannData24(size(points)[1],points)

end

function random_uniform_Grassmann23_data(
    np::Integer
)

    points = rand(Grassmanienne23,np)

    return points
end

function random_mixture_Grassmann23_data(
    np::Integer,
    sig = 1.0,
    centers_ = [hcat([1, 0, 0] .* 1.0,[0, 1, 0] .* 1.0),hcat([1, 0, 0] .* 1.0,[0, 0, 1] .* 1.0),hcat([0, 1, 0] .* 1.0,[0, 0, 1] .* 1.0)]
)
    M = Grassmann(3,2)
    normal_distrib = Normal(0,1)
    points = Vector{Matrix{Float64}}(undef,np)
    for i in 1:np
        G_ = centers_[rand(1:length(centers_),1)][1]
        σ = 1.0
        X = sig*rand(normal_distrib, 1)[1]*rand(M; σ, vector_at = G_)
        Y = exp(M, G_, X)
        points[i] = Y
    end

    return points
end

function random_uniform_Grassmann24_data(
    np::Integer
)

    points = rand(Grassmanienne24,np)

    return points
end

function random_mixture_Grassmann24_data(
    np::Integer,
    sig = 1.0,
    centers_ = [hcat([1, 0, 0, 0] .* 1.0,[0, 1, 0, 0] .* 1.0),hcat([1, 0, 0, 0] .* 1.0,[0, 0, 1, 0] .* 1.0),hcat([1, 0, 0, 0] .* 1.0,[0, 0, 0, 1] .* 1.0),hcat([0, 1, 0, 0] .* 1.0,[0, 0, 1, 0] .* 1.0),hcat([0, 1, 0, 0] .* 1.0,[0, 0, 0, 1] .* 1.0),hcat([0, 0, 1, 0] .* 1.0,[0, 0, 0, 1] .* 1.0)]
)
    M = Grassmann(4,2)
    normal_distrib = Normal(0,1)
    points = Vector{Matrix{Float64}}(undef,np)
    for i in 1:np
        G_ = centers_[rand(1:length(centers_),1)][1]
        σ = 1.0
        X = sig*rand(normal_distrib, 1)[1]*rand(M; σ, vector_at = G_)
        Y = exp(M, G_, X)
        points[i] = Y
    end

    return points
end

#--------------------------------------------------------------------------
# Function to Plot data
#--------------------------------------------------------------------------

"""
$(SIGNATURES)

- data_grassmann_23: data on the torus
- alpha : transparency parameter for points

Plot the data on the torus.
"""

function plot_data(
    data_grassmann_23::GrassmannData23,
    alpha = 0.2
)
# We plot the orthogonal direction x onto the sphere (2 points : x and -x for each point in the Grassmannian)
    N = 32
    u = range(0, stop=2π, length=N)
    v = range(0, stop=π, length=N)
    x = cos.(u) .* sin.(v)'
    y = sin.(u) .* sin.(v)'
    z = repeat(cos.(v)',outer=[N, 1])
    plot(x,y,z,color = "black",legend = false,grid = false, axis = false,alpha = 0.3)
    ortho_vectors = [[x[2,1]*x[3,2] - x[3,1]*x[2,2] , x[3,1]*x[1,2] - x[1,1]*x[3,2] , x[1,1]*x[2,2] - x[1,2]*x[2,1]] for x in data_grassmann_23.points]
    vectors = hcat(hcat([x/norm(x) for x in ortho_vectors],[-x/norm(x) for x in ortho_vectors])...)
    plot!(vectors[1,:],vectors[2,:],vectors[3,:],seriestype=:scatter,legend = false, alpha = alpha,showaxis = false, color = 2)
end

function plot_data(
    data_grassmann_24::GrassmannData24,
    dims = [1,2,3],
    alpha = 0.2
)

    @assert dims == [1,2,3] || dims == [1,2,4]  || dims == [1,3,4]  || dims == [2,3,4] "dims should be [1,2,3], [1,2,4], [1,3,4] or [2,3,4]."
    data_grassmann_23 = GrassmannData23(data_grassmann_24.np,[ x[dims,:] for x in data_grassmann_24.points])
    plot_data(data_grassmann_23,alpha)
end





#----------------------------------------------------------
# True Signature
#----------------------------------------------------------

export squared_distance_matrix

function squared_distance_matrix(
    data::GrassmannData23
)

    return [[distance(Grassmanienne23,x,y) for x in data.points] for y in data.points]

end

function squared_distance_matrix(
    data::GrassmannData24
)

    return [[distance(Grassmanienne24,x,y) for x in data.points] for y in data.points]

end

function squared_distance_matrix(
    data1::GrassmannData23,
    data2::GrassmannData23
)
    return [[distance(Grassmanienne23,x,y) for x in data1.points] for y in data2.points]

end

function squared_distance_matrix(
    data1::GrassmannData24,
    data2::GrassmannData24
)
    return [[distance(Grassmanienne24,x,y) for x in data1.points] for y in data2.points]

end
