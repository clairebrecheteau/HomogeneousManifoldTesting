#--------------------------------------------------------------------------
# Struct SphereData
#--------------------------------------------------------------------------

export SphereData

"""
$(TYPEDEF)

SphereData to store data on the sphere S^2.
- np: number of points
- points: array containing the coordinates of the np points (size = (3,np))
"""
struct SphereData{T} <: AbstractData

    np::Int
    points::Matrix{T}

end


#--------------------------------------------------------------------------
# Constructor for SphereData
#--------------------------------------------------------------------------

"""
$(SIGNATURES)

- np: sample size

Returns a sample of size np, on the sphere S^2.
"""

function SphereData(
    np::Integer,
    method::String,
    params...
)

    @assert method in ["uniform","FvM","MFvM"] "Parameter method should be in ['uniform', 'FvM', 'MFvM']"
    if method == "uniform"
        points = random_uniform_sphere_data(np)
    elseif method == "FvM"
        points = random_von_Mises_Fischer_sphere_data(np, params...)
    elseif method == "MFvM"
        points = random_mixture_von_Mises_Fischer_sphere_data(np, params...)
    end

    return SphereData(size(points)[2],points)

end


function random_uniform_sphere_data(
    np::Integer
)

    norm_dist = Normal(0, 1)
    norm_sample = rand(norm_dist,(3,np))
    points = norm_sample./sqrt.(sum((norm_sample.^2),dims = 1))

    return points

end


function random_von_Mises_Fischer_sphere_data(
    np::Integer,
    κ::Real,
    center = [1.0,0,0]
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
    sample = r_alt(n = np, p = 3, alt = "vMF",kappa = kappa, mu = center)[, , 1]
    """
    points = collect(transpose(@rget(sample)))

    return points
end


function random_mixture_von_Mises_Fischer_sphere_data(
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
    sample = r_alt(n = np, p = 3, alt = "MvMF",kappa = kappa)[, , 1]
    """
    points = collect(transpose(@rget(sample)))

    return points
end


#--------------------------------------------------------------------------
# Function to Plot data
#--------------------------------------------------------------------------

"""
$(SIGNATURES)

- data_sphere: data on the sphere
- alpha : transparency parameter for points

Plot the data on the sphere.
"""

function plot_data(
    data_sphere::SphereData,
    alpha = 0.2
)
    N = 32
    u = range(0, stop=2π, length=N)
    v = range(0, stop=π, length=N)
    x = cos.(u) .* sin.(v)'
    y = sin.(u) .* sin.(v)'
    z = repeat(cos.(v)',outer=[N, 1])
    plot(x,y,z,color = "black",legend = false,grid = false, axis = false,alpha = 0.3)
    plot!(data_sphere.points[1,:],data_sphere.points[2,:],data_sphere.points[3,:],seriestype=:scatter,legend = false, alpha = alpha,showaxis = false, color = 2)

end

export plot_data_clustering

function plot_data_clustering(
    data_sphere::SphereData,
    centers::Vector{<:Vector{<:Real}},
    colors::Vector{<:Integer},
    alpha = 0.4
)
    N = 32
    u = range(0, stop=2π, length=N)
    v = range(0, stop=π, length=N)
    x = cos.(u) .* sin.(v)'
    y = sin.(u) .* sin.(v)'
    z = repeat(cos.(v)',outer=[N, 1])
    plot(x,y,z,color = "black",legend = false,grid = false, axis = false,alpha = 0.3, ms = 4)
    points_x = [data_sphere.points[1,i] for i in 1:data_sphere.np if colors[i]!= 0]
    points_y = [data_sphere.points[2,i] for i in 1:data_sphere.np if colors[i]!= 0]
    points_z = [data_sphere.points[3,i] for i in 1:data_sphere.np if colors[i]!= 0]
    colors_ = [colors[i] for i in 1:data_sphere.np if colors[i]!= 0]
    plot!(points_x,points_y,points_z,seriestype=:scatter,legend = false ,showaxis = false, color = colors_, alpha = alpha, ms= 4)
    noise_x = [data_sphere.points[1,i] for i in 1:data_sphere.np if colors[i]== 0]
    noise_y = [data_sphere.points[2,i] for i in 1:data_sphere.np if colors[i]== 0]
    noise_z = [data_sphere.points[3,i] for i in 1:data_sphere.np if colors[i]== 0]
    plot!(noise_x,noise_y,noise_z,seriestype=:scatter,legend = false ,showaxis = false, color = "black", alpha = 0.3, ms= 4)
    c_aux = [ c./norm(c) for c in centers]
    c_aux_x = [c[1] for c in c_aux]
    c_aux_y = [c[2] for c in c_aux]
    c_aux_z = [c[3] for c in c_aux]
    plot!(c_aux_x ,c_aux_y ,c_aux_z ,seriestype=:scatter,legend = false, showaxis = false, color = "black", pch = 3, alpha = 1, ms = 6)
end

#----------------------------------------------------------
# True Signature
#----------------------------------------------------------
export squared_distance_matrix

function squared_distance_matrix(
    data1::SphereData,
    data2::SphereData
)

    @assert sum([!isapprox(sq_norm,1) for sq_norm in sum(data1.points.^2, dims = 1)])==0 "The points should be of norm 1."
    @assert sum([!isapprox(sq_norm,1) for sq_norm in sum(data2.points.^2, dims = 1)])==0 "The points should be of norm 1."

    return [[sq_distance_sphere(data1.points[:,i],data2.points[:,j]) for i in 1:data1.np] for j in 1:data2.np]

end


function squared_distance_matrix(#squared_distance_matrix
    data::SphereData
)

    @assert sum([!isapprox(sq_norm,1) for sq_norm in sum(data.points.^2, dims = 1)])==0 "The points should be of norm 1."

    return [sq_distance_sphere(view(data.points,:,i),data.points) for i in 1:data.np]

end


function sq_distance_sphere(
    x,
    y
)

    return (acos(min(max(sum(x.*y),-1),1)))^2


end

function sq_distance_sphere(
    x::SubArray{T, 1, Matrix{T}, Tuple{Base.Slice{Base.OneTo{S}}, S}, true},#Vector{<:Real},
    vect_y::Matrix{T}
) where {T<:Real} where {S<:Real}

    return vec((acos.(min.(max.(sum(x.*vect_y, dims = 1),-1),1))).^2)

end

"""
function squared_distance_matrix(
    data::SphereData
)

    @assert sum([!isapprox(sq_norm,1) for sq_norm in sum(data.points.^2, dims = 1)])==0 "The points should be of norm 1."

    return [[sq_distance_sphere(data.points[:,i],data.points[:,j]) for i in 1:data.np] for j in 1:data.np]

end
"""