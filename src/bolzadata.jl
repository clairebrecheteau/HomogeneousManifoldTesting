#--------------------------------------------------------------------------
# Struct BolzaData
#--------------------------------------------------------------------------

export BolzaData

"""
$(TYPEDEF)

BolzaData to store data on the Bolza surface
- np: number of points
- points: vector containing the coordinates of the np points (size = (2,np))
"""
struct BolzaData{T} <: AbstractData

    np::Int
    points::Matrix{T}

end


#--------------------------------------------------------------------------
# Constructor for BolzaData
#--------------------------------------------------------------------------

"""
$(SIGNATURES)

- np: sample size

Returns a sample of size np, on the Bolza surface.
"""

function BolzaData(
    np::Integer,
    method::String,
    params...
)

    @assert method in ["uniform","geodesic","large circle incremental","large circle uniform","brownian","mixture_brownian"] "Parameter method should be in ['uniform','geodesic','large circle incremental','large circle uniform','brownian','mixture_brownian']"
    if method == "uniform"
        points = random_uniform_bolza_data(np)
    elseif method == "geodesic"
        points = random_geodesic_bolza_data(np,params...)
    elseif method == "large circle incremental"
        points = random_large_circle_incremental_angle_bolza_data(np,params...)
    elseif method == "large circle uniform"
        points = random_large_circle_uniform_angle_bolza_data(np,params...)
    elseif method == "brownian"
        points = random_brownian_bolza_data(np,params...)
    elseif method == "mixture_brownian"
        points =  random_mixture_brownian_bolza_data(np, params...)
    end

    return BolzaData(size(points)[2],points)

end


function random_uniform_bolza_data(
    np::Integer
)

    R = 1/sqrt(sqrt(2))
    x_0 = (1+1/sqrt(2))/(2/sqrt(sqrt(2))*cos(pi/8))
    size = 0
    sample = (0+0*im).*zeros(np)
    while (size < np)
        U = rand(1)[1]
        V = rand(1)[1]
        Z = sqrt(U/(sqrt(2) - 1 + U))*exp(2*im*pi*V)
        dist_to_8_centers = [abs(Z - x_0*exp(2*im*pi/8*k)) for k in 1:8]
        if sum([ 1*(d  < sqrt(x_0^2 - 1)) for d in dist_to_8_centers ]) == 0 
            size +=1
            sample[size] = Z
        end
    end

    points = hcat([[real(sample[i]), imag(sample[i])] for i in 1:np]...)

    return points
    
end

### Auxiliary functions

function Fuchsian_elements(
)
    
    alpha_ = sqrt(sqrt(2) - 1)
    F = [hcat([1 + sqrt(2) ,alpha_*(2 + sqrt(2))*exp(- im*k*pi/4)], [alpha_*(2 + sqrt(2))*exp(im*k*pi/4), 1 + sqrt(2)]) for k in 1:4]
    F2 = [F[1], F[2], F[3], F[4], inv(F[1]), inv(F[2]), inv(F[3]), inv(F[4])]
    
    return F2
    
end

export fuch_elements_

const fuch_elements_ = Fuchsian_elements()

function Fuchsian_compositions( 
    )
# Returns the coefficients $α_{k,k'}$ and $β_{k,k'}$ of the composition of the $k$ and $k'$ Fuchsian elements.

    F = Fuchsian_elements()

    F_ = Vector{Matrix{ComplexF64}}(undef,9)
    F_[1:8] .= F
    F_[9] = hcat([1.0 + 0*im, 0], [0, 1])

    α_β =  Array{ComplexF64}(undef, 2, 81)

    for j in 1:9
        for i in 1:9
            αi = F_[i][1,1]
            βi =  F_[i][1,2]
            αj = F_[j][1,1]
            βj =  F_[j][1,2]
            α_β[1,(j-1)*9+i] = αi*αj + βi*conj(βj)
            α_β[2,(j-1)*9+i] = αi*βj + αj*βi
        end
    end

    indices = [1]

    for i in 2:81
        n_occ = 0
        for j in indices
            n_occ += 1*(isapprox(α_β[1:2,i] ,α_β[1:2,j]))
        end
        if n_occ == 0
            append!(indices, i)
        end
    end

    return α_β[1:2,indices]

end

function omega_(
    z::Complex
)
    
    return (z-im)/(z+im)
    
end
    
function geodesic_transform(
    z::Complex,
    eiθ::Complex,
    t::Real
)
    
    c_ = real(eiθ) 
    s_ = imag(eiθ) 
    x = real(z)
    y = imag(z)
    a = sqrt(y)*c_ - x/sqrt(y)*s_
    b = sqrt(y)*s_ + x/sqrt(y)*c_
    c = - s_ / sqrt(y)
    d = c_ / sqrt(y)
    z_aux = (im*a*exp(t) + b)/(im*c*exp(t) + d)
    eiθ_aux = (d - im*c*exp(t))/abs(d - im*c*exp(t))
    
    return [z_aux, eiθ_aux] 
        
end
    
function fuch_transform(
    z::Complex,
    eiθ::Complex,
    fuch::Matrix{<:Complex}
)
    
    a = fuch[1,1] + real(fuch[2,1])
    b = imag(fuch[2,1])
    c = b
    d = fuch[1,1] - real(fuch[2,1])
    z_aux = (a*z + b)/(c*z + d)
    eiθ_aux = (c*conj(z) + d)/abs(c*conj(z) + d)*eiθ
    
    return [z_aux, eiθ_aux]
    
end



### Brownian Motion

function random_mixture_brownian_bolza_data(
    np::Integer,
    z0::Vector{<:Complex},
    T::Vector{<:Real},
    nsteps::Vector{<:Integer},
    proba::Vector{<:Real}
)

    n_components = length(proba)
    n_per_comp = rand(Multinomial(np, proba))
    points = hcat([random_brownian_bolza_data(n_per_comp[i],z0[i],T[i],nsteps[i]) for i in 1:n_components]...)

    return points
end


function random_brownian_bolza_data(
    np::Integer,
    z0::Complex,
    T::Real,
    nsteps::Integer
)

    sample = [simulate_brownian_bolza(z0,T,nsteps)[1][nsteps+1] for i in 1:np]
    points = collect(transpose(hcat([real.(sample), imag.(sample)]...)))

    return points
end

export simulate_brownian_bolza

function simulate_brownian_bolza(
    z0::Complex,
    T::Real,
    nsteps::Integer,
    fuch = fuch_elements_
)

    # Il y aura des problèmes dans l'algo si le pas est trop gros car on peut sortir de la surface de Bolza trop loin.

    # z0 est le point initial (un complexe)
    # T le temps, on cherche à estimer $B_T$
    # nsteps le nombre d'itérations du MB.

    #fuch = Fuchsian_elements()
    x_0 = (1+1/sqrt(2))/(2/sqrt(sqrt(2))*cos(pi/8))
    r = sqrt(x_0^2 - 1)

    dist_to_8_centers = [abs(z0 - x_0*exp(2*im*pi/8*k + im*pi)) for k in 1:8]

    @assert (abs(z0)<1 && sum([ 1*(d  < r) for d in dist_to_8_centers ]) == 0) "z0 should be in the Bolza Surface"

    sq_h = 0.5*sqrt(T/nsteps)
    chi2_dist = Chisq(2)

    path = 0*im.*zeros(nsteps+1)
    z = z0
    path[1] = z0

    jumps = []

    for step in 2:(nsteps+1)
        
        zaux = z
        
        R = sqrt(rand(chi2_dist,1)[1])
        Θ = 2*pi*rand(1)[1]
        #ξ = sq_h*R*exp(im*Θ)

        move = tanh(sq_h*R)/2*exp(im*Θ) #" 0 vérifier"
        z = (move + zaux )/( conj(zaux)*move + 1)

        # On le renvoie sur la surface de Bolza si jamais ça sort.

        dist_to_8_centers = [abs(z - x_0*exp(2*im*pi/8*k + im*pi)) for k in 1:8] # Changement
        continue_ = (sum([ 1*(d  < r) for d in dist_to_8_centers ]) == 0)
        if !continue_
            push!(jumps, step-1)
            while !continue_
                f = fuch[argmin(dist_to_8_centers)]
                z = (f[1,1]*z + f[1,2])/(f[2,1]*z + f[2,2]) # Changement
                dist_to_8_centers = [abs(z - x_0*exp(2*im*pi/8*k + im*pi)) for k in 1:8]
                continue_ = (sum([ 1*(d  < r) for d in dist_to_8_centers ]) == 0)
            end
        end

        path[step] = z
        
    end

    return [path,jumps]
end


### Geodesics

function random_geodesic_bolza_data(
    np::Integer,
    epsilon::Real = 0.02
)

    Theta = [rand(1)[1]]
    geodesic = simulate_geodesic_bolza(1,epsilon,np,Theta)[2][:,1]
    points = collect(transpose([real.(geodesic) imag.(geodesic)]))

    return points
end

function random_large_circle_incremental_angle_bolza_data(
    np::Integer,
    epsilon::Real = 0.02,
    nsteps::Integer = 1000
)

    U0 = rand(1)[1]
    Theta = [2*pi*i*U0 for i in 1:np]
    points = real.(simulate_geodesic_bolza(np,epsilon,nsteps,Theta)[1])
    
    return points
end

function random_large_circle_uniform_angle_bolza_data(
    np::Integer,
    epsilon::Real = 0.02,
    nsteps::Integer = 1000
)

    Theta = (2*pi).*rand(np)
    points = real.(simulate_geodesic_bolza(np,epsilon,nsteps,Theta)[1])
    
    return points
end

function simulate_geodesic_bolza(
    np::Integer,
    epsilon::Real, # = 0.02
    nsteps::Integer, # = 1000
    Theta::Vector{<:Real},
    fuch = fuch_elements_

)

    @assert (epsilon>0)&&(epsilon<1) "radius should be positive, ...and should be small enough..."

    #fuch = Fuchsian_elements()
    x_0 = (1+1/sqrt(2))/(2/sqrt(sqrt(2))*cos(pi/8))
    r = sqrt(x_0^2 - 1)

    sample = 0*im.*zeros(np)
    dyn_syst = 0*im.*zeros(nsteps,np)    

    for i in 1:np

        z = im
        z0 = z
        eiθ = exp(im*Theta[i])
        eiθ0 = eiθ

        omega = omega_(z)
        dyn_syst[1,i] = omega # omega = 0

        t = 0
        if (nsteps > 1)
            for step_ in 2:nsteps

                t += epsilon
                tr_t = geodesic_transform(z0,eiθ0,t)
                z = tr_t[1]
                omega = omega_(z)
                eiθ = tr_t[2]
                dist_to_8_centers = [abs(omega - x_0*exp(2*im*pi/8*k + im*pi)) for k in 1:8]

                continue_ = (sum([ 1*(d  < r) for d in dist_to_8_centers ]) == 0)
                while !continue_
                    f = fuch[argmin(dist_to_8_centers)]
                    tr = fuch_transform(z,eiθ,f)
                    z = tr[1]
                    omega = omega_(z)
                    eiθ = tr[2]
                    z0 = z
                    eiθ0 = eiθ
                    t = 0
                    dist_to_8_centers = [abs(omega - x_0*exp(2*im*pi/8*k + im*pi)) for k in 1:8]
                    continue_ = (sum([ 1*(d  < r) for d in dist_to_8_centers ]) == 0)
                end

                dyn_syst[step_,i] = omega
            end
        end
        sample[i] = dyn_syst[nsteps,i]
    end

    trajectories = [hcat([[real(sample[i]), imag(sample[i])] for i in 1:np]...), dyn_syst]

    return trajectories

end

#--------------------------------------------------------------------------
# Function to Plot data
#--------------------------------------------------------------------------

"""
$(SIGNATURES)

- data_bolza: data on the Bolza surface
- alpha : transparency parameter for points

Plot the data on the Bolza surface.
"""

function plot_data(
    data_bolza::BolzaData,
    alpha = 0.2
)

    precision_ = 1000
    theta_disque = ((0:precision_)./precision_).*(2*pi)
    Bolza = Bolza_surface(precision_)
    plot(cos.(theta_disque),sin.(theta_disque), aspect_ratio=:equal,legend = false,grid = false,axis = false, color = :black)
    x_in = Bolza[3]
    y_in = Bolza[4]
    for k = 1:7
        plot!(x_in[k], y_in[k], linewidth=3, lc=:black)
    end
    plot!(x_in[8], y_in[8], linewidth=3, lc=:black)
    plot!(data_bolza.points[1,:], data_bolza.points[2,:],seriestype=:scatter, alpha = alpha, color = 2)

end


function Bolza_surface(
    precision_ = 1000
)

    x_0 = (1+1/sqrt(2))/(2/sqrt(sqrt(2))*cos(pi/8))
    radius = sqrt(x_0^2 - 1)
    theta_intersection_D = asin(sqrt(1 - 1/(x_0^2))/sqrt(x_0^2 - 1))
    theta_intersection_in = asin((1/sqrt(sqrt(2))*sin(pi/8))/sqrt(x_0^2 - 1))
    theta = pi .+(((0:precision_)./precision_).*(2*theta_intersection_D)) .- theta_intersection_D
    theta_in = pi .+(((0:precision_)./precision_).*(2*theta_intersection_in)) .- theta_intersection_in
    circle_0 = hcat(radius.*cos.(theta) .+ x_0, radius.*sin.(theta))
    x0 = circle_0[:,1]
    y0 = circle_0[:,2]
    x = [x0.*cos(2*pi*k/8) + y0.*sin(2*pi*k/8) for k in 0:7]
    y = [x0.*sin(2*pi*k/8) - y0.*cos(2*pi*k/8) for k in 0:7]
    circle_0_in = hcat(radius.*cos.(theta_in) .+ x_0, radius.*sin.(theta_in))
    x0_in = circle_0_in[:,1]
    y0_in = circle_0_in[:,2]
    x_in = [x0_in.*cos(2*pi*k/8) + y0_in.*sin(2*pi*k/8) for k in 0:7]
    y_in = [x0_in.*sin(2*pi*k/8) - y0_in.*cos(2*pi*k/8) for k in 0:7]
    return x,y, x_in, y_in

end

export plot_brownian_motion_Bolza

function plot_brownian_motion_Bolza(path,jumps,nsteps,alpha_)
    precision_ = 1000
    theta_disque = ((0:precision_)./precision_).*(2*pi)

    x_0 = (1+1/sqrt(2))/(2/sqrt(sqrt(2))*cos(pi/8))
    r = sqrt(x_0^2 - 1)

    Bolza = Bolza_surface(precision_)
    plot(cos.(theta_disque),sin.(theta_disque), aspect_ratio=:equal)
    x_in = Bolza[3]
    y_in = Bolza[4]
    for k = 1:8
        plot!(x_in[k], y_in[k], linewidth=3, lc=:black)
    end

    if size(jumps)[1] == 0
        plot!(real.(path),imag.(path), alpha = alpha_, linecolor = :black)
    else
        if jumps[1] > 1
            plot!(real.(path[1:jumps[1]]),imag.(path[1:jumps[1]]), alpha = alpha_, linecolor = :black)
        end
        if size(jumps)[1] > 0
            for t in 1:(size(jumps)[1] - 1)
                if jumps[t+1] > (jumps[t]+1)
                    plot!(real.(path[(jumps[t]+1):jumps[t+1]]),imag.(path[(jumps[t]+1):jumps[t+1]]), alpha = alpha_, linecolor = :black)
                end
            end
        end
        if nsteps > jumps[end]
            plot!(real.(path[(jumps[end]+1):end]),imag.(path[(jumps[end]+1):end]), alpha = alpha_, linecolor = :black)
        end
    end
    plot!(cos.(theta_disque),sin.(theta_disque), aspect_ratio=:equal, legend = false, linecolor = :black)
end

#----------------------------------------------------------
# True Signature
#----------------------------------------------------------


"""
function squared_distance_matrix(
    data::BolzaData
)

    x_0 = (1+1/sqrt(2))/(2/sqrt(sqrt(2))*cos(pi/8))
    r = sqrt(x_0^2 - 1)
    radius_insc = x_0 - r

    alpha_ = sqrt(sqrt(2) - 1)
    f11 = 1 + sqrt(2)
    f12 = - alpha_*(2 + sqrt(2))

    # On envoie le premier point que 0, le suivant avec la même transformation
    # On l'envoie sur [0,1] par une rotation
    # puis on le ramène sur la surface de Bolza.
    # Ensuite, on calcule très simplement la distance entre cette image et 0.

    # C'est pas encore tout le temps pareil que la fonction précédente mais bon... peut-être erreur numérique...

    sample_size = size(data.points)[2]
    
    points_complex = view(data.points,1,:) .+ im*view(data.points,2,:)
    
    dist_matrix = Array{eltype(data.points)}(undef,sample_size,sample_size)
    for j in 1:sample_size
        dist_matrix[j,j] = 0
        for i in 1:(j - 1)
            z = abs((points_complex[i] - points_complex[j])/(1 - points_complex[i]*conj(points_complex[j]))) # On l'envoie sur [0,1] par une rotation.
            while z > radius_insc
                z = abs((f11*z + f12)/(f12*z + f11))
            end
            dist_matrix[i,j] = (2*atanh(z))^2
            dist_matrix[j,i] = dist_matrix[i,j]
        end
    end
    
    return [dist_matrix[i,:] for i in 1:sample_size]

end

#
#
# ATTENTION : POUR L'INSTANT, NE FONCTIONNE PAS !!!!
# La transformation n'est peut-être pas une isométrie. Auquel cas il faudra peut-être aussi modifier la fonction géodésique et brownien...
#
#
function squared_distance_matrix_new(
    data::BolzaData,
    fuch = fuch_elements_
)

    x_0 = (1+1/sqrt(2))/(2/sqrt(sqrt(2))*cos(pi/8))
    r = sqrt(x_0^2 - 1)

    # On envoie le premier point que 0, le suivant avec la même transformation
    # On l'envoie sur [0,1] par une rotation
    # puis on le ramène sur la surface de Bolza.
    # Ensuite, on calcule très simplement la distance entre cette image et 0.

    sample_size = size(data.points)[2]
    
    points_complex = view(data.points,1,:) .+ im*view(data.points,2,:)
    
    dist_matrix = Array{eltype(data.points)}(undef,sample_size,sample_size)
    for j in 1:sample_size
        dist_matrix[j,j] = 0
        for i in 1:(j - 1)
            z = (points_complex[i] - points_complex[j])/(1 - points_complex[i]*conj(points_complex[j]))
            dist_to_8_centers = [abs(z - x_0*exp(2*im*pi/8*k + im*pi)) for k in 1:8] # Changement
            continue_ = (sum([ 1*(d  < r) for d in dist_to_8_centers ]) == 0)
            while !continue_
                f = fuch[argmin(dist_to_8_centers)]
                z = (f[1,1]*z + f[1,2])/(f[2,1]*z + f[2,2]) # Changement
                dist_to_8_centers = [abs(z - x_0*exp(2*im*pi/8*k + im*pi)) for k in 1:8]
                continue_ = (sum([ 1*(d  < r) for d in dist_to_8_centers ]) == 0)
            end

            dist_matrix[i,j] = (2*atanh(abs(z)))^2
            dist_matrix[j,i] = dist_matrix[i,j]
        end
    end
    
    return [dist_matrix[i,:] for i in 1:sample_size]

end


"""


function squared_distance_matrix(
    data::BolzaData
)
    # Distance minimale entre un point et les 81 transformations possibles de l'autre point (Fuchsian + identité)
    # Petite amélioration car on se compare à 81 points (il y a redondance mais ça va...).
    # On peut améliorer cette fonction en calculant d'abord les produit des matrices et en enlevant les matrices égales (moins de 81 calculs comme ça...) On gagnerait.
    # Créer la fonction en faisant trois produits... pour améliorer encore (mais ça risque de devenir assez long...)
    sample_size = size(data.points)[2]

    α_β = Fuchsian_compositions() # Matrix of size 2x65 that contains all pairs(alpha,beta) for the matrices (alpha, beta, bar beta, bar alpha) that are composition of two elements of Fuchsian + the identity.
    n_transfo = size(α_β)[2]

    points_complex = view(data.points,1,:) .+ im*view(data.points,2,:)

    points_transformed_complex = Array{eltype(points_complex)}(undef,n_transfo,sample_size)

    for i in 1:n_transfo
        points_transformed_complex[i,:] = ( (α_β[1,i] .*points_complex) .+ α_β[2,i])./( (conj(α_β[2,i]) .*points_complex) .+ conj(α_β[1,i]))
    end

    real_ptc = real.(points_transformed_complex)
    imag_ptc = imag.(points_transformed_complex)

    dist_matrix = Array{eltype(data.points)}(undef,sample_size,sample_size)
    for j in 1:sample_size
        dist_matrix[j,j] = 0
        for i in 1:(j - 1)
            dist_matrix[i,j] = minimum(distance_Poincare_disk(view(data.points,:,i),view(real_ptc,:,j),view(imag_ptc,:,j)))^2
            dist_matrix[j,i] = dist_matrix[i,j]
        end
    end

    return [dist_matrix[i,:] for i in 1:sample_size] #[view(dist_matrix,i,:) for i in 1:sample_size]
end


function  distance_Poincare_disk(
    p::SubArray{T, 1, Matrix{T}, Tuple{Base.Slice{Base.OneTo{S}}, S}, true},#Vector{<:Real},
    x_vect::SubArray{T, 1, Matrix{T}, Tuple{Base.Slice{Base.OneTo{S}}, S}, true},#Vector{<:Real},
    y_vect::SubArray{T, 1, Matrix{T}, Tuple{Base.Slice{Base.OneTo{S}}, S}, true}#Vector{<:Real}
) where {T <: Real} where {S <: Real}# Returns a vector that contains the distances of p to the points (x[i],y[i])  

    C = 1 - norm(p)^2
    norm_q_sq =  x_vect.^2 + y_vect.^2
    norm_pq_sq = (x_vect.-p[1]).^2 + (y_vect.-p[2]).^2
    
    return acosh.( 1 .+ (2 .* norm_pq_sq ./ (C.* (1 .- norm_q_sq)))) 
end


function distance_Poincare_disk(
    p,
    q
) # Validée ! C'est la formule de Manifold.jl, ça correspond aussi à la formule qu'on obtient avec int sqrt(g(gamma.,gamma.))
    
    return(acosh( 1 + 2 * norm(p .- q)^2 / ((1 - norm(p)^2) * (1 - norm(q)^2)))) # LA BONNE
 
end

function sq_distance_matrix_Bolza(points::Matrix{Float64})
    sample_size = size(points)[2]
    return [[distance_Poincare_disk(points[:,i],points[:,j])^2 for i in 1:sample_size] for j in 1:sample_size]
end


"""
function squared_distance_matrix_old(
    data::BolzaData
)
    # Distance minimale entre un point et les 81 transformations possibles de l'autre point (Fuchsian + identité)
    # Petite amélioration car on se compare à 81 points (il y a redondance mais ça va...).
    # On peut améliorer cette fonction en calculant d'abord les produit des matrices et en enlevant les matrices égales (moins de 81 calculs comme ça...) On gagnerait.
    # Créer la fonction en faisant trois produits... pour améliorer encore (mais ça risque de devenir assez long...)
    sample_size = size(data.points)[2]
    fuch = Fuchsian_elements()

    points_complex = data.points[1,:] .+ im*data.points[2,:]
    points_transformed_complex = [[(f[1,1]*q + f[1,2])/(f[2,1]*q + f[2,2]) for q in points_complex] for f in fuch]
    points_transformed_complex = [ points_transformed_complex ; [points_complex] ]
    points_transformed_complex_bis = vcat([[points_transformed_complex] ; [[[(f[1,1]*q + f[1,2])/(f[2,1]*q + f[2,2]) for q in pts] for pts in points_transformed_complex] for f in fuch]]...)
    points_transformed = [transpose(hcat([real.(pts),imag.(pts)]...)) for pts in points_transformed_complex_bis]

    return [[min([distance_Poincare_disk(data.points[:,i],points_transformed[nfuch][:,j]) for nfuch in 1:81]...)^2 for i in 1:sample_size] for j in 1:sample_size]
end

function squared_distance_matrix_old_bis(
    data::BolzaData
)
    # Distance minimale entre un point et les 81 transformations possibles de l'autre point (Fuchsian + identité)
    # Petite amélioration car on se compare à 81 points (il y a redondance mais ça va...).
    # On peut améliorer cette fonction en calculant d'abord les produit des matrices et en enlevant les matrices égales (moins de 81 calculs comme ça...) On gagnerait.
    # Créer la fonction en faisant trois produits... pour améliorer encore (mais ça risque de devenir assez long...)
    sample_size = size(data.points)[2]

    α_β = Fuchsian_compositions() # Matrix of size 2x65 that contains all pairs(alpha,beta) for the matrices (alpha, beta, bar beta, bar alpha) that are composition of two elements of Fuchsian + the identity.
    n_transfo = size(α_β)[2]

    points_complex = data.points[1,:] .+ im*data.points[2,:]

    points_transformed_complex = Array{eltype(points_complex)}(undef,n_transfo,sample_size)

    for i in 1:n_transfo
        points_transformed_complex[i,:] = ( (α_β[1,i] .*points_complex) .+ α_β[2,i])./( (conj(α_β[2,i]) .*points_complex) .+ conj(α_β[1,i]))
    end

    real_ptc = real.(points_transformed_complex)
    imag_ptc = imag.(points_transformed_complex)

    dist_matrix = Array{eltype(data.points)}(undef,sample_size,sample_size)
    for j in 1:sample_size
        dist_matrix[j,j] = 0
        for i in 1:(j - 1)
            dist_matrix[i,j] = minimum([distance_Poincare_disk(data.points[:,i],[real_ptc[nfuch,j],imag_ptc[nfuch,j]]) for nfuch in 1:n_transfo])^2
            dist_matrix[j,i] = dist_matrix[i,j]
        end
    end

    return [dist_matrix[i,:] for i in 1:sample_size]
end

 function distance_Bolza(p,q)
     # Pas précise car on ne se compare qu'aux images par les matrices de Bolza du point. On gagne à se comparer aux images des images.
     
     fuch = Fuchsian_elements()
     p = p[1] + im*p[2]
     q = q[1] + im*q[2]
     
     return min(distance_Poincare_disk(p,q),min([acosh( 1 + 2 * abs(p - (f[1,1]*q + f[1,2])/(f[2,1]*q + f[2,2]))^2 / ((1 - abs(p)^2) * (1 - abs((f[1,1]*q + f[1,2])/(f[2,1]*q + f[2,2]))^2))) for f in fuch]...))
 
 end
 """
 
 
 
 
 
 # Modified version of Pierre Navaro :
 
"""
function squared_distance_matrix(
    data::BolzaData
)
    # Distance minimale entre un point et les 81 transformations possibles de l'autre point (Fuchsian + identité)
    # Petite amélioration car on se compare à 81 points (il y a redondance mais ça va...).
    # On peut améliorer cette fonction en calculant d'abord les produit des matrices et en enlevant les matrices égales (moins de 81 calculs comme ça...) On gagnerait.
    # Créer la fonction en faisant trois produits... pour améliorer encore (mais ça risque de devenir assez long...)
    sample_size = size(data.points, 2)

    α_β = Fuchsian_compositions() # Matrix of size 2x65 that contains all pairs(alpha,beta) for the matrices (alpha, beta, bar beta, bar alpha) that are composition of two elements of Fuchsian + the identity.
    n_transfo = size(α_β, 2)

    points_complex = view(data.points,1,:) .+ im*view(data.points,2,:)

    points_transformed_complex = Array{eltype(points_complex)}(undef,n_transfo,sample_size)

    for i in 1:n_transfo
        points_transformed_complex[i,:] .= ( (α_β[1,i] .* points_complex) .+ α_β[2,i])./( (conj(α_β[2,i]) .*points_complex) .+ conj(α_β[1,i]))
    end

    real_ptc = real.(points_transformed_complex)
    imag_ptc = imag.(points_transformed_complex)

    dist_matrix = Array{eltype(data.points)}(undef,sample_size,sample_size)
    norm_q_sq =  zeros(eltype(data.points), sample_size)
    norm_pq_sq = zeros(eltype(data.points), sample_size)

    for j in axes(dist_matrix, 2)
        dist_matrix[j,j] = 0
        for i in 1:(j - 1)
            dist_matrix[i,j] = minimum(distance_Poincare_disk(norm_q_sq, norm_ps_sq, view(data.points,:,i),view(real_ptc,:,j),view(imag_ptc,:,j)))^2
            dist_matrix[j,i] = dist_matrix[i,j]
        end
    end

    return [dist_matrix[i,:] for i in 1:sample_size] #[view(dist_matrix,i,:) for i in 1:sample_size]
end

function  distance_Poincare_disk(norm_q_sq, norm_pq_sq,
    p::SubArray{T, 1, Matrix{T}, Tuple{Base.Slice{Base.OneTo{S}}, S}, true},#Vector{<:Real},
    x_vect::SubArray{T, 1, Matrix{T}, Tuple{Base.Slice{Base.OneTo{S}}, S}, true},#Vector{<:Real},
    y_vect::SubArray{T, 1, Matrix{T}, Tuple{Base.Slice{Base.OneTo{S}}, S}, true}#Vector{<:Real}
) where {T <: Real} where {S <: Real}# Returns a vector that contains the distances of p to the points (x[i],y[i])  

    C = 1 - norm(p)^2
    norm_q_sq .=  x_vect.^2 .+ y_vect.^2
    norm_pq_sq .= (x_vect.-p[1]).^2 .+ (y_vect.-p[2]).^2
    
    return acosh.( 1 .+ (2 .* norm_pq_sq ./ (C.* (1 .- norm_q_sq)))) 
end
"""
 
 

