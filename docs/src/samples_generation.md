# Generating samples on the circle, the sphere, the torus, the Grassmannian and the Bolza Surface

```@example samplegen
using HomogeneousManifoldTesting
using LinearAlgebra: Diagonal
```

## 1) Circle

```@example samplegen
a = CircleData(500,"uniform")
plot_data(a)
```

```@example samplegen
κ = 2
center_ = [sqrt(2)/2,sqrt(2)/2]
a = CircleData(500,"FvM",κ,center_)
plot_data(a)
```

```@example samplegen
κ = 10
a = CircleData(500,"MFvM",κ)
plot_data(a)
```

```@example samplegen
a = CircleData(100,"normal",6,0.3)
plot_data(a)
```

```@example samplegen
a = CircleData(100,"mixture normal",[0,1],[0.1,0.2],[0.3,0.7])
plot_data(a)
```

```@example samplegen
a = CircleData(100,"dynamical system",0.2,4)
plot_data(a)
```

```@example samplegen
a = CircleData(100,"dynamical system 2",1000)
plot_data(a)
```

## 2) Sphere

```@example samplegen
a = SphereData(500,"uniform")
plot_data(a)
```

```@example samplegen
κ = 2
center_ = [0,sqrt(2)/2,sqrt(2)/2]
a = SphereData(500,"FvM",κ,center_)
plot_data(a,0.25)
```

```@example samplegen
κ = 10
a = SphereData(500,"MFvM",κ)
plot_data(a)
```

## 3) Torus

```@example samplegen
a = TorusData(500,"uniform")
plot_data(a)
```

```@example samplegen
μ = [0,0.5]
Σ = Matrix(Diagonal([0.3,0.03]))
a = TorusData(1000,"normal",μ,Σ)
plot_data(a)
```

```@example samplegen
vect_μ = [[0,0.5],[0.5,0.3]]
vect_Σ = [Matrix(Diagonal([0.2^2,0.3^2])),Matrix(Diagonal([0.5^2,0.1^2]))]
proba = [0.3,0.7]
a = TorusData(1000,"mixture normal",vect_μ,vect_Σ,proba)
plot_data(a)
```

```@example samplegen
a = TorusData(100,"large circle incremental",200,[0.5,0.5])
plot_data(a)
```

```@example samplegen
a = TorusData(100,"large circle uniform",200,[0.5,0.5])
plot_data(a)
```

```@example samplegen
a = TorusData(100,"large circle equidistant",200,[0.5,0.5])
plot_data(a)
```

```@example samplegen
p0 = [0.5,0.5]
r1 = 0.4
r2 = 1
angle = 0.4
a = TorusData(100,"large ellipse incremental",p0,r1,r2,angle)
plot_data(a)
```

```@example samplegen
p0 = [0.5,0.5]
r1 = 0.4
r2 = 1
angle = 0.4
a = TorusData(100,"large ellipse uniform",p0,r1,r2,angle)
plot_data(a)
```

```@example samplegen
p0 = [0.5,0.5]
r = 0.4
angle = 2.4
angle1 = 4.5
a = TorusData(100,"large circle line incremental",p0,r,angle,angle1)
plot_data(a)
```

```@example samplegen
p0 = [0.5,0.5]
r = 0.4
angle = 2.4
angle1 = 4.5
a = TorusData(100,"large circle line uniform",p0,r,angle,angle1)
plot_data(a)
```

```@example samplegen
alp = 0.241
bet = 0.834
a = TorusData(1000,"dynamical system",alp,bet)
plot_data(a)
```

```@example samplegen
p0 = [0.5,0.5]
a = TorusData(1000,"dynamical system2",p0)
plot_data(a)
```

```@example samplegen
p0 = [0.5,0.5]
a = TorusData(1000,"dynamical system3",p0)
plot_data(a)
```

## 4) Grassmann (2,3)

```@example samplegen
a = GrassmannData23(100,"uniform")
plot_data(a)
```

```@example samplegen
a = GrassmannData23(100,"mixture_normal",0.1)
plot_data(a)
```

## 5) Grassmann (2,4)


The plot function is not relevant. It consists in projecting the two basis vectors onto the coordinates $[1,2,3]$ and then to compute the orthogonal of the spanned plane in $\mathbb{R}^3$.

```@example samplegen
a = GrassmannData24(100,"uniform")
plot_data(a,[1,2,3])
```

```@example samplegen
a = GrassmannData24(100,"mixture_normal",1)
plot_data(a,[1,2,3])
```

## 6) Bolza surface

```@example samplegen
a = BolzaData(500,"uniform")
plot_data(a)
```

```@example samplegen
a = BolzaData(10000,"geodesic",0.02)
plot_data(a)
```

```@example samplegen
a = BolzaData(100,"large circle incremental",200/1000,1000)
plot_data(a)
```

```@example samplegen
a = BolzaData(100,"large circle uniform",200/1000,1000)
plot_data(a)
```

```@example samplegen
a = BolzaData(100,"brownian",0+0.3*im,30,5000)
plot_data(a)
```

```@example samplegen
a = BolzaData(500,"mixture_brownian",[0+0.3*im, 0 - 0.4*im],[1,1],[10,10],[0.3,0.7])
plot_data(a)
```

```@example samplegen
z0,T,nsteps = 0*im, 30, 5000
a,jumps = simulate_brownian_bolza(z0,T,nsteps)
plot_brownian_motion_Bolza(a,jumps,5000,0.2)
```
