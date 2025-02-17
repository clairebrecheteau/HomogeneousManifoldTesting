# Signatures computation

```@example signatures
using HomogeneousManifoldTesting
using Plots
```

```@example signatures
sample_size = 100
mc_size = 10000
vect_prop_near_neigh = [0.1];
```

## 1) Circle

```@example signatures
test_circle = Test_Homogeneity(sample_size,vect_prop_near_neigh,mc_size,CircleData);
```

```@example signatures
test_circle.true_signatures
```

```@example signatures
test_circle.barycenter_signatures
```

```@example signatures
x = (1:length(test_circle.barycenter_signatures[1]))/length(test_circle.barycenter_signatures[1])
plot(x,test_circle.barycenter_signatures[1], label = "barycenter signature", title = "Uniform distribution on the circle \n prop. nearest neighbours h = $(vect_prop_near_neigh[1]), sample size = $sample_size")
plot!(x,test_circle.true_signatures .* ones(length(test_circle.barycenter_signatures[1])), label =  "true signature")
```

## 2) Grassmann (2,3)

```@example signatures
mc_size = 100
```

```@example signatures
test_grassmann23 = Test_Homogeneity(sample_size,vect_prop_near_neigh,mc_size,GrassmannData23);
```

```@example signatures
# Not correct
test_grassmann23.true_signatures
```

```@example signatures
test_grassmann23.barycenter_signatures
```

```@example signatures
x = (1:length(test_grassmann23.barycenter_signatures[1]))/length(test_circle.barycenter_signatures[1])
plot(x,test_grassmann23.barycenter_signatures[1], label = "barycenter signature", title = "Uniform distribution on the Grassmannian (2,3) \n prop. nearest neighbours h = $(vect_prop_near_neigh[1]), sample size = $sample_size")
plot!(x,test_grassmann23.true_signatures .* ones(length(test_circle.barycenter_signatures[1])), label =  "true signature")
```

!!! warning "CAUTION"
    The true signature is not available yes for Grassmannian (2,3) and  Grassmannian (2,4). The default value 0 is provided. Therefore the test option "homogeneity" is not correct. Only the test option "iidness" should be used. 
