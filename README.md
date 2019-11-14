# RCWA.jl

## Modeling

Modeling is carried out layer by layer. 

### Materials

So far, only isotropic materials are supported. RCWA is a frequency domain method, so spectroscopic data can be taken directly and arbitrary permittivities can be entered. One can implement constant permittivities or interpolated spectroscopic data.

```julia
M1=ConstantPerm(4+2im)#Independent of wavelength, the permittivity has a value of 4+2i
```
### Geometry

One can specify the distribution of materials within each layer with simple geometric shapes. Currently, the package implements rectangular and elliptic inclusions. Rotation and translation in the plane is also possible. 

```julia
R=Rectangle(.5,.5) 

### Custom
