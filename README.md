# RCWA.jl

## Modeling

Modeling is carried out layer by layer. 

### Materials

So far, only isotropic materials are supported. RCWA is a frequency domain method, so spectroscopic data can be taken directly and arbitrary permittivities can be entered. One can implement constant permittivities or interpolated spectroscopic data.

```julia
M1=ConstantPerm(4+2im)#Independent of wavelength, the permittivity has a value of 4+2i
```

### Custom
