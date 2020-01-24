# RCWA.jl

This implements both the scattering matrix and the enhanced transmission matrix RCWA algorithms in julia for periodic multilayer structures in nano-optics and RF.

## Modeling

Modeling is carried out layer by layer. 

### Materials

So far, only isotropic materials with are supported. RCWA is a frequency domain method, so spectroscopic data can be taken directly and arbitrary permittivities can be entered. One can implement constant permittivities or interpolated spectroscopic data.

```julia
Air=ConstantPerm(1)#air layer with relative permittivity of 1
M1=ConstantPerm(4+2im)#Independent of wavelength, the permittivity has a value of 4+2i
wavelength=600:100:1600
e=[3,4,5,6,5,4,3,2,3,4,5]
E=interpolate((wavelength,), e, Gridded(Linear())
M2=InterpolPerm(E)
```
### Geometry

One can specify the distribution of materials within each layer with simple geometric shapes. Currently, the package implements rectangular and elliptic inclusions. Rotation and translation in the plane is also possible. All coordinates are relative to the cell size in the respective direction.

```julia
R=Rectangle(.2,.2) #create rectangle with width and height one fifth of the cell size
E=Ellipse(.1,.3) #create ellipse with relative radius 0.1 along the x axis and 0.3 along the y axis
R=Rotation(R,pi/4) #rotate the rectangle in the plane by 45 degrees
E=Shift(E,.8,.1) #shift the ellipse in x-direction by 0.8 and in y-direction by 0.1
Geo=Combination([R,E]) #combine the ellipse and rectangle into one geometry object
```

There are structures implemented for plain (homogenous) layers and patterned layers. A plain layer requires just a thickness and a material object. A patterned layer is defined by its thickness, an array of materials (at least one) and an array of geometry objects (size should be one smaller than that of the materials array). 

```julia
L1=PlainLayer(100,M1) #homogenous layer of thickness 100
L2=PatternedLayer(200,[M1,M2],[Geo]) #patterned layer of thickness 100 with inclusion of material M1 in a background of M2
```

A model object requires an arrays of layers (sorted in direction of the light propagation) and materials for the superstrate and substrate halfspaces.

```julia
Mdl=RCWAModel([L1,L2],Air,M2) #device with two layers in air, M2 is the substrate
```
## References

1. D. M. Whittaker and I. S. Culshaw, Scattering-matrix treatment of patterned multilayer photonic structures, Phys. Rev. B60(1999), 2610–2618.1

2. Marco Liscidini, Dario Gerace, Lucio Claudio Andreani, and J. E. Sipe, Scattering-matrix analysis of periodically patterned multilayers with asymmetric unit cells and birefringent media,Phys. Rev. B77(2008), 035324.1

3. Raymond Rumpf, Improved formulation of scattering matrices for semi-analytical methods thatis consistent with convention, Progress In Electromagnetics Research B35(2011), 241–261.1
