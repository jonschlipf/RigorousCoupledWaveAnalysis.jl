# RigorousCoupledWaveAnalysis.jl - Rigorous Coupled-Wave Analysis (RCWA)

This implements both the scattering matrix and the enhanced transmission matrix RCWA algorithms in julia for periodic multilayer structures in nano-optics and RF.

| Status | Coverage |
| :----: | :----: |
| [![Build Status](https://travis-ci.com/jonschlipf/RigorousCoupledWaveAnalysis.jl.svg?branch=master)](https://travis-ci.com/jonschlipf/RigorousCoupledWaveAnalysis.jl) | [![codecov.io](http://codecov.io/gh/jonschlipf/RigorousCoupledWaveAnalysis.jl/coverage.svg?branch=master)](http://codecov.io/gh/jonschlipf/RigorousCoupledWaveAnalysis.jl?branch=master) |

## Modeling


### Materials

The package supports isotropic materials along with diagonal and in-plane anisotropy so far. RCWA is a frequency domain method, so spectroscopic data can be taken directly and arbitrary permittivities can be entered. One can implement constant permittivities, dispersion formulas or interpolated spectroscopic data.

```julia
Air=ConstantPerm(1) #air material with relative permittivity of 1
M1=ConstantPerm(4+2im) #independent of wavelength, the permittivity has a value of 4+2i
wavelength=600:100:1600  #wavelength axis for interpolated permittivity data
e=[3,4,5,6,5,4,3,2,3,4,5] #permittivity data to be interpolated
using Interpolations
E=interpolate((wavelength,), e, Gridded(Linear()))
M2=InterpolPerm(E) #model with interpolated permittivity
```
Some literature materials are already included by default. More can be incorporated on request.
```julia
Si=InterpolPerm(RigorousCoupledWaveAnalysis.si_schinke) #Si from interpolated literature values
Ge=InterpolPerm(RigorousCoupledWaveAnalysis.ge_nunley) #Ge from interpolated literature values
SiO2=ModelPerm(RigorousCoupledWaveAnalysis.sio2_malitson) #SiO2 from literature dispersion formula
Al=ModelPerm(RigorousCoupledWaveAnalysis.al_rakic) #Al from literature dispersion formula
```
### Geometry

One can specify the distribution of materials within each layer with simple geometric shapes. Currently, the package analytically implements rectangular and elliptic inclusions. Rotation and translation in the plane is also possible. All coordinates are relative to the cell size in the respective direction. For simple verification of the geometry design, the RigorousCoupledWaveAnalysis.drawable function yields the x and y coordinates of the outline. 

```julia
R=Rectangle(.2,.2) #create rectangle with width and height one fifth of the cell size
E=Ellipse(.1,.3) #create ellipse with relative radius 0.1 along the x axis and 0.3 along the y axis
R=Rotation(R,pi/4) #rotate the rectangle in the plane by 45 degrees
E=Shift(E,.8,.1) #shift the ellipse in x-direction by 0.8 and in y-direction by 0.1
Geo=Combination([R,E]) #combine the ellipse and rectangle into one geometry object
Cir=Circle(480/950) #a circle with a diameter of 480 nm in a unit cell with a pitch of 950 nm
x,y=RigorousCoupledWaveAnalysis.drawable(Geo)
```

It is also possible to compute the RCWA for arbitrary structures defined by a bit mask using the Fourier transform. Here, a reciprocal space grid is required before modeling. See the section on grids below for guidelines how to choose the grid order.

```julia
N=4
nx,ny,dnx,dny=ngrid(N,N) #define a grid of reciprocal space, with maximum spatial frequency N
using Random
f=bitrand(100,100) #define the geometry of the unit cell as a random 10x10 bit mask
F=real2recip(dnx,dny,f) #Fourier transform the geometry to the reciprocal space grid
Geo=Custom(F) #custom geometry object with the defined structure
```
### Layers

There are structures implemented for simple (=homogenous) layers and patterned layers. A plain layer requires just a thickness and a material object. A patterned layer is defined by its thickness, an array of materials (at least one) and an array of geometry objects (the size of the geometry array should be one smaller than that of the materials array).

```julia
nha=PatternedLayer(100,[Al,Air],[Cir]) # the a nanohole in the Al layer is filled with air
spa=SimpleLayer(50,SiO2)
nsi=SimpleLayer(20,Si)
nge=SimpleLayer(20,Ge)
ige=SimpleLayer(480,Ge)
```
### Building models

A model object requires an arrays of layers (sorted in direction of the light propagation) and the materials for the superstrate and substrate halfspaces.

```julia
Mdl=RCWAModel([nha,spa,nsi,nge,ige],Air,Si) # a nanohole array device with the layers defined as in the previous section on a Si substrate
```

### Computation grid

RCWA computations are carried out in discretized reciprocal space, so a grid is required. The grid can be discretized by 2π/a, where a is the lattice constant of the 2D unit cell. One has to specify the direction of the impinging plane wave by θ and α, as well as the lattice constant in x and y ax and ay. The only parameter in RCWA that affects accuracy vs performance is N, the maximum spatial frequency to be considered. N=4 is normally a good value for all-dielectric metasurfaces, plasmonics requires higher N.

```julia
N=4 #maximum spatial frequency, same for x and y
λ=1000 #nm wavelength
θ=1E-5 #elevation angle, zero will yield a singularity inversion error
α=0 #azimuth angle
ax=ay=500 #500 nm square cell
λ=1000 #wavelength
Grd=rcwagrid(N,N,ax,ay,θ,α,λ,Air) #create the grid, superstrate is air
```

### Solution

One can employ the enhanced transmission matrix (etm) approach to solve the Maxwell equations for their system. This will yield the reflected and transmitted power.

```julia
ste,stm=rcwasource(Grd)    #create a source object
Rte,Tte=etm_reftra(ste,Mdl,Grd,λ) #run the etm algorithm for TE polarization
Rtm,Ttm=etm_reftra(stm,Mdl,Grd,λ) #run the etm algorithm for TM polarization
```
The scatter matrix method can be called in the same manner.
```julia
ste,stm=rcwasource(Grd)    #create a source object
Rte,Tte=srcwa_reftra(ste,Mdl,Grd,λ) #run the srcwa algorithm for TE polarization
Rtm,Ttm=srcwa_reftra(stm,Mdl,Grd,λ) #run the srcwa algorithm for TM polarization
```
Absorptions within the layers can be computed from the power flows between them
```julia
ste,stm=rcwasource(Grd)    #create a source object
Rte,Tte,fte=etm_reftra_flows(ste,Mdl,Grd,λ) #TE
Ate=-fte[end]-Tte #absorption in lowest layer
Rtm,Ttm,ftm=etm_reftra_flows(stm,Mdl,Grd,λ) #TM
Ate=-ftm[end]-Ttm #absorption in lowest layer
```
### Local fields

Local electric and magnetic fields are obtainable via the propagating amplitudes as well:
```julia
em=RigorousCoupledWaveAnalysis.eigenmodes(Grd,λ,ly)      #get the eigenmodes of propagation in the first layer (this is the nanohole array)
a,b=etm_amplitudes(ste,Mdl,Grd,λ) #get propagating wave amplitudes inside layer
points=[100,100,10]               #set the number of points to compute in x,y,z
E,H=RigorousCoupledWaveAnalysis.getfields(a[1],b[1],Mdl.layers[1].thickness,em,Grd,points,λ) #compute the electric and magnetic field
```
Or via scattering matrix algorithm:
```julia
using LinearAlgebra
em=RigorousCoupledWaveAnalysis.eigenmodes(Grd,λ,ly)        #get the eigenmodes of propagation in the first layer (this is the nanohole array)
a,b=srcwa_amplitudes(ste,Mdl,Grd,λ) #get propagating wave amplitudes outside layer
ain,bout=slicehalf(.5*[em.W\I+em.V\Grd.V0 em.W\I-em.V\Grd.V0;em.W\I-em.V\Grd.V0 em.W\I+em.V\Grd.V0]*[a[:,1];b[:,1]]) #get propagating wave amplitudes inside layer
points=[100,100,10]               #set the number of points to compute in x,y,z
E,H=RigorousCoupledWaveAnalysis.getfields(a[1],b[1],Mdl.layers[1].thickness,em,Grd,points,λ) #compute the electric and magnetic field
```


## References

D. M. Whittaker and I. S. Culshaw, Scattering-matrix treatment of patterned multilayer photonic structures, Phys. Rev. B60(1999), 2610–2618.1

Marco Liscidini, Dario Gerace, Lucio Claudio Andreani, and J. E. Sipe, Scattering-matrix analysis of periodically patterned multilayers with asymmetric unit cells and birefringent media,Phys. Rev. B77(2008), 035324.1

Raymond Rumpf, Improved formulation of scattering matrices for semi-analytical methods thatis consistent with convention, Progress In Electromagnetics Research B35(2011), 241–261.1
