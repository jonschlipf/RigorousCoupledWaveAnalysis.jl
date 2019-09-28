using Pkg
Pkg.activate(".")
Pkg.update()
using LayeredPhotonics
Nx=1
Ny=1
model=Model([PatternedLayer(1,[Constant(4),Constant(1)],[Circle(.5)])],Constant(4),Constant(1))
grd,mtr=srcwa_matrices(model,Nx,Ny,1,.000001,0,1,2)
a0te,a0tm=srcwa_source(grd.kin,Nx,Ny)
R,T=srcwa_reftra(a0te,grd::Srcwa_grid,mtr::Srcwa_matrices)
