using Pkg
Pkg.activate(".")
Pkg.update()
using RCWA
Nx=3
Ny=3
λ=1000
model=Model([PlainLayer(100,Constant(2)),PatternedLayer(200,[Constant(3+1im),Constant(2)],[Circle(.5)])],Constant(1),Constant(4))
grd=srcwa_grid(model,Nx,Ny,λ,1E-5,0,1000,1000)
mtr=srcwa_matrices(model,grd,λ)
a0te,a0tm=srcwa_source(grd.kin,Nx,Ny)
R,T=srcwa_reftra(a0te,grd::Srcwa_grid,mtr::Srcwa_matrices)
println(R)
println(T)
a,b=srcwa_amplitudes(a0te,grd::Srcwa_grid,mtr::Srcwa_matrices)
A=srcwa_abs(a,b,grd)
println(A)
println(A[1]-A[end])
