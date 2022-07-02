
using LinearAlgebra
ftype=Float64
export ngrid,kgrid,rcwagrid,modes_freespace,RCWAGrid,rcwasource
"""
    RCWAGrid(Nx,Ny,nx,ny,dnx,dny,px,py,Kx,Ky,k0,V0,Kz0)
Structure to store a grid for RCWA computation
# Attributes
* `Nx` : Maximum order in x
* `Ny` : Maximum order in y
* `nx` : Order in x
* `ny` : Order in y
* `dnx` : Differential order in x
* `dny` : Differential order in y
* `px` : Structure pitch in x
* `py` : Structure pitch in y
* `Kx` : Wave vector in x
* `Ky` : Wave vector in y
* `k0` : Free-space 0th-order wavevector
* `V0` : Magnetic eigenstate of free space
* `Kz0` : Wave vector in y, free space
"""
struct RCWAGrid
	Nx::Int64
	Ny::Int64
    nx::Array{Int64,1}
    ny::Array{Int64,1}
    dnx::Array{Int64,2}
    dny::Array{Int64,2}
    px::ftype
    py::ftype
    Kx::Diagonal{Complex{ftype},Array{Complex{ftype},1}}
    Ky::Diagonal{Complex{ftype},Array{Complex{ftype},1}}
    k0::Array{ftype,1}
    V0::AbstractArray{Complex{ftype},2}
    Kz0::Diagonal{Complex{ftype},Array{Complex{ftype},1}}
end
"""
	ngrid(Nx,Ny)
Computes the orders n of the plane wave expansion
# Arguments
* `Nx` : Maximum order in x
* `Ny` : Maximum order in y
# Outputs
* `nx` : Order in x
* `ny` : Order in y
* `dnx` : Differential order in x
* `dny` : Differential order in y
"""
function ngrid(Nx::Int64,Ny::Int64)
    #reciprocal space coordinate indices
    nx=vec([r  for r in -Nx:Nx, c in -Ny:Ny])
    ny=vec([c  for r in -Nx:Nx, c in -Ny:Ny])
    #difference matrix for 2dft
    dnx=[nx[a]-nx[b] for a in 1:length(nx),b in 1:length(nx)]
    dny=[ny[a]-ny[b] for a in 1:length(ny),b in 1:length(ny)]
    return nx,ny,dnx,dny
end
"""
    kgrid(nx,ny,px,py,θ,α,λ)
Computes the wave vectors of the plane wave expansion
# Arguments
* `nx` : Order in x
* `ny` : Order in y
* `px` : Structure pitch in x
* `py` : Structure pitch in y
* `θ` : inclination angle of incoming wave
* `α` : azimuth angle of incoming wave
* `λ` : wavelength
# Outputs
* `Kx` : Wave vector in x
* `Ky` : Wave vector in y
* `k0` : Free-space 0th-order wavevector
"""
function kgrid(nx::Array{Int64,1},ny::Array{Int64,1},px::Real,py::Real,θ::Real,α::Real,λ::Real,sup)
    #all k vectors are generally normalized to k0 here
    #The incoming wave, transformed from spherical coordinates to normalized cartesian coordinates, |k0|=1
    k0=[sin(θ*π/180)*cos(α*π/180),sin(θ*π/180)*sin(α*π/180),cos(θ*π/180)]*real(sqrt(get_permittivity(sup,λ)))
    #the spatial vectors are the sum of the incoming vector and the reciprocal lattice vectors
    kx=k0[1].+real(λ)*nx/px;
    ky=k0[2].+real(λ)*ny/py;
    #need matrix for later computations
    Kx=Diagonal(kx)
    Ky=Diagonal(ky)
    return Complex.(Kx),Complex.(Ky),k0
end
"""
    modes_freespace(Kx,Ky)
Computes the eigenmodes of propagation through free space, for normalization
# Arguments
* `Kx`: x-axis component of the propagation vector
* `Ky`: y-axis component of the propagation vector
# Outputs
* `V0`: coordinate transform between free space eigenmode amplitude and magnetic field
* `Kz0`: z-axis component of the propagation vector in free space
"""
function modes_freespace(Kx::Diagonal{Complex{ftype},Array{Complex{ftype},1}},Ky::Diagonal{Complex{ftype},Array{Complex{ftype},1}})
    #just because |k|=1
    Kz0=sqrt.(Complex.(I-Kx*Kx-Ky*Ky))
	Kz0[imag.(Kz0).<0].*=-1

    #P0 is identity
    Q0=[Kx*Ky I-Kx*Kx;Ky*Ky-I -Ky*Kx]
    #propagation


    q0=Diagonal(Matrix([1im*Kz0 0I;0I 1im*Kz0]))

    #Free space, so W is identity
    try 
        Q0/q0
    catch SingularException
        throw(error("A singular exception occured, since q0 can not be inverted as it contains elements that equal zero. This is likely due to the z-component of the momentum being zero for some reciprocal lattice vectors. A likely solution is to slightly vary the size of the unit cell (e.g. 700.000001 nm instead of 700 nm)."))
    end
    V0=Q0/q0
    return V0,Kz0
end
"""
	rcwagrid(Nx,Ny,px,py,θ,α,λ,sup)
Create a reciprocal space grid for RCWA simulation
# Arguments
* `nx` : Maximum order in x
* `ny` : Maximum order in y
* `px` : Structure pitch in x
* `py` : Structure pitch in y
* `θ` : inclination angle of incoming wave
* `α` : azimuth angle of incoming wave
* `λ` : wavelength
* `sup` : superstrate material
# Outputs
* `grd`: RCWA grid struct
"""
function rcwagrid(Nx::Int64,Ny::Int64,px::Real,py::Real,θ::Real,α::Real,λ::Real,sup)
    if θ==0
        @warn "The incident angle θ is set to zero, which will result in a division by zero, and thus an invalid result. For vertical incidence simulation, consider a value close to zero, e.g. 1e-5."
    end
    nx,ny,dnx,dny=ngrid(Nx,Ny)
    Kx,Ky,k0=kgrid(nx,ny,px,py,θ,α,λ,sup)
    V0,Kz0=modes_freespace(Kx,Ky)
	return RCWAGrid(Nx,Ny,nx,ny,dnx,dny,px,py,Kx,Ky,k0,V0,Kz0)
end
"""
	rcwasource(grd)
Create a reciprocal space grid for RCWA simulation
# Arguments
* `grd`: RCWA grid struct
# Outputs
* `ψte`: Amplitude vector describing an incoming TE wave
* `ψtm`: Amplitude vector describing an incoming TM wave
"""
function rcwasource(grd,nsup=1) #nsup specification is deprecated here
    #incoming wave vector
	kin=grd.k0
    #the total number of scattering states
    width=(grd.Nx*2+1)*(grd.Ny*2+1)
    #vertical
    normal=[0,0,1]
    #te polarization E-field is perpendicular with z-axis and propagation direction (so, parallel with surface)
    kte=cross(normal,kin)/norm(cross(normal,kin)
)
    #tm polarization E-field is perpendicular with te and propagation direction (so, not necessarily parallel with surface)
    ktm=cross(kin,kte)/norm(cross(kin,kte))
	#initialize array
    ψ0te=zeros(width*2)*1im
	#prefill the 0th scattering order elements for x and y
    ψ0te[convert(Int64,(width+1)/2)]=kte[1]
    ψ0te[convert(Int64,(width+1)/2)+width]=kte[2]
	#same for tm
    ψ0tm=zeros(width*2)*1im
    ψ0tm[convert(Int64,(width+1)/2)]=ktm[1]
    ψ0tm[convert(Int64,(width+1)/2)+width]=ktm[2]
    return ψ0te,ψ0tm
end

