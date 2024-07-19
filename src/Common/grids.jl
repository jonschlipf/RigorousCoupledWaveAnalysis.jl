using LinearAlgebra,CUDA
export ngrid,kgrid,rcwagrid,modes_freespace,RCWAGrid,rcwasource,rcwagrid_pml
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
    Nx::Integer
    Ny::Integer
    nx::DenseVector{<:Integer}
    ny::DenseVector{<:Integer}
    dnx::AbstractArray{<:Integer,2}
    dny::AbstractArray{<:Integer,2}
    px::Real
    py::Real
    Kx::Diagonal#{<:Number,DenseVector{<:Number}}
    Ky::Diagonal#{<:Number,DenseVector{<:Number}}
    k0::Vector{<:Number}
    V0::AbstractArray{<:Number,2}
    Kz0::Diagonal#{<:Number,Vector{<:Number}}
end
"""
	ngrid(Nx,Ny,use_cude=false)
Computes the orders n of the plane wave expansion
# Arguments
* `Nx` : Maximum order in x
* `Ny` : Maximum order in y
* `use_cuda` : optional, switch to CUDA GPU solver
# Outputs
* `nx` : Order in x
* `ny` : Order in y
* `dnx` : Differential order in x
* `dny` : Differential order in y
"""
function ngrid(Nx::Integer,Ny::Integer,use_cuda=false)
    #reciprocal space coordinate indices
    nx=vec([r  for r in -Nx:Nx, c in -Ny:Ny])
    ny=vec([c  for r in -Nx:Nx, c in -Ny:Ny])
    #difference matrix for 2dft
    dnx=[nx[a]-nx[b] for a in eachindex(nx),b in eachindex(nx)]
    dny=[ny[a]-ny[b] for a in eachindex(ny),b in eachindex(ny)]
    dnxr=use_cuda ? CuArray(dnx) : dnx
    dnyr=use_cuda ? CuArray(dny) : dny
    nxr=use_cuda ? CuArray(nx) : nx
    nyr=use_cuda ? CuArray(ny) : ny
    return nxr,nyr,dnxr,dnyr
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
function kgrid(nx::AbstractVector{<:Integer},ny::AbstractVector{<:Integer},px::Real,py::Real,θ::Real,α::Real,λ::Real,sup)
    #all k vectors are generally normalized to k0 here
    #The incoming wave, transformed from spherical coordinates to normalized cartesian coordinates, |k0|=1
    k0=[sin(θ*π/180)*cos(α*π/180),sin(θ*π/180)*sin(α*π/180),cos(θ*π/180)]*real(sqrt(get_permittivity(sup,λ)))
    #the spatial vectors are the sum of the incoming vector and the reciprocal lattice vectors
    kx=k0[1].+real(λ)*nx/px;
    ky=k0[2].+real(λ)*ny/py;
    #need matrix for later computations
    Kx=Diagonal(Complex.(kx))
    Ky=Diagonal(Complex.(ky))
    return Kx,Ky,k0
end
function pmlgrid(n,pml_fraction::Real,γ::Number=1/(1-1im))
    f=0.0im*n
    for i in eachindex(f)
        f[i]=-pml_fraction*(-1)^n[i]*((1+γ/4)*sinc(pml_fraction*n[i])+1/2*sinc(n[i]*pml_fraction-1)+.5*sinc(pml_fraction*n[i]+1)-γ/8*sinc(n[i]*pml_fraction-2)-γ/8*sinc(n[i]*pml_fraction+2))
        if n[i]==0
            f[i]+=1
        end
    end
    return f
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
function modes_freespace(Kx::Diagonal{<:Number, <:DenseVector{<:Number}},Ky::Diagonal{<:Number, <:DenseVector{<:Number}})
    IM=Diagonal(Kx*0 .+1)
    #just because |k|=1
    Kz0=sqrt.(Complex.(IM-Kx*Kx-Ky*Ky))
	Kz0[imag.(Kz0).<0].*=-1

    #P0 is identity
    Q0=[Kx*Ky IM-Kx*Kx;Ky*Ky-IM -Ky*Kx]
    #propagation

    q0=Diagonal([1im*Kz0 0I;0I 1im*Kz0])

    #Free space, so W is identity
    V0=Q0/q0
    return V0,Kz0
end
"""
	rcwagrid(Nx,Ny,px,py,θ,α,λ,sup,use_cuda=false)
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
* `use_cuda` : optional, switch to CUDA GPU solver
# Outputs
* `grd`: RCWA grid struct
"""
function rcwagrid(Nx::Integer,Ny::Integer,px::Real,py::Real,θ::Real,α::Real,λ::Real,sup,use_gpu=false)
    if use_gpu&&!CUDA.functional()        
        @warn "CUDA not functional, fallback to CPU."
        use_gpu=false
    end
    nx,ny,dnx,dny=ngrid(Nx,Ny,use_gpu)
    Kx,Ky,k0=kgrid(nx,ny,px,py,θ,α,λ,sup)
    V0,Kz0=modes_freespace(Kx,Ky)
	return RCWAGrid(Nx,Ny,nx,ny,dnx,dny,px,py,Kx,Ky,k0,V0,Kz0)
end
function rcwagrid_pml(Nx::Integer,Ny::Integer,px::Real,py::Real,θ::Real,α::Real,λ::Real,sup,pml_fraction::Real,γ::Number=1/(1-1im),use_gpu=false)
    if use_gpu&&!CUDA.functional()        
        @warn "CUDA not functional, fallback to CPU."
        use_gpu=false
    end
    nx,ny,dnx,dny=ngrid(Nx,Ny,use_gpu)
    Kx,Ky,k0=kgrid(nx,ny,px,py,θ,α,λ,sup)
    Kx=pmlgrid(dnx,pml_fraction,γ)*Kx
    Ky=pmlgrid(dny,pml_fraction,γ)*Ky

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
    if grd.dnx isa CuArray
        ψ0tm=CuArray(ψ0tm)
        ψ0te=CuArray(ψ0te)
    end
    return ψ0te,ψ0tm
end

