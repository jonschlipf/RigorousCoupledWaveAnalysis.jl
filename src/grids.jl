module grids
using LinearAlgebra
using ..models
using ..materials
export Meshgrid,ngrid,kgrid,meshgrid,Rcwagrid,rcwagrid,modes_freespace,grid,RcwaGrid

struct Meshgrid
    x::Array{Float64,2}
    y::Array{Float64,2}
end

struct RcwaGrid
    dnx::Array{Float64,2}
    dny::Array{Float64,2}
    k0::Float64
    Kx::Array{Complex{Float64},2}
    Ky::Array{Complex{Float64},2}
    kin::Array{Float64,1}
    V0::Array{Complex{Float64},2}
    Kz0::Array{Complex{Float64},2}
    nx::Array{Float64,1}
    ny::Array{Float64,1}
end

function ngrid(Nx,Ny)
    #reciprocal space coordinate indices
    ny=[c  for r in -Nx:Nx, c in -Ny:Ny]
    nx=[r  for r in -Nx:Nx, c in -Ny:Ny]
    nx=vec(nx)
    ny=vec(ny)
    #difference matrix for 2dft
    dnx=[nx[a]-nx[b] for a in 1:length(nx),b in 1:length(nx)]
    dny=[ny[a]-ny[b] for a in 1:length(ny),b in 1:length(ny)]
    return nx,ny,dnx,dny
end

function kgrid(nx,ny,θ,α,λ,ax,ay,epsilon_ref)
    #straightforward
    #all k vectors are generally normalized to k0 here
    k0=2*π/real(λ)
    #The incoming wave, transformed from spherical coordinates to normalized cartesian coordinates, |kin|=1
    kin=[sin(θ*π/180)*cos(α*π/180),sin(θ*π/180)*sin(α*π/180),cos(θ*π/180)]*real(sqrt(epsilon_ref))
    #the spatial vectors are the sum of the incoming vector and the reciprocal lattice vectors
    kx=kin[1].+(2*π/ax)*nx/k0;
    ky=kin[2].+(2*π/ay)*ny/k0;
    #need matrix for later computations
    Kx=Diagonal(kx)
    Ky=Diagonal(ky)
    return k0,Kx,Ky,kin
end

function meshgrid(acc)
    #creating a square meshgrid
    x=[r  for r in -acc/2+.5:acc/2-.5, c in -acc/2+.5:acc/2-.5]/acc
    y=[c  for r in -acc/2+.5:acc/2-.5, c in -acc/2+.5:acc/2-.5]/acc
    return Meshgrid(x,y)
end

function rcwagrid(model::RCWAModel,Nx,Ny,λ,θ,α,ax,ay)
    nx,ny,dnx,dny=ngrid(Nx,Ny)
    k0,Kx,Ky,kin=kgrid(nx,ny,θ,α,λ,ax,ay,get_permittivity(model.εsup,λ))
    V0,Kz0=modes_freespace(Kx,Ky)
    return RcwaGrid(dnx,dny,k0,Kx,Ky,kin,V0,Kz0,nx,ny)
end

"""
    modes_freespace(Kx,Ky)
    Computes the eigenmodes of propagation through free space, for normalization
    Kx: x-axis component of the propagation vector
    Ky: y-axis component of the propagation vector
    returns
    V0: coordinate transform between free space eigenmode amplitude and magnetic field
    Kz0: z-axis component of the propagation vector in free space
"""
function modes_freespace(Kx,Ky)
    #just because |k|=1
    Kz0=sqrt.(Complex.(I-Kx*Kx-Ky*Ky))
    #P0 is identity
    Q0=[Kx*Ky I-Kx*Kx;Ky*Ky-I -Ky*Kx]
    #propagation
    q0=1im*Kz0
    q0=[q0 q0*0;0*q0 q0]
    #Free space, so W is identity
    #W0=I+0*Q0
    V0=Q0/Diagonal(q0)
    return V0,Kz0
end

end
