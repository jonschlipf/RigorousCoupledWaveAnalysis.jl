module srcwa
using LinearAlgebra
export srcwa_reftra,srcwa_source,srcwa_matrices,Srcwa_grid,Srcwa_matrices,srcwa_absorptions
using ..models
using ..materials
using ..grid
using ..scatterMatrices

struct Srcwa_grid
    k0::Float64
    Kx::AbstractArray{Complex{Float64},2}
    Ky::AbstractArray{Complex{Float64},2}
    kin::AbstractArray{Float64,1}
    V0::AbstractArray{Complex{Float64},2}
    Kz0::AbstractArray{Complex{Float64},2}
end

struct Srcwa_matrices
    Kzref::AbstractArray{Complex{Float64},2}
    Kztra::AbstractArray{Complex{Float64},2}
    Slayers::Array{ScatterMatrix,1}
end

function srcwa_source(kinc,Nx,Ny)
    #the total number of scattering states
    width=(Nx*2+1)*(Ny*2+1)
    #vertical
    normal=[0,0,1]
    #te polarization E-field is perpendicular with z-axis and propagation direction (so, parallel with surface)
    kte=cross(normal,kinc)/norm(cross(normal,kinc))
    #tm polarization E-field is perpendicular with te and propagation direction (so, not necessarily parallel with surface)
    ktm=cross(kinc,kte)/norm(cross(kinc,kte))
    esource=zeros(width*2)*1im
    esource[convert(Int64,(width+1)/2)]=kte[1]
    esource[convert(Int64,(width+1)/2)+width]=kte[2]
    a0te=esource
    esource=zeros(width*2)*1im
    esource[convert(Int64,(width+1)/2)]=ktm[1]
    esource[convert(Int64,(width+1)/2)+width]=ktm[2]
    a0tm=esource#/sqrt(epsref)
    return a0te,a0tm
end
function a2p(a,Kx,Ky,Kz,kz0)
    ex,ey,ez=a2e(a,Kx,Ky,Kz)
    return e2p(ex,ey,ez,Kz,kz0)
end
function a2e(a,Kx,Ky,Kz)
    ex,ey=slicehalf(a)
    ez=-Kz\(Kx*ex+Ky*ey)
    return ex,ey,ez
end
function a2e2d(a,W)
    e=W*a
    ex,ey=slicehalf(e)
    return ex,ey
end
#just slices a vector e in half
function slicehalf(e)
    mylength=convert(Int64,length(e)/2)
    return e[1:mylength],e[mylength+1:end]
end
function e2p(ex,ey,ez,Kz,kz0)
    P=abs.(ex).^2+abs.(ey).^2+abs.(ez).^2
    P=sum(real.(Kz)*P/real(kz0))
    return P
end

function srcwa_matrices(model,Nx,Ny,λ,θ,α,ax,ay)
    nx,ny,dnx,dny=ngrid(Nx,Ny)
    k0,Kx,Ky,kin=kgrid(nx,ny,θ,α,λ,ax,ay,get_permittivity(model.εsup,λ))
    V0,Kz0=modes_freespace(Kx,Ky)

    Slayers=Array{ScatterMatrix,1}(undef,length(model.layers)+2)
    for ct=2:length(Slayers)-1
        Slayers[ct]=scattermatrix_layer(dnx,dny,Kx,Ky,k0,λ,model.layers[ct-1],V0)
    end
    Slayers[1],Kzref=scattermatrix_ref(Kx,Ky,get_permittivity(model.εsup,λ),V0)
    Slayers[end],Kztra=scattermatrix_tra(Kx,Ky,get_permittivity(model.εsub,λ),V0)
    return Srcwa_grid(k0,Kx,Ky,kin,V0,Kz0),Srcwa_matrices(Kzref,Kztra,Slayers)
end
function srcwa_reftra(a0,grd::Srcwa_grid,mtr::Srcwa_matrices)
    S=concatenate(mtr.Slayers)
    #a0te,a0tm=srcwa_source(grd.kin,Nx,Ny)
    aRte=S.S11*a0
    R=a2p(S.S11*a0,grd.Kx,grd.Ky,mtr.Kzref,grd.kin[3])
    T=a2p(S.S21*a0,grd.Kx,grd.Ky,mtr.Kztra,grd.kin[3])
    return R,T
end
function srcwa_absorptions(a0,grd::Srcwa_grid,mtr::Srcwa_matrices)
    power=zeros(length(mtr.Slayers)-1)
    for ct=1:length(power)
        Sbefore=concatenate(mtr.Slayers[1:ct])
        Safter=concatenate(mtr.Slayers[ct+1:end])
        a=(I-Sbefore.S22*Safter.S11)\(Sbefore.S21*a0)
        b=(I-Safter.S11*Sbefore.S22)\(Safter.S11*Sbefore.S21*a0)
        ex,ey=a2e2d(a+b,I)
        hx,hy=a2e2d(-a+b,grd.V0)
        power[ct]=imag(sum(ex.*conj.(hy)-ey.*conj.(hx)))/grd.kin[3]
    end
    return power[1:end-1]-power[2:end]
end




end
