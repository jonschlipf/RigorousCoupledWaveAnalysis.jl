module srcwa
using LinearAlgebra
export srcwa_reftra,a2p,slicehalf,scatterSource,srcwa_matrices,Srcwa_matrices,srcwa_amplitudes,srcwa_abs
using ..models
using ..materials
using ..common
using ..grids
using ..scatterMatrices

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

function scatterSource(kinc,Nx,Ny)
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
#just slices a vector e in half
function slicehalf(e)
    mylength=convert(Int64,size(e,1)/2)
    return e[1:mylength,:],e[mylength+1:end,:]
end
function e2p(ex,ey,ez,Kz,kz0)
    P=abs.(ex).^2+abs.(ey).^2+abs.(ez).^2
    P=sum(real.(Kz)*P/real(kz0))
    return P
end

function srcwa_reftra(a0,model::Model,grd::Grid,λ)
    ref=halfspace(grd.Kx,grd.Ky,get_permittivity(model.εsup,λ))
    tra=halfspace(grd.Kx,grd.Ky,get_permittivity(model.εsub,λ))
    S=scattermatrix_ref(ref,grd.V0)
    for ct=1:length(model.layers)
        S=concatenate(S,scattermatrix_layer(eigenmodes(grd.dnx,grd.dny,grd.Kx,grd.Ky,grd.k0,λ,model.layers[ct]),grd.V0))
    end
    S=concatenate(S,scattermatrix_tra(tra,grd.V0))
    R=a2p(S.S11*a0,grd.Kx,grd.Ky,ref.Kz,grd.kin[3])
    T=a2p(S.S21*a0,grd.Kx,grd.Ky,tra.Kz,grd.kin[3])
    return R,T
end

function srcwa_amplitudes(a0,grd::Grid,mtr::Array{ScatterMatrix,1})
    a=zeros(length(a0),length(mtr)-1)*1im
    b=zeros(length(a0),length(mtr)-1)*1im
    for ct=1:size(a,2)
        Sbefore=concatenate(mtr[1:ct])
        Safter=concatenate(mtr[ct+1:end])
        a[:,ct]=(I-Sbefore.S22*Safter.S11)\(Sbefore.S21*a0)
        b[:,ct]=(I-Safter.S11*Sbefore.S22)\(Safter.S11*Sbefore.S21*a0)
    end
    return a,b
end

function srcwa_abs(a,b,grd::Grid)
    p=zeros(size(a,2))
    for ct=1:size(a,2)
        ex,ey=a2e2d(a[:,ct]+b[:,ct],I)
        hx,hy=a2e2d(-a[:,ct]+b[:,ct],grd.V0)
        p[ct]=imag(sum(ex.*conj.(hy)-ey.*conj.(hx)))/grd.kin[3]
    end
    return p
end



end
