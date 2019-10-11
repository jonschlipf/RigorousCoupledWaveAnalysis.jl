module srcwa
using LinearAlgebra
export srcwa_reftra,a2p,slicehalf,scatterSource,srcwa_matrices,Srcwa_matrices,srcwa_amplitudes,srcwa_abs
using ..models
using ..materials
using ..common
using ..grids
using ..scatterMatrices



function scatterSource(kinc,Nx,Ny)
    #the total number of scattering states
    width=(Nx*2+1)*(Ny*2+1)
    #vertical
    normal=[0,0,1]
    #te polarization E-field is perpendicular with z-axis and propagation direction (so, parallel with surface)
    kte=cross(normal,kinc)/norm(cross(normal,kinc))
    #tm polarization E-field is perpendicular with te and propagation direction (so, not necessarily parallel with surface)
    ktm=cross(kinc,kte)/norm(cross(kinc,kte))
    a0te=zeros(width*2)*1im
    a0te[convert(Int64,(width+1)/2)]=kte[1]
    a0te[convert(Int64,(width+1)/2)+width]=kte[2]

    a0tm=zeros(width*2)*1im
    a0tm[convert(Int64,(width+1)/2)]=ktm[1]
    a0tm[convert(Int64,(width+1)/2)+width]=ktm[2]

    return a0te,a0tm
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
