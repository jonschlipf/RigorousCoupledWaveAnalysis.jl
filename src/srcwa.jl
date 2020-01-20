module srcwa
using LinearAlgebra
export srcwa_reftra,a2p,slicehalf,scatterSource,srcwa_matrices,Srcwa_matrices,srcwa_amplitudes,srcwa_abs
export srcwa_fields
using ..models
using ..materials
using ..common
using ..grids
using ..scatterMatrices
using ..ft2d



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


function srcwa_reftra(a0,model::RCWAModel,grd::RcwaGrid,λ)
    ref=halfspace(grd.Kx,grd.Ky,get_permittivity(model.εsup,λ))
    tra=halfspace(grd.Kx,grd.Ky,get_permittivity(model.εsub,λ))
    S=scattermatrix_ref(ref,grd.V0)
    for ct=1:length(model.layers)
        S=concatenate(S,scattermatrix_layer(eigenmodes(grd.dnx,grd.dny,grd.Kx,grd.Ky,grd.k0,λ,model.layers[ct]),grd.V0))
    end
    S=concatenate(S,scattermatrix_tra(tra,grd.V0))
    R=a2p(S.S11*a0,I,grd.Kx,grd.Ky,ref.Kz,grd.kin[3])
    T=a2p(S.S21*a0,I,grd.Kx,grd.Ky,tra.Kz,grd.kin[3])
    return R,T
end

function srcwa_amplitudes(a0,grd::RcwaGrid,mtr::Array{ScatterMatrix,1})
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

function srcwa_abs(a,b,grd::RcwaGrid)
    p=zeros(size(a,2))
    for ct=1:size(a,2)
        ex,ey=a2e2d(a[:,ct]+b[:,ct],I)
        hx,hy=a2e2d(-a[:,ct]+b[:,ct],grd.V0)
        p[ct]=imag(sum(ex.*conj.(hy)-ey.*conj.(hx)))/grd.kin[3]
    end
    return p
end

function srcwa_fields(a,b,em,grd,sz)
    W0=I+0*grd.V0
    x=[r  for r in -sz[1]/2+.5:sz[1]/2-.5, c in -sz[2]/2+.5:sz[2]/2-.5]/sz[1]
    y=[c  for r in -sz[1]/2+.5:sz[1]/2-.5, c in -sz[2]/2+.5:sz[2]/2-.5]/sz[2]
    efield=zeros(size(x,1),size(y,2),sz[3],3)*1im
    hfield=zeros(size(x,1),size(y,2),sz[3],3)*1im

    outside=Matrix([W0 W0;grd.V0 -grd.V0])
    inside=Matrix([em.W+0*em.V em.W+0*em.V;em.V -em.V])
    ain,bout=slicehalf(inside\outside*[a[:,1];b[:,1]])
    aout,bin=slicehalf(inside\outside*[a[:,2];b[:,2]])

    for zind=1:sz[3]
        #propagation of the waves
        a=(em.X^((zind-.5)/sz[3]))*ain
        #b=exp(-Matrix(q)*k0*(zind-1))*bout
        b=(em.X^((sz[3]+.5-zind)/sz[3]))*bin
        #convert amplitude vectors to electric fields

        ex,ey,ez=a2e(a+b,em.W,grd.Kx,grd.Ky,grd.Kz0)
        hx,hy,hz=a2e(-a+b,em.V,grd.Kx,grd.Ky,grd.Kz0)
        #convert from reciprocal lattice vectors to real space distribution
        efield[:,:,zind,1]=recipvec2real(grd.nx,grd.ny,ex,x,y)
        efield[:,:,zind,2]=recipvec2real(grd.nx,grd.ny,ey,x,y)
        efield[:,:,zind,3]=recipvec2real(grd.nx,grd.ny,ez,x,y)

        hfield[:,:,zind,1]=recipvec2real(grd.nx,grd.ny,hx,x,y)
        hfield[:,:,zind,2]=recipvec2real(grd.nx,grd.ny,hy,x,y)
        hfield[:,:,zind,3]=recipvec2real(grd.nx,grd.ny,hz,x,y)
    end
    return efield,hfield
end


end
