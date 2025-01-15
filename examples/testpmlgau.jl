
using RigorousCoupledWaveAnalysis,FFTW
Si=InterpolPerm(RigorousCoupledWaveAnalysis.si_schinke)
air=ConstantPerm(1.0)
Si=ConstantPerm(2.0+1e-4im)
Si=ConstantPerm(9.0)


function tilted_gaussian_distr(θ,λ,w0,Ex,Ey,dx,dy,px,py)
    k=2π/λ
    xa=-px/2:dx:px/2-dx
    ya=-py/2:dy:py/2-dy
    x=[x for x in xa,y in ya]
    y=[y for x in xa,y in ya]

    xg=x*cosd(θ)
    zg=x*sind(θ).+1e-7
    rsq=y.^2+xg.^2

    zR=π*w0*w0/λ
    R=zg.*(1 .+(zR./zg).^2)
    psi=atan.(zg/zR)
    w=w0*sqrt.(1 .+(zg./zR).^2)

    E=w0./w.*exp.(-rsq.*w.^-2).*exp.(-1im*(k*zg+k*rsq./(2R)-psi))
    return E
end

function tilted_gaussian(θ,λ,w0,Ex,Ey,dx,dy,px,py,grd)
    # there is a mismatch when comparing with recipvec2real
    # problematic?
    E=tilted_gaussian_distr(θ,λ,w0,Ex,Ey,dx,dy,px,py)
    xa=-px/2:dx:px/2-dx
    ya=-py/2:dy:py/2-dy
    Nx=Int(floor(length(xa)/2))
    Ny=Int(floor(length(ya)/2))
    nx=[r for r in -Nx:Nx-1, c in -Ny:Ny-1]
    ny=[c for r in -Nx:Nx-1, c in -Ny:Ny-1]
    x=[x for x in xa,y in ya]
    y=[y for x in xa,y in ya]

    Eshift=E.*exp.(sind(θ)*1im*x*2π/λ)
    Er=fftshift(fft(ifftshift(Eshift)))
    Er2=.0*grd.nx*1im
        nx.*Er
    for index in eachindex(nx)
        Er2[(grd.nx.==nx[index]).&(grd.ny.==ny[index])].=Er[index]
    end
    return cat(Er2*Ex,Er2*Ey,dims=1)
end
function tilted_gaussian_distr_2D(θ,λ,w0,Ex,Ey,dx,dy,px,py)
    k=2π/λ
    xa=-px/2:dx:px/2-dx
    ya=-py/2:dy:py/2-dy
    x=[x for x in xa,y in ya]
    y=[y for x in xa,y in ya]

    xg=x*cosd(θ)
    zg=x*sind(θ).+1e-7
    rsq=xg.^2

    zR=π*w0*w0/λ
    R=zg.*(1 .+(zR./zg).^2)
    psi=atan.(zg/zR)
    w=w0*sqrt.(1 .+(zg./zR).^2)

    E=w0./w.*exp.(-rsq.*w.^-2).*exp.(-1im*(k*zg+k*rsq./(2R)-psi))
    return E
end

function tilted_gaussian_2D(θ,λ,w0,Ex,Ey,dx,dy,px,py,grd)
    # there is a mismatch when comparing with recipvec2real
    # problematic?
    E=tilted_gaussian_distr_2D(θ,λ,w0,Ex,Ey,dx,dy,px,py)
    xa=-px/2:dx:px/2-dx
    ya=-py/2:dy:py/2-dy
    Nx=Int(floor(length(xa)/2))
    Ny=Int(floor(length(ya)/2))
    nx=[r for r in -Nx:Nx-1, c in -Ny:Ny-1]
    ny=[c for r in -Nx:Nx-1, c in -Ny:Ny-1]
    x=[x for x in xa,y in ya]
    y=[y for x in xa,y in ya]

    Eshift=E.*exp.(sind(θ)*1im*x*2π/λ)
    Er=fftshift(fft(ifftshift(Eshift)))
    Er2=.0*grd.nx*1im
        nx.*Er
    for index in eachindex(nx)
        Er2[(grd.nx.==nx[index]).&(grd.ny.==ny[index])].=Er[index]
    end
    return cat(Er2*Ex,Er2*Ey,dims=1)
end

N=10
θ=1e-5+50
px=10000000
py=px
dair=10000*3
hlay=5000000
hlay=px
mdl=RCWAModel([PatternedLayer(hlay/3,[air,Si],[Rectangle(1e-5,1e-5)]),PatternedLayer(hlay/3,[Si,Si],[Rectangle(.5,.5)]),PatternedLayer(hlay/3,[air,Si],[Rectangle(1e-5,1e-5)])],air,air)
mdl=RCWAModel([AnisotropicLayer(hlay/3,air),AnisotropicLayer(hlay/3,Si),AnisotropicLayer(hlay/3,air)],air,air)
grd=rcwagrid_pml(N,N,px,py,θ,0,1001,mdl.εsup,.9px,.9py,5/(1-1im))
grd2=rcwagrid(N,N,px,py,θ,0,1001,mdl.εsup)
ste,stm=rcwasource(grd)
sgau=tilted_gaussian(θ,1001,40000,1,0,10000,10000,px,py,grd)
sgau=stm

xypoints=[100,100]
zpoints=0:hlay/100:hlay
E,H=RigorousCoupledWaveAnalysis.etm_getfields_stack(mdl,grd,xypoints,zpoints,1001.0,sgau)
function poynting(E,H)
    return (E[:,:,:,2].*conj.(H[:,:,:,3])-E[:,:,:,3].*conj.(H[:,:,:,2])),
    (E[:,:,:,3].*conj.(H[:,:,:,1])-E[:,:,:,1].*conj.(H[:,:,:,3])),
    (E[:,:,:,1].*conj.(H[:,:,:,2])-E[:,:,:,2].*conj.(H[:,:,:,1]))
    end
S=poynting(E,H)
using Plots


n=grd.dnx
    nsec=grd.dny
    sz=.6px
    pitch=px
    f=0.0im*n
    f2=0.0im*n
    f3=0.0im*n
γ=1/(1-1im)
    for i in eachindex(f)
        q=pitch-sz
    f[i]=sinc(n[i])*px-q*(-1)^n[i]*sinc(n[i]*q/pitch)*(.5+γ/8)+ # DC
            -.25*q*(-1)^n[i]*sinc(n[i]*q/pitch-1)+
            -.25*q*(-1)^n[i]*sinc(n[i]*q/pitch+1)+ #cos squared
            γ/16*q*(-1)^n[i]*sinc(n[i]*q/pitch-2)+
            γ/16*q*(-1)^n[i]*sinc(n[i]*q/pitch+2) #cos squared times sin squared
        f2[i]=-(-1)^n[i]*q*sinc(n[i]*q/pitch)+ # part x<e/2
            (pitch*sinc(n[i])+q*(-1)^n[i]*sinc(n[i]*q/pitch))*(.5+γ/8)+ # DC
            -.25*q*(-1)^n[i]*sinc(n[i]*q/pitch-1)+
            -.25*q*(-1)^n[i]*sinc(n[i]*q/pitch+1)+ #cos squared
            γ/16*q*(-1)^n[i]*sinc(n[i]*q/pitch-2)+
            γ/16*q*(-1)^n[i]*sinc(n[i]*q/pitch+2) #cos squared times sin squared
    f3[i]=sinc(nsec[i])*px-q*(-1)^nsec[i]*sinc(nsec[i]*q/pitch)*(.5+γ/8)+ # DC
            -.25*q*(-1)^nsec[i]*sinc(nsec[i]*q/pitch-1)+
            -.25*q*(-1)^nsec[i]*sinc(nsec[i]*q/pitch+1)+ #cos squared
            γ/16*q*(-1)^nsec[i]*sinc(nsec[i]*q/pitch-2)+
            γ/16*q*(-1)^nsec[i]*sinc(nsec[i]*q/pitch+2) #cos squared times sin squared
end
        f=f.*sinc.(nsec)
        f2=f2.*sinc.(nsec)
        f3=f3.*sinc.(n)

λ=1001
source=sgau
em=RigorousCoupledWaveAnalysis.eigenmodes(grd,λ,mdl.layers)
ref=RigorousCoupledWaveAnalysis.halfspace(grd.Kx,grd.Ky,mdl.εsup,λ)
tra=RigorousCoupledWaveAnalysis.halfspace(grd.Kx,grd.Ky,mdl.εsub,λ)
b0,a0,b,a=RigorousCoupledWaveAnalysis.etm_propagate(ref,tra,em,source,grd)
ai=exp( em[1].q*2π/λ*px/5)*a[1]
bi=exp(-em[1].q*2π/λ*px/5)*b[1]
