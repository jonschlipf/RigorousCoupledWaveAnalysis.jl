using RCWA,LinearAlgebra,ProgressMeter

m1=ConstantPerm(1.0)

wls=300:10:800
pout=zeros(length(wls))
prad=zeros(length(wls))

pitch=1200

N=3
py=ones((2N+1)^2)

px=pz=0py

@showprogress for i=1:length(wls)
#i=1
λ=wls[i]

mdl=RCWAModel([],m1,m1)

grd=rcwagrid(mdl,N,N,λ,1E-10,90,pitch,pitch)


py1=(I+exp(-2π*.5*1im*Diagonal(grd.nx))+exp(-2π*.5*1im*Diagonal(grd.ny))+exp(-2π*.7*1im*Diagonal(grd.nx))*exp(-2π*.5*1im*Diagonal(grd.ny)))*py/2


sup=halfspace(grd.Kx,grd.Ky,mdl.εsup,λ)
sub=halfspace(grd.Kx,grd.Ky,mdl.εsub,λ)

Ssup=scattermatrix_ref(sup,grd.V0)
Ssub=scattermatrix_tra(sub,grd.V0)

a0,b0=RCWA.dipoleSource(grd,px,py1,pz,Ssup,Ssub)

up=.5RCWA.a2p(Ssup.S12*b0,I,grd.Kx,grd.Ky,sup.Kz,grd.kin[3])
dn=.5RCWA.a2p(Ssub.S21*a0,I,grd.Kx,grd.Ky,sub.Kz,grd.kin[3])

prad[i]=real.(RCWA.dipoleRad(a0,b0,Ssub,grd,sub,px,py1,pz,m1))
pout[i]=up+dn
end
