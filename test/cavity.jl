using RCWA
n1=1+1rand()+1im*rand()
n2=n1+rand()+1im*rand() #n2>n1 against TIR
n3=n1+rand()+1im*rand() #n3>n1 against TIR
n1=1
n2=2
n3=3
n2=real(n2)
n1=real(n1)
n3=real(n3) 
#n3=n2
θ=88rand()+1
θ=1E-5
θ=45
θ2=asind.(real.(n1)./real.(n2).*sind.(θ))
rs=zeros(length(0:.01:3))
λ=200:1000
for i=1:length(d)
thickness=λ/cosd(θ2)/real(n2)
#thickness=d[i]*n1
#thickness=0
thickness*=d[i]
mdl=RCWAModel([SimpleLayer(thickness,ConstantPerm(n2^2))],ConstantPerm(
n1^2),ConstantPerm(n3^2))
grd=rcwagrid(0,0,100,100,θ,rand()*360,λ)
ψte,ψtm=rcwasource(grd,n1)
Rte,Tte=srcwa_reftra(ψte,mdl,grd,λ)
rs[i]=Rte
end
