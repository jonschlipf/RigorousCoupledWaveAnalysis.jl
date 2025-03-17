using RigorousCoupledWaveAnalysis,LinearAlgebra
using Test
include("analytical.jl")
@testset "VerticalIncidence" begin
n1=rand()*10
n2=rand()*10
mdl=RCWAModel([],ConstantPerm(n1^2),ConstantPerm(n2^2))
grd=rcwagrid(0,0,100,100,1E-3,rand()*360,1000,ConstantPerm(n1^2))
ψte,~=rcwasource(grd,n1)
R,T=srcwa_reftra(ψte,mdl,grd,1000)
R0,~,T0,~=fresnelPower(n1,n2,0)
@test isapprox(R,R0,atol=1E-5)
@test isapprox(T,T0,atol=1E-5)
@test isapprox(R+T,1,atol=1E-5)
end

@testset "GrazingIncidence" begin
n1=rand()*10
n2=n1+rand()*10
θ=1+88rand()
mdl=RCWAModel([],ConstantPerm(n1^2),ConstantPerm(n2^2))
grd=rcwagrid(0,0,100,100,θ,rand()*360,1000,ConstantPerm(n1^2))
ψte,ψtm=rcwasource(grd,n1)
Rte,Tte=srcwa_reftra(ψte,mdl,grd,1000)
Rtm,Ttm=srcwa_reftra(ψtm,mdl,grd,1000)
R0s,R0p,T0s,T0p=fresnelPower(n1,n2,θ)
@test isapprox(Rte,R0s,atol=1E-5)
@test isapprox(Tte,T0s,atol=1E-5)
@test isapprox(Rtm,R0p,atol=1E-5)
@test isapprox(Ttm,T0p,atol=1E-5)
@test isapprox(Rte+Tte,1,atol=1E-5)
@test isapprox(Rtm+Ttm,1,atol=1E-5)
end

@testset "Fabryperot" begin
n1=1+1rand()
n2=n1+rand()+1im*rand() #n2>n1 against TIR
n3=n1+rand() #n3>n1 against TIR
θ=88rand()+1
λ=1000rand()
thickness=100rand()
#θ2=asind.(real.(n1)./real.(n2).*sind.(θ))
#thickness=λ/real(n2)/cosd(θ2)
mdl=RCWAModel([SimpleLayer(thickness,ConstantPerm(n2^2))],ConstantPerm(n1^2),ConstantPerm(n3^2))
grd=rcwagrid(0,0,100,100,θ,rand()*360,λ,ConstantPerm(n1^2))
ψte,ψtm=rcwasource(grd)
Rte,Tte=srcwa_reftra(ψte,mdl,grd,λ)
Rtm,Ttm=srcwa_reftra(ψtm,mdl,grd,λ)
R0s,R0p,T0s,T0p=fabryperot(n1,n2,n3,2π/λ,thickness,θ)
@test isapprox(Rte,R0s,atol=1E-5)
@test isapprox(Rtm,R0p,atol=1E-5)
@test isapprox(Tte,T0s,atol=1E-5)
@test isapprox(Ttm,T0p,atol=1E-5)
end

@testset "unitySource" begin
n1=1+10rand()
λ=1000rand()
θ=88rand()+1
α=360rand()
grd=rcwagrid(0,0,100rand(),100rand(),θ,α,λ,ConstantPerm(n1^2))
ste,stm=rcwasource(grd,(n1))
ref=RigorousCoupledWaveAnalysis.halfspace(grd.Kx,grd.Ky,ConstantPerm(n1^2),λ)
P1=-RigorousCoupledWaveAnalysis.a2p(ste,0ste,ref.V,I,grd.k0[3])
P2=-RigorousCoupledWaveAnalysis.a2p(stm,0ste,ref.V,I,grd.k0[3])
@test P1≈1
@test P2≈1
end

@testset "RT_ScatVSETM_simple" begin
eps1=10rand()+10im*rand()
eps2=10rand()+10im*rand()
eps3=10rand()+10im*rand()
eps4=10rand()+10im*rand()
λ=1000rand()
mdl=RCWAModel([PatternedLayer(100rand(),[ConstantPerm(eps2),ConstantPerm(eps3)],[Circle(rand())])],ConstantPerm(eps1),ConstantPerm(eps4))
grd=rcwagrid(1,1,100rand(),100rand(),88rand()+.1,360rand(),λ,ConstantPerm(eps1))
ste,stm=rcwasource(grd,(√eps1))
Rte1,Tte1=etm_reftra(ste,mdl,grd,λ) 
Rtm1,Ttm1=etm_reftra(stm,mdl,grd,λ) 
Rte2,Tte2=srcwa_reftra(ste,mdl,grd,λ) 
Rtm2,Ttm2=srcwa_reftra(stm,mdl,grd,λ) 
@test isapprox(Rte1,Rte2,atol=1E-5)
@test isapprox(Rtm1,Rtm2,atol=1E-5)
@test isapprox(Tte1,Tte2,atol=1E-5)
@test isapprox(Ttm1,Ttm2,atol=1E-5)
end


@testset "RT_ScatVSETM_conservation_simple" begin
eps1=10rand()
eps2=10rand()
eps3=10rand()
eps4=10rand()
λ=1000rand()
mdl=RCWAModel([PatternedLayer(100rand(),[ConstantPerm(eps2),ConstantPerm(eps3)],[Circle(rand())])],ConstantPerm(eps1),ConstantPerm(eps4))
grd=rcwagrid(1,1,100rand(),100rand(),88rand()+.1,360rand(),λ,ConstantPerm(eps1))
ste,stm=rcwasource(grd,(√eps1))
Rte1,Tte1=etm_reftra(ste,mdl,grd,λ) 
Rtm1,Ttm1=etm_reftra(stm,mdl,grd,λ) 
Rte2,Tte2=srcwa_reftra(ste,mdl,grd,λ) 
Rtm2,Ttm2=srcwa_reftra(stm,mdl,grd,λ) 
@test isapprox(Rte1,Rte2,atol=1E-5)
@test isapprox(Rtm1,Rtm2,atol=1E-5)
@test isapprox(Tte1,Tte2,atol=1E-5)
@test isapprox(Ttm1,Ttm2,atol=1E-5)
@test Rte1+Tte1≈1
@test Rtm1+Ttm1≈1
@test Rte2+Tte2≈1
@test Rtm2+Ttm2≈1
end


@testset "RT_TaylorVSETM_simple" begin
eps1=10rand()+10im*rand()
eps2=10rand()+10im*rand()
eps3=10rand()+10im*rand()
eps4=10rand()+10im*rand()
λ=1000rand()
px=100rand()
py=100rand()
    mdl=RCWAModel([PatternedLayer(minimum([px,py])*rand(),[ConstantPerm(eps2),ConstantPerm(eps3)],[Circle(rand())])],ConstantPerm(eps1),ConstantPerm(eps4)) # thickness > pitch for sufficiently accurate Taylor
grd=rcwagrid(1,1,px,py,88rand()+.1,360rand(),λ,ConstantPerm(eps1))
ste,stm=rcwasource(grd,(√eps1))
Rte1,Tte1=etm_reftra(ste,mdl,grd,λ) 
Rtm1,Ttm1=etm_reftra(stm,mdl,grd,λ) 
Rte2,Tte2=taylor_reftra(ste,mdl,grd,λ) 
Rtm2,Ttm2=taylor_reftra(stm,mdl,grd,λ) 
@test isapprox(Rte1,Rte2,atol=1E-5)
@test isapprox(Rtm1,Rtm2,atol=1E-5)
@test isapprox(Tte1,Tte2,atol=1E-5)
@test isapprox(Ttm1,Ttm2,atol=1E-5)
end

