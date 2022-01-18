#include("mats.jl")
struct Oscillator
    f::Float64
    Γ::Float64
    ω::Float64
end
struct Metal
    ωp::Float64
    f0::Float64
    Γ0::Float64
    oscillators::Array{Oscillator,1}
end

function LorentzDrude(m::Metal,λ)
    Ωp=sqrt(m.f0)*m.ωp
    ω=1239.8/λ
    epsilon=1 .-(Ωp^2)/ω/(ω+1im*m.Γ0)
    for o in m.oscillators
        epsilon+=o.f*m.ωp^2/(o.ω^2-ω^2-1im*ω*o.Γ)
    end
    return epsilon
end

Al=Metal(14.98,0.523,0.047,[
    Oscillator(.227,.333,.162),
    Oscillator(.050,.312,1.544),
    Oscillator(.166,1.351,1.188),
    Oscillator(.030,3.382,3.473)])
Ag=Metal(9.010,0.845,0.048,[
    Oscillator(0.065,3.886,0.816),
    Oscillator(0.124,0.452,4.481), 
    Oscillator(0.011,0.065,8.185), 
    Oscillator(0.840,0.916,9.083), 
    Oscillator(5.646,2.419,20.290)])
Au=Metal(9.030, 0.760,0.053, [
    Oscillator(0.024,0.241,0.415), 
    Oscillator(0.010,0.345,0.830), 
    Oscillator(0.071,0.870,2.969), 
    Oscillator(0.601,2.494,4.304), 
    Oscillator(4.384,2.214,13.320)])
Ti=LorentzDrudeMetal(7.290,0.148,0.082, [
    Oscillator(0.899,2.276,0.777), 
    Oscillator(0.393,2.518,1.545), 
    Oscillator(0.187,1.663,2.509), 
    Oscillator(0.001,1.762,19.430)])

al_rakic(λ)=LorentzDrude(Al,λ)
ag_rakic(λ)=LorentzDrude(Ag,λ)
au_rakic(λ)=LorentzDrude(Au,λ)
ti_rakic(λ)=LorentzDrude(Ti,λ)
#A. D. Rakić, A. B. Djurišic, J. M. Elazar, and M. L. Majewski. Optical properties of metallic films for vertical-cavity optoelectronic devices, Appl. Opt. 37, 5271-5283 (1998)
