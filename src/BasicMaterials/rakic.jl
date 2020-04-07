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
    ω=1.2398/λ
    epsilon=1 .-(Ωp^2)/ω/(ω-1im*m.Γ0)
    for o in m.oscillators
        epsilon+=o.f*m.ωp^2/(o.ω^2-ω^2+1im*ω*o.Γ)
    end
    return epsilon
end

Al=Metal(14.98,0.523,0.047,[
    Oscillator(.227,.333,.162),
    Oscillator(.050,.312,1.544),
    Oscillator(.166,1.351,1.188),
    Oscillator(.030,3.382,3.473)])

al_rakic(λ)=LorentzDrude(Al,λ)
#A. D. Rakić, A. B. Djurišic, J. M. Elazar, and M. L. Majewski. Optical properties of metallic films for vertical-cavity optoelectronic devices, Appl. Opt. 37, 5271-5283 (1998)
