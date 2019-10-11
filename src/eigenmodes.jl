module  eigenmodes
using LinearAlgebra
using ..materials
using ..models
using ..grid
export Eigenmode,Halfspace
export eigenmode,halfspace,eigmodes
struct Eigenmode
    V::AbstractArray{Complex{Float64},2}
    W::AbstractArray{Complex{Float64},2}
    X::AbstractArray{Complex{Float64},2}
end
struct Halfspace
    Kz::AbstractArray{Complex{Float64},2}
    V::AbstractArray{Complex{Float64},2}
end
function eigmodes(model::Model,grd::Rcwagrid,λ)
    ems=Array{Eigenmode,1}(undef,length(model.layers))
    for ct=1:length(ems)
        ems[ct]=eigenmode(grd.dnx,grd.dny,grd.Kx,grd.Ky,grd.k0,λ,model.layers[ct])
    end
    ref=halfspace(grd.Kx,grd.Ky,get_permittivity(model.εsup,λ))
    tra=halfspace(grd.Kx,grd.Ky,get_permittivity(model.εsub,λ))
    return ref,tra,ems
end
function eigenmode(dnx,dny,Kx,Ky,k0,λ,l::PatternedLayer)
    ε=Kx*0+get_permittivity(l.materials[1],λ)*I
    for ct=1:length(l.geometries)
        ε+=reciprocal(l.geometries[ct],dnx,dny)*(get_permittivity(l.materials[ct+1],λ)-get_permittivity(l.materials[ct],λ))
    end
    η=I/ε
    P=[Kx*η*Ky I-Kx*η*Kx;Ky*η*Ky-I -Ky*η*Kx]
    Q=[Kx*Ky ε-Kx*Kx;Ky*Ky-ε -Ky*Kx]
    #eigenmodes
    ev=eigen(Matrix(P*Q))
    q=Diagonal(sqrt.(Complex.(ev.values)))
    q[real.(q).>0].*=-1
    #W is transform between amplitude vector and E-Field
    W=ev.vectors
    #V is transform between amplitude vector and H-Field
    V=Q*W/Diagonal(q)
    #X the factor applied to the amplitudes when propagatin through the layer
    X=exp(q*k0*l.thickness)
    return Eigenmode(V,W,X)
end

function eigenmode(dnx,dny,Kx,Ky,k0,λ,l::PlainLayer)
    ε=get_permittivity(l.material,λ)
    Kz=sqrt.(Complex.(ε*I-Kx*Kx-Ky*Ky))
    Q=[Kx*Ky ε*I-Kx*Kx;Ky*Ky-ε*I -Ky*Kx]
    q=[1im*Kz zeros(size(Kz));zeros(size(Kz)) 1im*Kz]
    q[real.(q).>0].*=-1
    #W is identity
    V=Q/Diagonal(q)
    W=I+0*V
    X=exp(Matrix(q*k0*l.thickness))
    return Eigenmode(V,W,X)
end

function halfspace(Kx,Ky,ε)
    Kz=sqrt.(Complex.(ε*I-Kx*Kx-Ky*Ky))
    Q=[Kx*Ky ε*I-Kx*Kx;Ky*Ky-ε*I -Ky*Kx]
    q0=[1im*Kz zeros(size(Kz));zeros(size(Kz)) 1im*Kz]
    V=Q/Diagonal(q0)
    return Halfspace(Kz,V)
end

end  # module  eigenmodes
