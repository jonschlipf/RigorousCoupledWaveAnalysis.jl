module  common
using LinearAlgebra
using ..materials
using ..models
using ..grids
export Eigenmodes,Halfspace
export eigenmodes,halfspace
export a2e,a2e2d,a2p,e2p,slicehalf
struct Eigenmodes
    V::AbstractArray{Complex{Float64},2}
    W::AbstractArray{Complex{Float64},2}
    X::AbstractArray{Complex{Float64},2}
end
struct Halfspace
    Kz::AbstractArray{Complex{Float64},2}
    V::AbstractArray{Complex{Float64},2}
end
function eigenmodes(dnx,dny,Kx,Ky,k0,λ,l::PatternedLayer)
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
    return Eigenmodes(Matrix(V),Matrix(W),X)
end

function eigenmodes(dnx,dny,Kx,Ky,k0,λ,l::PlainLayer)
    ε=get_permittivity(l.material,λ)
    Kz=sqrt.(Complex.(ε*I-Kx*Kx-Ky*Ky))
    Q=[Kx*Ky ε*I-Kx*Kx;Ky*Ky-ε*I -Ky*Kx]
    q=[1im*Kz zeros(size(Kz));zeros(size(Kz)) 1im*Kz]
    q[real.(q).>0].*=-1
    #W is identity
    V=Q/Diagonal(q)
    W=I+0*V
    X=exp(Matrix(q*k0*l.thickness))
    return Eigenmodes(Matrix(V),Matrix(W),X)
end
function eigenmodes(g::RcwaGrid,λ,l::PlainLayer)
    return eigenmodes(g.dnx,g.dny,g.Kx,g.Ky,g.k0,λ,l)
end
function eigenmodes(g::RcwaGrid,λ,l::PatternedLayer)
    return eigenmodes(g.dnx,g.dny,g.Kx,g.Ky,g.k0,λ,l)
end

function halfspace(Kx,Ky,ε)
    Kz=sqrt.(Complex.(ε*I-Kx*Kx-Ky*Ky))
    Q=[Kx*Ky ε*I-Kx*Kx;Ky*Ky-ε*I -Ky*Kx]
    q0=[1im*Kz zeros(size(Kz));zeros(size(Kz)) 1im*Kz]
    V=Q/Diagonal(q0)
    return Halfspace(Kz,V)
end

function eigenmodes(g::RcwaGrid,λ,l::Array{Layer,1})
    rt=Array{Eigenmodes,1}(undef,length(l))
    for cnt=1:length(l)
        rt[cnt]=eigenmodes(g,λ,l[cnt])
    end
    return rt
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
function a2p(a,W,Kx,Ky,Kz,kz0)
    ex,ey,ez=a2e(a,W,Kx,Ky,Kz)
    return e2p(ex,ey,ez,Kz,kz0)
end
function a2e(a,W,Kx,Ky,Kz)
    e=W*a
    ex,ey=slicehalf(e)
    ez=-Kz\(Kx*ex+Ky*ey)
    return ex,ey,ez
end
function a2e2d(a,W)
    e=W*a
    ex,ey=slicehalf(e)
    return ex,ey
end


end  # module  eigenmodes
