module scatterMatrices
using LinearAlgebra
using ..materials
using ..models
export ScatterMatrix,scattermatrix_ref,scattermatrix_tra,scattermatrix_layer,concatenate,modes_freespace
struct ScatterMatrix
    S11::Array{Complex{Float64},2}
    S12::Array{Complex{Float64},2}
    S21::Array{Complex{Float64},2}
    S22::Array{Complex{Float64},2}
end
"""
    modes_freespace(Kx,Ky)
    Computes the eigenmodes of propagation through free space, for normalization
    Kx: x-axis component of the propagation vector
    Ky: y-axis component of the propagation vector
    returns
    V0: coordinate transform between free space eigenmode amplitude and magnetic field
    Kz0: z-axis component of the propagation vector in free space
"""
function modes_freespace(Kx,Ky)
    #just because |k|=1
    Kz0=sqrt.(Complex.(I-Kx*Kx-Ky*Ky))
    #P0 is identity
    Q0=[Kx*Ky I-Kx*Kx;Ky*Ky-I -Ky*Kx]
    #propagation
    q0=1im*Kz0
    q0=[q0 q0*0;0*q0 q0]
    #Free space, so W is identity
    #W0=I+0*Q0
    V0=Q0/Diagonal(q0)
    return V0,Kz0
end
function scattermatrix_ref(Kx,Ky,ε,V0)
    Kz=sqrt.(Complex.(ε*I-Kx*Kx-Ky*Ky))
    Q=[Kx*Ky ε*I-Kx*Kx;Ky*Ky-ε*I -Ky*Kx]
    q0=[1im*Kz zeros(size(Kz));zeros(size(Kz)) 1im*Kz]
    V=Q/Diagonal(q0)
    A=I+V0\Matrix(V)
    B=I-V0\Matrix(V)
    Ai=I/A
    S11=-Ai*B
    S12=2*Ai
    S21=.5*(A-B*Ai*B)
    S22=B*Ai
    return ScatterMatrix(S11,S12,S21,S22),Kz
end
function scattermatrix_tra(Kx,Ky,ε,V0)
    Kz=sqrt.(Complex.(ε*I-Kx*Kx-Ky*Ky))
    Q=[Kx*Ky ε*I-Kx*Kx;Ky*Ky-ε*I -Ky*Kx]
    q0=[1im*Kz zeros(size(Kz));zeros(size(Kz)) 1im*Kz]
    V=Q/Diagonal(q0)
    A=I+V0\Matrix(V)
    B=I-V0\Matrix(V)
    Ai=I/A
    S11=B*Ai
    S12=.5*(A-B*Ai*B)
    S21=2*Ai
    S22=-Ai*B
    return ScatterMatrix(S11,S12,S21,S22),Kz
end
function scattermatrix_layer(dnx,dny,Kx,Ky,k0,λ,l::PatternedLayer,V0)
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
    A=Matrix(W)\I+(Matrix(V)\I)*V0
    B=Matrix(W)\I-(Matrix(V)\I)*V0
    Ai=I/A
    S11=S22=(A-X*B*Ai*X*B)\(X*B*Ai*X*A-B)
    S12=S21=(A-X*B*Ai*X*B)\X*(A-B*Ai*B)
    return ScatterMatrix(S11,S12,S21,S22)
end

function scattermatrix_layer(Kx,Ky,k0,λ,l::PlainLayer,V0)
    ε=get_permittivity(l.material,λ)
    Kz=sqrt.(Complex.(ε*I-Kx*Kx-Ky*Ky))
    Q=[Kx*Ky ε*I-Kx*Kx;Ky*Ky-ε*I -Ky*Kx]
    q=[1im*Kz zeros(size(Kz));zeros(size(Kz)) 1im*Kz]
    q[real.(q).>0].*=-1
    #W is identity
    V=Q/Diagonal(q)
    X=exp(Matrix(q*k0*l.thickness))
    A=I+(Matrix(V)\I)*V0
    B=I-(Matrix(V)\I)*V0
    Ai=I/A
    S11=S22=(A-X*B*Ai*X*B)\(X*B*Ai*X*A-B)
    S12=S21=(A-X*B*Ai*X*B)\X*(A-B*Ai*B)
    return ScatterMatrix(S11,S12,S21,S22)
end

function concatenate(S11a,S12a,S21a,S22a,S11b,S12b,S21b,S22b)
    S11=S11a+(S12a/(I-S11b*S22a))*S11b*S21a
    S12=(S12a/(I-S11b*S22a))*S12b
    S21=(S21b/(I-S22a*S11b))*S21a
    S22=S22b+(S21b/(I-S22a*S11b))*S22a*S12b
    return S11,S12,S21,S22
end
function concatenate(S1,S2)
    S11,S12,S21,S22=concatenate(S1.S11,S1.S12,S1.S21,S1.S22,S2.S11,S2.S12,S2.S21,S2.S22)
    return ScatterMatrix(S11,S12,S21,S22)
end
function concatenate(Sin::Array{ScatterMatrix,1})
    Sout=Sin[1]
    for i=2:length(Sin)
        Sout=concatenate(Sout,Sin[i])
    end
    return Sout
end
end
