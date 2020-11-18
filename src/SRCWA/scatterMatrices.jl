
using LinearAlgebra

export ScatterMatrix,scattermatrix_ref,scattermatrix_tra,scattermatrix_layer,concatenate,scatMatrices
struct ScatterMatrix
    S11::Array{Complex{Float64},2}
    S12::Array{Complex{Float64},2}
    S21::Array{Complex{Float64},2}
    S22::Array{Complex{Float64},2}
end
function scatMatrices(m::RCWAModel,g::RcwaGrid,λ)
    s=Array{ScatterMatrix,1}(undef,length(m.layers)+2)
    s[1]=scattermatrix_ref(halfspace(g.Kx,g.Ky,m.εsup,λ),g.V0)
    s[end]=scattermatrix_tra(halfspace(g.Kx,g.Ky,m.εsub,λ),g.V0)
    for cnt=2:length(m.layers)+1
        s[cnt]=scattermatrix_layer(eigenmodes(g,λ,m.layers[cnt-1]),g.V0)
    end
    return s
end

function scattermatrix_ref(r::Halfspace,V0)
    A=I+V0\Matrix(r.V)
    B=I-V0\Matrix(r.V)
    Ai=I/A
    S11=-Ai*B
    S12=2*Ai
    S21=.5*(A-B*Ai*B)
    S22=B*Ai
    return ScatterMatrix(S11,S12,S21,S22)
end
function scattermatrix_tra(t::Halfspace,V0)
    A=I+V0\Matrix(t.V)
    B=I-V0\Matrix(t.V)
    Ai=I/A
    S11=B*Ai
    S12=.5*(A-B*Ai*B)
    S21=2*Ai
    S22=-Ai*B
    return ScatterMatrix(S11,S12,S21,S22)
end
function scattermatrix_layer(e::Eigenmodes,V0)
    A=e.W\I+(e.V\I)*V0
    B=e.W\I-(e.V\I)*V0
    Ai=I/A
    common=(A-e.X*B*Ai*e.X*B)\I
    S11=S22=common*(e.X*B*Ai*e.X*A-B)
    S12=S21=common*(A-e.X*B*Ai*e.X*B)\e.X*(A-B*Ai*B)
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
