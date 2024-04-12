
using LinearAlgebra

export ScatterMatrix,scattermatrix_ref,scattermatrix_tra,scattermatrix_layer,concatenate,scatMatrices
"""
    ScatterMatrix(S11,S12,S21,S22)

Structure to store the scattering matrix of a layer or halfspace
# Attributes
* `S11` : reflection matrix of port 1
* `S12` : transmission matrix port 2 to port 1
* `S21` : transmission matrix port 1 to port 2
* `S22` : reflection matrix of port 2
"""
struct ScatterMatrix
    S11::AbstractArray{<:Number,2}
    S12::AbstractArray{<:Number,2}
    S21::AbstractArray{<:Number,2}
    S22::AbstractArray{<:Number,2}
end
"""
    scatMatrices(m::RCWAModel,g::RCWAGrid,λ)
Computes the scattering matrices of the device (all layers, superstrate, and substrate)
# Arguments
* `m` : RCWA model object 
* `g` : RCWA grid object 
* `λ` : free-space wavelength 
# Outputs
* `s` : array of scattering matrices
"""
function scatMatrices(m::RCWAModel,g::RCWAGrid,λ)
    s=Array{ScatterMatrix,1}(undef,length(m.layers)+2) #preallocate
    s[1]=scattermatrix_ref(halfspace(g.Kx,g.Ky,m.εsup,λ),g.V0) #superstrate
    s[end]=scattermatrix_tra(halfspace(g.Kx,g.Ky,m.εsub,λ),g.V0) #substrate
    for cnt=2:length(m.layers)+1
        s[cnt]=scattermatrix_layer(eigenmodes(g,λ,m.layers[cnt-1]),g.V0) #layers in between
    end
    return s
end
"""
    scattermatrix_ref(sup::Halfspace,V0)
Computes the scattering matrix of the superstrate halfspace
# Arguments
* `sup` : superstrate halfspace eigenmode object 
* `V0` :  Magnetic eigenmodes of free space
# Outputs
* `S` : scattering matrix
"""
function scattermatrix_ref(sup::Halfspace,V0)
    A=I+V0\sup.V #boundary conditions, W=W0=I
    B=I-V0\sup.V
    Ai=I/A #precompute inverse for speed
    S11=-Ai*B
    S12=2*Ai
    S21=.5*(A-B*Ai*B)
    S22=B*Ai
    return ScatterMatrix(S11,S12,S21,S22)
end
"""
    scattermatrix_tra(sub::Halfspace,V0)
Computes the scattering matrix of the substrate halfspace
# Arguments
* `sub` : superstrate halfspace eigenmode object 
* `V0` :  Magnetic eigenmodes of free space
# Outputs
* `S` : scattering matrix
"""
function scattermatrix_tra(sub::Halfspace,V0)
    A=I+V0\sub.V #boundary conditions, W=W0=I
    B=I-V0\sub.V
    Ai=I/A #precompute inverse of A for speed
    S11=B*Ai
    S12=.5*(A-B*Ai*B)
    S21=2*Ai
    S22=-Ai*B
    return ScatterMatrix(S11,S12,S21,S22)
end
"""
    scattermatrix_layer(e::Eigenmodes,V0)
Computes the scattering matrix of a layer
# Arguments
* `e` : layer eigenmode object 
* `V0` :  Magnetic eigenmodes of free space
# Outputs
* `S` : scattering matrix
"""
function scattermatrix_layer(e::Eigenmodes,V0)
    A=e.W\I+e.V\V0 #boundary conditions, W0=I
    B=e.W\I-e.V\V0
    Ai=inv(A) #precompute inverse of A for speed
    common=inv(A-e.X*B*Ai*e.X*B) #precompute this part for speed
    S11=S22=common*(e.X*B*Ai*e.X*A-B)
    S12=S21=common*e.X*(A-B*Ai*B)
    return ScatterMatrix(S11,S12,S21,S22)
end
"""
	concatenate(S11a,S12a,S21a,S22a,S11b,S12b,S21b,S22b)
	concatenate(S1::ScatterMatrix,S2::ScatterMatrix)
	concatenate(Sin::Array{ScatterMatrix,1})

Computes the total scattering matrix for combined layers through concatenation
# Arguments
* `S11a` : S11 component of the first scattering matrix
* `S12a` : S12 component of the first scattering matrix
* `S21a` : S21 component of the first scattering matrix
* `S22a` : S22 component of the first scattering matrix
* `S11b` : S11 component of the second scattering matrix
* `S12b` : S12 component of the second scattering matrix
* `S21b` : S21 component of the second scattering matrix
* `S22b` : S22 component of the second scattering matrix
* `S1` : first scattering matrix
* `S2` : second scattering matrix
* `Sin` : array of scattering matrices (chain concatenation)
# Outputs
* `Sout` : total scattering matrix
"""
function concatenate(S11a,S12a,S21a,S22a,S11b,S12b,S21b,S22b)
    S11=S11a+(S12a/(I-S11b*S22a))*S11b*S21a #direct reflection (S11a) plus infinite internal passes
    S12=(S12a/(I-S11b*S22a))*S12b #infinite internal passes -> geometric series
    S21=(S21b/(I-S22a*S11b))*S21a
    S22=S22b+(S21b/(I-S22a*S11b))*S22a*S12b
    return S11,S12,S21,S22
end
function concatenate(S1::ScatterMatrix,S2::ScatterMatrix)
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
