module RCWA
export Layer,PlainLayer,PatternedLayer,Material,Constant,Circle,Model
export Meshgrid,ngrid,kgrid,meshgrid
export ScatterMatrix,scattermatrix_ref,scattermatrix_tra,scattermatrix_layer,concatenate
export srcwa_reftra,scatterSource,srcwa_matrices,rcwagrid,Srcwa_matrices,srcwa_amplitudes,rcwagrid,srcwa_abs
export Circle,Rectangle,Ellipse,etm1,etmsource,eigenmode,eigmodes,Eigenmode
include("ft2d.jl")
include("materials.jl")
include("models.jl")
include("grid.jl")
include("eigenmodes.jl")
include("scatterMatrices.jl")
include("srcwa.jl")
include("etm.jl")


using .models
using .materials
using .grid
using .eigenmodes
using .scatterMatrices
using .srcwa
using .ft2d
using .etm

greet() = print("Hello World!")

end # module
