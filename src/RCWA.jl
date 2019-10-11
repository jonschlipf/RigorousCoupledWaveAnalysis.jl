module RCWA
export Layer,PlainLayer,PatternedLayer,Material,Constant,Circle,Model
export Meshgrid,ngrid,kgrid,meshgrid
export ScatterMatrix,scattermatrix_ref,scattermatrix_tra,scattermatrix_layer,concatenate
export srcwa_reftra,srcwa_source,srcwa_matrices,Srcwa_grid,Srcwa_matrices,srcwa_amplitudes,srcwa_grid,srcwa_abs
export Circle,Rectangle,Ellipse
include("ft2d.jl")
include("materials.jl")
include("models.jl")
include("grid.jl")
include("scatterMatrices.jl")
include("srcwa.jl")


using .models
using .materials
using .grid
using .scatterMatrices
using .srcwa
using .ft2d

greet() = print("Hello World!")

end # module
