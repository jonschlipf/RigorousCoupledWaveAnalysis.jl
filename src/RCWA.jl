module RCWA
export Layer,PlainLayer,PatternedLayer,Material,ConstantPerm,Circle,Model
export Meshgrid,ngrid,kgrid,meshgrid,halfspace,etm_reftra
export ScatterMatrix,scattermatrix_ref,scattermatrix_tra,scattermatrix_layer,concatenate
export srcwa_reftra,scatterSource,srcwa_matrices,rcwagrid,Srcwa_matrices,srcwa_amplitudes,RcwaGrid,srcwa_abs
export Circle,Rectangle,Ellipse,etm1,etmsource,eigenmodes,Eigenmode,scatMatrices,get_permittivity
export srcwa_fields,recip2real,real2recip,InterpolPerm,Custom
export Combination,Rotation,Shift,ModelPerm
include("ft2d.jl")
include("materials.jl")
include("models.jl")
include("grids.jl")
include("common.jl")
include("scatterMatrices.jl")
include("srcwa.jl")
include("etm.jl")


using .models
using .materials
using .grids
using .common
using .scatterMatrices
using .srcwa
using .ft2d
using .etm

greet() = print("Hello World!")

end # module
