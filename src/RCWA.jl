module RCWA
export Layer,PlainLayer,PatternedLayer,Material,ConstantPerm,Circle,RCWAModel
export Meshgrid,ngrid,kgrid,meshgrid,halfspace,etm_reftra
export ScatterMatrix,scattermatrix_ref,scattermatrix_tra,scattermatrix_layer,concatenate
export srcwa_reftra,scatterSource,srcwa_matrices,rcwagrid,Srcwa_matrices,srcwa_amplitudes,RcwaGrid,srcwa_abs
export Circle,Rectangle,Ellipse,etm1,etmsource,eigenmodes,Eigenmode,scatMatrices,get_permittivity
export srcwa_fields,recip2real,real2recip,InterpolPerm,Custom
export Combination,Rotation,Shift,ModelPerm

export ge_nunley,si_schinke,sio2_malitson,zno_bond

include("Common/Common.jl")
include("SRCWA/SRCWA.jl")
include("ETM/ETM.jl")
include("BasicMaterials/BasicMaterials.jl")


using .Common
using .SRCWA
using .ETM
using .BasicMaterials


end # module
