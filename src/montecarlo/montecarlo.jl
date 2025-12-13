include("core/abstracttypes.jl")
include("core/alias.jl")
include("core/observers.jl")
include("core/interfaces.jl")
include("core/localupdates.jl")

include("model/ising.jl")
include("model/Pottsmodel.jl")
include("model/XYmodel.jl")

const Ising = IsingModel()
const Potts = PottsModel()
const XY = XYModel()
export Ising, Potts, XY