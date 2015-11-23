
module JNeuron

using ArrayViews, DistributedArrays, Yeppp, PowerSeries

#types
export Neuron, Network, Intracellular, Extracellular, Stim

#functions
export input, instantiate, set_nsegs!, add!, add

include("types.jl")
include("Section.jl")
include("Neuron.jl")
include("import_3d.jl")
include("channels.jl")
include("mod2j.jl")
include("solver.jl")
include("extracellular.jl")
include("network.jl")
include("visualization.jl")
include("speedy.jl")

end
