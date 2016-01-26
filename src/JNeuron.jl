
module JNeuron

using ArrayViews, DistributedArrays, PowerSeries, Distributions

#types
export Neuron, Network, Intracellular, Extracellular, Stim

#functions
export input, instantiate, set_nsegs!, add!, add, run!

include("types.jl")
include("mutex.jl")
include("Node.jl")
include("import_3d.jl")
include("channels.jl")
include("mod2j.jl")
include("solver.jl")
include("extracellular.jl")
include("network.jl")
include("visualization.jl")
include("speedy.jl")
include("parallel.jl")

end
