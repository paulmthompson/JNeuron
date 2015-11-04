
module JNeuron

using ArrayViews

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

end
