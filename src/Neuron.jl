
#=
Main data container for JNeuron
=#


export Neuron

type Neuron
    secstack::Array{Section,1}
end

function Neuron()
    Neuron(Array(Section,0))
end

function add_sec(neuron::Neuron, sec::Section) #Add section to workspace
    append!(neuron.secstack,sec)
end
