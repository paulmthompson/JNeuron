
#=
Main data container for JNeuron
=#


export Neuron



function Neuron()
    Neuron(Array(Section,0),Array(Float64,0,0),Array(Float64,0),Array(Float64,0))
end

function add_sec(neuron::Neuron, sec::Section) #Add section to workspace
    push!(neuron.secstack,sec)
end