

function add_neuron!(network::Network,neuron::Neuron)
    
    push!(network.neur,neuron)

    nothing
end

function add_electrode!(network::Network,extra::Extracellular)

    coeffs=Array(Extra_coeffs,0)
    
    for j=1:length(network.neur)
        push!(coeffs,Extra_coeffs(extracellular(network.neur[j],extra.xyz,0.3)))
    end
    
    extra=typeof(extra)(extra.xyz,coeffs,zeros(Float64,length(network.t)))
    
    push!(network.extra,extra)

    nothing
end

function add_stim!(network::Network,Is::Float64,neur::Int64,node::Int64,tstart::Float64,tstop::Float64)
    myis=zeros(Float64,length(network.t))

    startind=findfirst(network.t.>tstart)
    endind=findfirst(network.t.>tstop)

    myis[startind:endind]=Is

    push!(network.stim,Stim(myis,neur,node))

    nothing
    
end

function run(network::Network)

    #get initial conditions
    for i=1:length(network.neur)
        JNeuron.initialcon!(network.neur[i],-65.0,network.t.step/network.t.divisor)
    end
                   
    for i=1:length(network.t)

        #stimulate if appropriate
        for j=1:length(network.stim)
            network.neur[network.stim[j].neur].rhs[network.stim[j].node]+=network.stim[j].Is[i]
        end
        
        for j=1:length(network.neur)

            #run simulation
            main(network.neur[j])

            #get extracellular potential if needed
            for k=1:length(network.extra)
                network.extra[k].v[i]+=sum(network.extra[k].coeffs[j].c.*network.neur[j].i_vm[1:size(network.extra[k].coeffs[j].c,1)])
            end
            
        end
                                                                                       
        #get intracellular potentials
        for j=1:length(network.intra)
            network.intra[j].v[i]=network.neur[network.intra[j].neur].nodes[network.intra[j].node].v
        end
               
    end

    nothing
        
end
