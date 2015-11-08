

#=
Add neuron to network
=#

function add!(network::Network,neuron::Neuron)
    
    push!(network.neur,neuron)

    nothing
end

#=
Add extracellular electrode to network
=#

function add!(network::Network,extra::Extracellular)

    coeffs=Array(Extra_coeffs,0)
    
    for j=1:length(network.neur)
        push!(coeffs,Extra_coeffs(extracellular(extra,network.neur[j],0.3)))
    end
    
    extra=typeof(extra)(extra.xyz,coeffs,zeros(Float64,length(network.t)))
    
    push!(network.extra,extra)

    nothing
end

#=
Add Intracellular Stimulation to network
=#

function add!(network::Network,stim::Stim)

    myis=zeros(Float64,length(network.t))

    startind=findfirst(network.t.>stim.tstart)
    endind=findfirst(network.t.>stim.tstop)

    myis[startind:endind]=stim.Is[1]

    stim.Is=myis

    push!(network.stim,stim)

    nothing
    
end

function run(network::Network)

    #get initial conditions

    map(initialcon!,network.neur)

    cur_inds=get_current(network.neur)
                             
    for i=1:length(network.t)

        #stimulate if appropriate
        for j=1:length(network.stim)
            network.neur[network.stim[j].neur].rhs[network.stim[j].node]+=network.stim[j].Is[i]
        end

        map(main,network.neur)
        
        for j=1:length(network.neur)

            #get extracellular potential if needed
            for k=1:length(network.extra)
                network.extra[k].v[i]+=sum(network.extra[k].coeffs[j].c.*fetch_current(cur_inds[j]))
            end
            
        end
                                                                                       
        #get intracellular potentials
        for j=1:length(network.intra)
            network.intra[j].v[i]=network.neur[network.intra[j].neur].nodes[network.intra[j].node].v
        end
               
    end

    nothing
        
end

function get_current(neur::DArray{Neuron,1})

    myind=Array(RemoteRef,length(neur))

    count=1
    
    for i=1:length(neur.indexes)
        for j in neur.indexes[i][1]
            myind[count] = @spawnat neur.pids[i] view(neur[j].i_vm,1:neur[j].internal_nodes)
            count+=1
        end

    end

    myind

end

function fetch_current(myind::RemoteRef)
    fetch(myind)
end

function get_current(neur::Array{Neuron,1})
    
    [view(neur[i].i_vm,1:neur[i].internal_nodes) for i=1:length(neur)]
       
end

function fetch_current(myind::ContiguousView)
    myind
end

                   
