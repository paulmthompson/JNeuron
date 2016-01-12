
#=
Add extracellular electrode to network
=#

function add!(network::Network,extra::Extracellular)

    coeffs=Array(Extra_coeffs,0)
    
    for j=1:length(fieldnames(network.neur))
        for k=1:length(getfield(network.neur,j))
            (mycoeffs,inds)=extracellular(extra,getfield(network.neur,j)[k],0.3)
            push!(coeffs,Extra_coeffs(mycoeffs,inds))
        end
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

function add!(network::Network,intra::Intracellular)

    intra.v=zeros(Float64,length(network.t))

    push!(network.intra,intra)

    nothing
    
end

function run!{T<:NeuronPool}(network::Network{T},init=true)

    #get initial conditions if uninitialized
    if init==true
        for j=1:length(fieldnames(network.neur))
            for k=1:length(getfield(network.neur,j))
                initialcon!(getfield(network.neur,j)[k])
            end
        end
    end

    #set up flags
    if length(network.extra)==0
        extra=false
    else
        extra=true
        cur_inds=get_current(network)
    end

    if length(network.intra)==0
        intra=false
    else
        intra=true
        v_inds=get_voltage(network)
    end

    if length(network.stim)==0
        stim=false
    else
        stim=true
        stim_inds=get_stim(network)
    end
                         
    for i=1:length(network.t)

        #stimulate if appropriate
        if stim==true
            for j=1:length(network.stim)
                add_stim(stim_inds[j],network.stim[j].Is[i])
            end
        end

        for j=1:length(fieldnames(network.neur))
            for k=1:length(getfield(network.neur,j))
                main(getfield(network.neur,j)[k])
            end
        end
 
        if extra==true

            for j=1:length(cur_inds)
                for l=1:length(network.extra)
                    network.extra[l].v[i]+=sum(network.extra[l].coeffs[j].c.*fetch_current(cur_inds[j]))
                end
            end

        end
        
        #get intracellular potentials
        for j=1:length(network.intra)
            network.intra[j].v[i]=fetch_voltage(v_inds[j])
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

fetch_current(myind::RemoteRef)=fetch(myind)

function get_current(network::Network)

    count=0
    for j=1:length(fieldnames(network.neur))
        for k=1:length(getfield(network.neur,j))
            count+=1
        end
    end
    
    cur=Array(SubArray{Float64,1},count)
    count=1
    
    for j=1:length(fieldnames(network.neur))
        for k=1:length(getfield(network.neur,j))
            cur[count]=sub(getfield(network.neur,j)[k].i_vm,network.extra[1].coeffs[count].inds)
            count+=1
        end
    end

    cur
       
end

fetch_current(myind::SubArray{Float64,1})=myind

function get_stim{T <: DArray{Neuron,1}}(network::Network{T})

    myind=Array(RemoteRef,length(network.stim))
    
    #stimulate if appropriate
    for j=1:length(network.stim)
        
        thisneur=network.stim[j].neur
        neurid=findfirst(thisneur.>=network.neur.cuts[1])
        myind[j] = @spawnat network.neur.pids[neurid] sub(network.neur[thisneur].rhs,network.stim[j].node:network.stim[j].node)
        
    end

    myind
    
end

function add_stim(myind::RemoteRef,Is::Float64)
    
    @spawn fetch(myind)[1]+=Is

    nothing
end

function get_stim(network::Network)

    count=0
    for j=1:length(fieldnames(network.neur))
        for k=1:length(getfield(network.neur,j))
            count+=1
        end
    end

    stim=Array(SubArray{Float64,1},count)
    count=1

    myind=[sub(getfield(network.neur,network.stim[i].mtype)[network.stim[i].neur].rhs, network.stim[i].node:network.stim[i].node) for i=1:length(network.stim)]
    
end

add_stim(myind::SubArray,Is::Float64)=myind[1]+=Is

function get_voltage{T <: DArray{Neuron,1}}(network::Network{T})

    myind=Array(RemoteRef,length(network.intra))

    #get voltage at particular location
    for j=1:length(network.intra)
        
        thisneur=network.intra[j].neur
        neurid=findfirst(thisneur.>=network.neur.cuts[1])
        myind[j] = @spawnat network.neur.pids[neurid] sub(network.neur[thisneur].v, network.intra[j].node:network.intra[j].node)
        
    end

    myind
    
end

fetch_voltage(myind::RemoteRef)=fetch(myind[1])

function get_voltage(network::Network)
    myind=[sub(getfield(network.neur,network.intra[i].mtype)[network.intra[i].neur].v, network.intra[i].node:network.intra[i].node) for i=1:length(network.intra)]
end

fetch_voltage(myind::SubArray)=myind[1]

