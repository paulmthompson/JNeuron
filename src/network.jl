

function add!(n::NetworkS,extra::Extracellular)

    coeffs=Array(Extra_coeffs,0)
    
    for j=1:length(fieldnames(n.neur)),k in getfield(n.neur,j)
        (mycoeffs,inds)=extracellular(extra,k,0.3)
        push!(coeffs,Extra_coeffs(mycoeffs,inds))
    end
    
    e=typeof(extra)(extra.xyz,coeffs,zeros(Float64,length(n.t)))
    push!(n.extra,e)

end

function add!(network::NetworkS,stim::Stim)

    myis=zeros(Float64,length(network.t))

    startind=findfirst(network.t.>stim.tstart)
    endind=findfirst(network.t.>stim.tstop)

    myis[startind:endind]=stim.Is[1]

    stim.Is=myis

    push!(network.stim,stim)

    nothing
    
end

add!(n::NetworkS,intra::Intracellular)=(intra.v=zeros(Float64,length(n.t));push!(n.intra,intra))

function run!(n::NetworkS,init=false)

    #get initial conditions if uninitialized
    if init==true
        init!(n)
    end
                         
    for i=1:length(n.t)

        #stimulate if appropriate
        #=
        if stim==true
            for j=1:length(network.stim)
                add_stim(stim_inds[j],network.stim[j].Is[i])
            end
        end
        =#

        if n.helper.flags[1]
            count=1
            for j=1 : @num_fields(n.neur)
                for k=1:length(getfield(n.neur,j))
                    @inbounds main(getfield(n.neur,j)[k])
                    for l=1:length(n.extra)
                        n.extra[l].v[i]+=a_mult_b(n.extra[l].coeffs[count],getfield(n.neur,j)[k])
                    end
                    count+=1
                end
            end
        else
            for j=1 : @num_fields(n.neur)
                for k=1:length(getfield(n.neur,j))
                    @inbounds main(getfield(n.neur,j)[k])
                end
            end
        end


        #=
        #get intracellular potentials
        for j=1:length(network.intra)
            network.intra[j].v[i]=fetch_voltage(v_inds[j])
        end
        =#
               
    end

    nothing
        
end

function init!(n::NetworkS)

    for j=1 : @num_fields(n.neur)
        for k=1:length(getfield(n.neur,j))
            initialcon!(getfield(n.neur,j)[k])
        end
    end

    #set up flags
    if length(n.extra)==0
        n.helper.flags[1]=false
    else
        n.helper.flags[1]=true
        #cur_inds=get_current(n)
    end

    if length(n.intra)==0
        n.helper.flags[2]=false
    else
        n.helper.flags[2]=true
        v_inds=get_voltage(n)
    end

    if length(n.stim)==0
        n.helper.flags[3]=false
    else
        n.helper.flags[3]=true
        stim_inds=get_stim(n)
    end
    
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

function get_current(n::NetworkS)

    count=0
    for j=1:length(fieldnames(n.neur))
        for k=1:length(getfield(n.neur,j))
            count+=1
        end
    end
    
    n.helper.extra=Array(SubArray{Float64,1},count)
    count=1
    
    for j=1:length(fieldnames(n.neur))
        for k=1:length(getfield(n.neur,j))
            n.helper.extra[count]=sub(getfield(n.neur,j)[k].i_vm,n.extra[1].coeffs[count].inds)
            count+=1
        end
    end

    nothing
       
end

fetch_current(myind::RemoteRef)=fetch(myind)

fetch_current(myind::SubArray{Float64,1})=myind
#=
function get_stim(network::Network)

    myind=Array(RemoteRef,length(network.stim))
    
    #stimulate if appropriate
    for j=1:length(network.stim)
        
        thisneur=network.stim[j].neur
        neurid=findfirst(thisneur.>=network.neur.cuts[1])
        myind[j] = @spawnat network.neur.pids[neurid] sub(network.neur[thisneur].rhs,network.stim[j].node:network.stim[j].node)
        
    end

    myind
    
end
=#
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

#=
function get_voltage(network::Network)

    myind=Array(RemoteRef,length(network.intra))

    #get voltage at particular location
    for j=1:length(network.intra)
        
        thisneur=network.intra[j].neur
        neurid=findfirst(thisneur.>=network.neur.cuts[1])
        myind[j] = @spawnat network.neur.pids[neurid] sub(network.neur[thisneur].v, network.intra[j].node:network.intra[j].node)
        
    end

    myind
    
end
=#
fetch_voltage(myind::RemoteRef)=fetch(myind[1])

function get_voltage(network::NetworkS)
    myind=[sub(getfield(network.neur,network.intra[i].mtype)[network.intra[i].neur].v, network.intra[i].node:network.intra[i].node) for i=1:length(network.intra)]
end

fetch_voltage(myind::SubArray)=myind[1]

function a_mult_b(x::Extra_coeffs,n::Neuron)
    c=0.0
    for i=1:length(x.c)
        @inbounds c+=x.c[i]*n.v[x.inds[i]]
    end
    c
end
