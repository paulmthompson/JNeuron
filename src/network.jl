

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
end

add!(n::NetworkS,intra::Intracellular)=(intra.v=zeros(Float64,length(n.t));push!(n.intra,intra))

function run!(n::NetworkS,init=false)

    #get initial conditions if uninitialized
    if init==true
        init!(n)
    end
                         
    for i=1:length(n.t)

        for j=1:length(n.stim)
            @inbounds getfield(n.neur,n.stim[j].mtype)[n.stim[j].neur].rhs[n.stim[j].node]+=n.stim[j].Is[i]
        end

        if n.helper.flags[1]
            count=1
            for j=1 : length(fieldnames(n.neur))
                for k=1:length(getfield(n.neur,j))
                    @inbounds main(getfield(n.neur,j)[k])
                    for l=1:length(n.extra)
                        @inbounds n.extra[l].v[i]+=a_mult_b(n.extra[l].coeffs[count],getfield(n.neur,j)[k])
                    end
                    count+=1
                end
            end
        else
            @inbounds main(n.neur)
        end

        for j=1:length(n.intra)
            @inbounds n.intra[j].v[i]=getfield(n.neur,n.intra[j].mtype)[n.intra[j].neur].v[n.intra[j].node]
        end             
    end

    nothing       
end

function run!(n::NetworkP,init=false)

    #get initial conditions if uninitialized
    if init==true
        init!(n)
    end
                         
    for i=1:length(n.t)

        #=
        for j=1:length(n.stim)
            @inbounds getfield(n.neur,n.stim[j].mtype)[n.stim[j].neur].rhs[n.stim[j].node]+=n.stim[j].Is[i]
        end
        =#

        @inbounds main(n.neur)
        #=
        if n.helper.flags[1]
            count=1
            for j=1 : length(fieldnames(n.neur))
                for k=1:length(getfield(n.neur,j))
                    @inbounds main(getfield(n.neur,j)[k])
                    for l=1:length(n.extra)
                        @inbounds n.extra[l].v[i]+=a_mult_b(n.extra[l].coeffs[count],getfield(n.neur,j)[k])
                    end
                    count+=1
                end
            end
        else
            @inbounds main(n.neur)
        end
        =#

        #=
        for j=1:length(n.intra)
            @inbounds n.intra[j].v[i]=getfield(n.neur,n.intra[j].mtype)[n.intra[j].neur].v[n.intra[j].node]
        end   
        =#          
    end

    nothing  
end

function init!(n::NetworkS)

    @inbounds initialcon!(n.neur)

    #set up flags
    if length(n.extra)==0
        n.helper.flags[1]=false
    else
        n.helper.flags[1]=true
    end

    if length(n.intra)==0
        n.helper.flags[2]=false
    else
        n.helper.flags[2]=true
    end

    if length(n.stim)==0
        n.helper.flags[3]=false
    else
        n.helper.flags[3]=true
    end
end

function init!(n::NetworkP)
    
    @inbounds initialcon!(n.neur)

    #set up flags
    if length(n.extra)==0
        n.helper.flags[1]=false
    else
        n.helper.flags[1]=true
    end

    if length(n.intra)==0
        n.helper.flags[2]=false
    else
        n.helper.flags[2]=true
    end

    if length(n.stim)==0
        n.helper.flags[3]=false
    else
        n.helper.flags[3]=true
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

fetch_current(myind::RemoteRef)=fetch(myind)

function add_stim(myind::RemoteRef,Is::Float64)
    
    @spawn fetch(myind)[1]+=Is

    nothing
end

fetch_voltage(myind::RemoteRef)=fetch(myind[1])

function a_mult_b(x::Extra_coeffs,n::Neuron)
    c=0.0
    for i=1:length(x.c)
        @inbounds c+=x.c[i]*n.i_vm[x.inds[i]]
    end
    c
end
