

function add!(n::NetworkS,extra::Extracellular)

    extra.v=zeros(Float64,length(n.t))
    coeffs=Array(Extra_coeffs,0)
    
    for j=1:length(fieldnames(n.neur)),k in getfield(n.neur,j)
        mycoeffs=extracellular(extra,k,0.3)
        push!(coeffs,mycoeffs)
    end

    add_extra_(n,extra,coeffs)
end

add_extra_(n::Network,e::Extracellular,c::Array{Extra_coeffs,1})=(e=typeof(e)(e.xyz,c,e.v);push!(n.extra,e))

add!(n::NetworkS,s::Stim)=add_stim_(n,s)

function add_stim_(n::Network,s::Stim)
    myis=zeros(Float64,length(n.t))

    startind=findfirst(n.t.>s.tstart)
    endind=findfirst(n.t.>s.tstop)

    myis[startind:endind]=s.Is[1]

    s.Is=myis

    push!(n.stim,s)
end

add!(n::NetworkS,intra::Intracellular)=add_intra_(n,intra)

add_intra_(n::Network,i::Intracellular)=(i.v=zeros(Float64,length(n.t));push!(n.intra,i))

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

function a_mult_b(x::Extra_coeffs,n::Neuron)
    c=0.0
    for i=1:length(x.c)
        @inbounds c+=x.c[i]*n.i_vm[x.inds[i]]
    end
    c
end
