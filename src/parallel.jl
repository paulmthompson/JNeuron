
function add!(n::NetworkP,extra::Extracellular)

    extra.v=zeros(Float64,length(n.t))
    coeffs=Array(Extra_coeffs,0)
    
    for j=1:length(fieldnames(n.neur))-4
        mycoeffs=add_extra(getfield(n.neur,j),extra)
        append!(coeffs,mycoeffs)
    end

    add_extra_(n,extra,coeffs)
end

function add!(n::NetworkP,intra::Intracellular)
    add_intra_(n,intra)
    (ind,newx)=neuron_ind(intra.neur,n.neur.c[intra.mtype])
    l=add_intra(ind,newx,intra.node,length(n.t),getfield(n.neur,intra.mtype))
    push!(n.neur.i,l)
    nothing
end

function add!(n::NetworkP,s::Stim)   
    add_stim_(n,s)
    (ind,newx)=neuron_ind(s.neur,n.neur.c[s.mtype])
    l=add_stim(ind,newx,s.node,s.Is,getfield(n.neur,s.mtype))
    push!(n.neur.s,l)
end

function run!(n::NetworkP,init=false)

    #get initial conditions if uninitialized
    if init==true
        init!(n)
    end
                         
    for i=1:length(n.t)
        @inbounds main(n.neur,i)             
    end

    for j=1:length(n.extra)
        for k=1:length(fieldnames(n.neur))-4
            for l=1:length(n.neur.c[k])
                n.extra[j].v+=fetch_extra(getfield(n.neur,k),l,j)
            end
        end
    end
    
    for j=1:length(n.intra)
        @inbounds n.intra[j].v=fetch_intra(getfield(n.neur,n.intra[j].mtype),n.intra[j].mtype,n.neur.i[j])   
    end 

    nothing  
end

init!(n::NetworkP)= @inbounds initialcon!(n.neur,0)

main{T<:Neuron}(n::Array{T,1})=(@inbounds for i=1:length(n);main(n[i]);end;nothing)

main{T<:Puddle}(p::Array{T,1},j::Int64)=(@inbounds for i=1:length(p);main(p[i],j);end;nothing)

function main(p::Puddle,i::Int64)

    for j=1:length(p.h.s)
        @inbounds p.n[p.h.s[j].neur].rhs[p.h.s[j].node]+=p.h.s[j].Is[i]
    end
    
    main(p.n)

    for j=1:length(p.h.e)
        for k=1:length(p.n)
            @inbounds p.h.e[j].v[i]+=a_mult_b(p.h.e[j].coeffs[k],p.n[k])
        end
    end

    for j=1:length(p.h.i)
        @inbounds p.h.i[j].v[i]=p.n[p.h.i[j].neur].v[p.h.i[j].node]
    end

    nothing
end

initialcon!{T<:Neuron}(n::Array{T,1})=(@inbounds for i=1:length(n);initialcon!(n[i]);end;nothing)

initialcon!{T<:Puddle}(p::Array{T,1},i::Int64)=(@inbounds for i=1:length(p);initialcon!(p[i],i);end;nothing)

initialcon!(p::Puddle,i::Int64)=initialcon!(p.n)

function neuron_ind(x::Int64,y::Array{UnitRange{Int64},1})

    ind=0
    newx=0
    mysum=0

    for i in y
    
        if (x<=i[end])&&(x>=i[1])
            newx=x-mysum
            ind+=1
            break
        end
    
        mysum+=length(i[1])
        ind+=1
    end

    (ind,newx)

end

function get_dims(dims)

    y=Array(UnitRange{Int64},0)

    for i in dims
        push!(y,i[1])
    end
    y
end

function add_intra(p::Int64,neur::Int64,node::Int64,t::Int64,n::DArray)
    
    remotecall_fetch(((x,neur,node,t)->(push!(localpart(x)[1].h.i,(Intracellular(1,neur,node,zeros(Float64,t))));length(localpart(x)[1].h.i))),p+1,n,neur,node,t)

end

function add_stim(p::Int64,neur::Int64,node::Int64,s::Array{Float64,1},n::DArray)

    remotecall_fetch(((x,neur,node,t)->(push!(localpart(x)[1].h.s,(Stim(t,1,neur,node,0.0,0.0)));length(localpart(x)[1].h.s))),p+1,n,neur,node,s)

end

function fetch_intra(nd::DArray,p::Int64,ind::Int64)
    remotecall_fetch(((n,j)->localpart(n)[1].h.i[j].v),p+1,nd,ind)
end

function add_extra(n::DArray,e::Extracellular)

    mycoeffs=Array(Extra_coeffs,0)
    
    for p in workers()
        append!(mycoeffs,add_extra(n,e,p))
    end

    mycoeffs  
end

function add_extra(n::DArray,e::Extracellular,p::Int64)
    remotecall_fetch(((x,y)->(c=extracellular(localpart(x)[1].n,y);push!(localpart(x)[1].h.e,typeof(e)(e.xyz,c,zeros(Float64,length(e.v))));c)),p,n,e)
end

function fetch_extra(nd::DArray,p::Int64,ind::Int64)
    remotecall_fetch(((n,j)->localpart(n)[1].h.e[j].v),p+1,nd,ind)
end
