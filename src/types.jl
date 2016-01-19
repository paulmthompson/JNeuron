
global num_neur = 0
global num_prop = 0
global num_pool = 0

abstract Channel
abstract Prop
abstract Neuron
abstract Source
abstract NeuronPool
abstract Network
abstract NetworkS <: Network
abstract NetworkP <: Network
abstract Helper
abstract Puddle

immutable Pt3d
    x::Float64
    y::Float64
    z::Float64
    d::Float64
    arc::Float64 #normalized distance from 0 to end
end

Pt3d(pt::Pt3d,xyz::Array{Float64,1})=Pt3d(xyz[1],xyz[2],xyz[3],pt.d,pt.arc)

add(pt::Pt3d,x::Float64,y::Float64,z::Float64)=Pt3d(pt.x+x,pt.y+y,pt.z+z,pt.d,pt.arc)

immutable SegArea
    l::Float64
    r::Float64
    t::Float64
end

SegArea()=SegArea(50.0,50.0,100.0)

immutable SegRi
    l::Float64
    r::Float64
end

SegRi()=SegRi(0.0,0.0)

immutable Node
    area::SegArea 
    ri::SegRi 
    parent::Int64 
    children::Array{Int64,1}
    internal::Bool
    pt3d::UnitRange{Int64}
end

Node(a::SegArea,r::SegRi,p::Int64,c::Array{Int64},m::UnitRange{Int64})=Node(a,r,p,c,true,m)

Node(p::Int64,c::Array{Int64,1},m::Int64)=Node(SegArea(),SegRi(),p,c,false,m:m)

Node(n::Node)=Node(n.area,n.ri,n.parent,n.children,n.internal,n.pt3d)

add_node(n::Neuron,nd::Node)=push!(n.nodes,nd)

type Prop0<:Prop
end

make_prop(a::Prop0,b::Int64)=Prop0()

function gen_prop_check(prop_array)

    global num_prop::Int64
    
    if method_exists(make_prop,(typeof(prop_array),Int))
    else
        num_prop+=1
        gen_prop(prop_array,num_prop)
    end

    nothing   
end

function gen_prop(a,k::Int64)

    if issubtype(typeof(a),Tuple)
        num_fields=sum([length(fieldnames(a[i]))-1 for i=1:length(a)])
    else    
        num_fields=sum(length(fieldnames(a)))
    end

    @eval begin
        type $(symbol("Prop_$k")) <: Prop
            p::Array{Float64,2}
        end

        make_prop(b::$(typeof(a)),n::Int64)=$(symbol("Prop_$k"))(zeros(Float64,$(num_fields),n))     
    end

    c=make_prop(a,0)

    @eval begin       
        make_prop(a::$(typeof(c)),n::Int64)=$(symbol("Prop_$k"))(zeros(Float64,$(num_fields),n))
    end
    
    nothing  
end

immutable Section
    mtype::Int64 
    pnode::UnitRange{Int64} 
    child::Array{Int64,1}
    parent::Int64
    pt3d::Array{Pt3d,1}
    length::Float64
    prop::DataType
end

Section(s::Section,p::Int64)=Section(s.mtype,s.pnode,s.child,p,s.pt3d,s.length,s.prop)

Section(s::Section,p::DataType)=Section(s.mtype,s.pnode,s.child,s.parent,s.pt3d,s.length,p)

Section(s::Section,p::UnitRange{Int64})=Section(s.mtype,p,s.child,s.parent,s.pt3d,s.length,s.prop)

Section!(s::Array{Section,1},ind::Int64,p::Int64)=(s[ind]=Section(s[ind],p); nothing)

Section!(n::Neuron,ind::Int64,p::Int64)=Section!(n.secs,ind,p)

Section!(s::Array{Section,1},ind::Int64,p::UnitRange{Int64})=(s[ind]=Section(s[ind],p); nothing)

Section!(n::Neuron,ind::Int64,p::UnitRange{Int64})=Section!(n.secs,ind,p)

add_sec(neuron::Neuron, sec::Section)=push!(neuron.secs,sec)

type Extra_coeffs
    c::Array{Float64,1}
    inds::Array{Int64,1}
end

type Point <: Source
end

type Line <: Source
end

type Mixed <: Source
end

type Extracellular{S<:Source}
    xyz::Array{Float64,1}
    coeffs::Array{Extra_coeffs,1}
    v::Array{Float64,1}
end

Extracellular(xyz::Array{Float64,1})=Extracellular(Line(),xyz)

Extracellular(s::Source,xyz)=Extracellular{typeof(s)}(xyz,Array(Extra_coeffs,0),Array(Float64,0))

type Stim
    Is::Array{Float64,1}
    mtype::Int64
    neur::Int64
    node::Int64
    tstart::Float64
    tstop::Float64
end

Stim(Is,mtype,neur,node,tstart,tstop)=Stim([Is],mtype,neur,node,tstart,tstop)

type Intracellular
    mtype::Int64
    neur::Int64
    node::Int64
    v::Array{Float64,1}
end

Intracellular(mtype::Int64,neur::Int64,node::Int64)=Intracellular(mtype,neur,node,Array(Float64,0))

type HelperS <: Helper
    flags::Array{Bool,1}
end

HelperS()=HelperS(falses(4))

type HelperP <: Helper
    flags::Array{Bool,1}
    dims::UnitRange{Int64}
    l::Int64
    i::Array{Intracellular,1}
    s::Array{Stim,1}
end

HelperP()=HelperP(falses(4),0:0,0,Array(Intracellular,0),Array(Stim,0))

HelperP(dims::UnitRange{Int64},l::Int64)=HelperP(falses(4),dims,l,Array(Intracellular,0),Array(Stim,0))

make_neuron()=nothing

function gen_neur_check(myprop)

    global num_neur::Int64
    
    if method_exists(make_neuron,(typeof(myprop),Array{Node,1},Array{Section,1},Array{Array{Int64,1}}))
    else     
        num_neur+=1
        gen_neuron(myprop,num_neur)
    end

    nothing
end

function gen_neuron(prop,k::Int64)

    if issubtype(typeof(prop),Prop)
        myprop=(prop,prop,prop,prop)
    else
        myprop=prop
    end
    
    @eval begin
        type $(symbol("Neuron_$k")) <: Neuron
            soma::$(typeof(myprop[1]))
            axon::$(typeof(myprop[2]))
            dendrite::$(typeof(myprop[3]))
            apical::$(typeof(myprop[4]))
            secs::Array{Section,1}
            v::Array{Float64,1} 
            a::Array{Float64,1}
            b::Array{Float64,1}
            d::Array{Float64,1}
            rhs::Array{Float64,1}
            Ra::Float64
            Cm::Float64
            dt::Float64
            nodes::Array{Node,1}
            i_vm::Array{Float64,1}
            divm::Array{Float64,1}
            diag_old::Array{Float64,1}
            internal_nodes::Array{Array{Int64,1},1}
            par::Array{Int64,1}
            v1::Array{Float64,1}
            i2::Array{Float64,1}
        end

        type $(symbol("Puddle_$k")) <: Puddle
            n::Array{$(symbol("Neuron_$k")),1}
            h::HelperP
        end

        Puddle(n::Array{$(symbol("Neuron_$k")),1})=$(symbol("Puddle_$k"))(n,HelperP())
        Puddle(n::Array{$(symbol("Neuron_$k")),1},h::HelperP)=$(symbol("Puddle_$k"))(n,h)
        Puddle(n::Array{$(symbol("Neuron_$k")),1},dims::UnitRange{Int64},l::Int64)=$(symbol("Puddle_$k"))(n,HelperP(dims,l))

        function make_neuron(prop::$(typeof(prop)),n::Array{Node,1},s::Array{Section,1},inter::Array{Array{Int64,1}})

            len=length(n)

            secs=deepcopy(s)
            nodes=deepcopy(n)
            inter_n=deepcopy(inter)
            
            (soma,axon,dendrite,apical)=make_props(prop,inter_n)
            
            n2 = $(symbol("Neuron_$k"))(soma,axon,dendrite,apical,secs,zeros(Float64,len),zeros(Float64,len),zeros(Float64,len),zeros(Float64,len),zeros(Float64,len),35.4,1.0,.025,nodes,zeros(Float64,len),zeros(Float64,len),zeros(Float64,len),inter_n,zeros(Float64,len),zeros(Float64,len),zeros(Float64,len))

            fillA!(n2)
            
            n2           
        end        
    end   
end

add_inode(n::Neuron,i::Int64)=push!(n.internal_nodes[n.secs[i].mtype],length(n.nodes))

find_parents!(n::Neuron)=(for i=1:length(n.secs), j in n.secs[i].child; Section!(n,j,i); end)

function add(neuron::Neuron,prop_array...)

    num_arg=length(prop_array)

    if num_arg==1

        prop_array=prop_array[1]
        
        gen_prop_check(prop_array)
        
        myprop=make_prop(prop_array,0)
        gen_current(prop_array,myprop)
        
    elseif num_arg==4
        
        for i=1:4
            gen_prop_check(prop_array[i])
        end

        myprop=(make_prop(prop_array[1],0),make_prop(prop_array[2],0),make_prop(prop_array[3],0),make_prop(prop_array[4],0))

        for i=1:4
            gen_current(prop_array[i],myprop[i])
        end

    else
        #error checking
    end

    gen_neur_check(myprop)
     
    n1=make_neuron(myprop,neuron.nodes,neuron.secs,neuron.internal_nodes)
    
end

gen_neuron(Prop0(),0)

function make_props(prop::Prop,inter_n::Array{Array{Int64,1}})

    soma=make_prop(prop,length(inter_n[1]))
    axon=make_prop(prop,length(inter_n[2]))
    dendrite=make_prop(prop,length(inter_n[3]))
    apical=make_prop(prop,length(inter_n[4]))

    (soma,axon,dendrite,apical)  
end

function make_props(prop::Tuple,inter_n::Array{Array{Int64,1}})

    soma=make_prop(prop[1],length(inter_n[1]))
    axon=make_prop(prop[2],length(inter_n[2]))
    dendrite=make_prop(prop[3],length(inter_n[3]))
    apical=make_prop(prop[4],length(inter_n[4]))

    (soma,axon,dendrite,apical)  
end

function gen_pool_check(neur,par,ts)

    global num_pool::Int64

    if par==true
        if method_exists(make_ppool,(typeof(neur),))
        else
            gen_npool(neur,par)
            num_pool+=1
        end
        p=make_ppool(neur)
    else
        if method_exists(make_spool,(typeof(neur),))
        else
            gen_npool(neur,par)
            num_pool+=1
        end
        p=make_spool(neur)
    end

    make_network(p,ts)   
end

make_spool()=nothing
make_ppool()=nothing

findmyid(d)=parse(Int,split(string(d),"_")[end])

function findid(d)
    id=parse(Int,split(string(d),"_")[end])
    symbol("Puddle_$id")
end

function findin(d)
    id=parse(Int,split(string(d),"_")[end])
    "Neuron_$id"
end

function gen_npool(neur, par::Bool)

    global num_pool::Int64

    if issubtype(typeof(neur),Neuron)
        myfields=[:($(symbol("N_1"))::Array{($(typeof(neur))),1})]
    else 
        typeinds=Array(DataType,0)
        mtypes=[typeof(neur[i]) for i=1:length(neur)]
        #Determine types of neurons
        for i=1:length(mtypes)
            if sum(mtypes[i].==typeinds)<1
                push!(typeinds,mtypes[i])
            end
        end

        if par==false
            myfields=[:($(symbol("N_$i"))::Array{($(typeinds[i])),1}) for i=1:length(typeinds)]
        else
            myfields=[:($(symbol("N_$i"))::DArray{($(findid(typeinds[i]))),1}) for i=1:length(typeinds)]
        end
    end
        
    if par==false
      
        @eval begin
            type $(symbol("Pool_$num_pool")) <: NeuronPool
                $(myfields...)
            end

            type $(symbol("Network_$num_pool")) <: NetworkS
                neur::$(symbol("Pool_$num_pool"))
                t::FloatRange{Float64} 
                extra::Array{Extracellular,1} 
                #extracellular Stimulation
                intra::Array{Intracellular,1} 
                stim::Array{Stim,1}
                helper::HelperS
            end
            
            function make_spool(neur::$(typeof(neur)))

                if issubtype(typeof(neur),Neuron)
                    $(symbol("Pool_$num_pool"))([neur])
                else
                    inds=Array(Array{Int64,1},0)
                    mtypes=[typeof(neur[i]) for i=1:length(neur)]

                    for i=1:length(fieldnames($(symbol("Pool_$num_pool"))))
                        push!(inds,find(mtypes.==eltype(fieldtype($(symbol("Pool_$num_pool")),i))))
                    end
                    
                    $(symbol("Pool_$num_pool"))([typeof(neur[inds[i]][1])[neur[inds[i]]...] for i=1:length(inds)]...)
                end
            end

            gen_net_func($(symbol("Pool_$num_pool")),"initialcon!")
            gen_net_func($(symbol("Pool_$num_pool")),"main")

            function make_network(p::$(symbol("Pool_$num_pool")),ts::Float64)
                $(symbol("Network_$num_pool"))(p,0.0:.025:ts,Array(Extracellular,0),Array(Intracellular,0),Array(Stim,0),HelperS())
            end
            
        end
    else     
        @eval begin
            type $(symbol("Pool_$num_pool")) <: NeuronPool
                $(myfields...)
                c::Array{Array{UnitRange{Int64},1},1}
                t::Array{Int64,1}
                i::Array{Int64,1}
                s::Array{Int64,1}
            end

            type $(symbol("Network_$num_pool")) <: NetworkP
                neur::$(symbol("Pool_$num_pool"))
                t::FloatRange{Float64} 
                extra::Array{Extracellular,1} 
                #extracellular Stimulation
                intra::Array{Intracellular,1} 
                stim::Array{Stim,1}
            end

            function make_ppool(neur::$(typeof(neur)))

                inds=Array(Array{Int64,1},0)
                mtypes=[findmyid(typeof(neur[i])) for i=1:length(neur)]

                for i=1:length(fieldnames($(symbol("Pool_$num_pool"))))-4
                    nt=findmyid(eltype(fieldtype(($(symbol("Pool_$num_pool"))),i)))
                    push!(inds,find(mtypes .== nt))
                end

                c=[Array(UnitRange{Int64},0) for i=1:length(inds)]
                t=zeros(Int64,length(inds))
                a=[typeof(neur[inds[i][1]])[neur[inds[i]]...] for i=1:length(inds)]
                
                b=Array(Any,length(inds))

                w=length(workers())
                
                for i=1:length(b)
      
                    if length(a[i])>w
                        dims=DistributedArrays.chunk_idxs(length(a[i]),w)[1]
                        l=length(a[i])
                        t[i]=l
                        c[i]=get_dims(dims)
                        mtype=typeof(Puddle([a[i][1]]))
                        b[i]=distribute(mtype[(k=j[1];Puddle(a[i][k],k,l)) for j in dims])
                    else
                        #b=distribute([Puddle(a[i])])
                    end

                end
                
                $(symbol("Pool_$num_pool"))(b...,c,t,Array(Int64,0),Array(Int64,0))
            end

            gen_net_func_p($(symbol("Pool_$num_pool")),"initialcon!")
            gen_net_func_p($(symbol("Pool_$num_pool")),"main")
            
            function make_network(p::$(symbol("Pool_$num_pool")),ts::Float64)
                $(symbol("Network_$num_pool"))(p,0.0:.025:ts,Array(Extracellular,0),Array(Intracellular,0),Array(Stim,0))
            end
        end        
    end
end

Network(n::Neuron,ts::Float64; par=false)=(gen_pool_check(n,par,ts))

Network(n::Tuple,ts::Float64; par=false)=(gen_pool_check(n,par,ts))

function gen_net_func(n::DataType,func::ASCIIString)
       
    a=length(fieldnames(n))
    
    @eval begin
        function $(symbol("$func"))(n::$(n))
            for j = 1 : $a
                for k=1:length(getfield(n,j))
                    $(symbol("$func"))(getfield(n,j)[k])
                end
            end
        end
    end   
    nothing    
end

function gen_net_func_p(n::DataType,func::ASCIIString)

    a=length(fieldnames(n))-4
    
    @eval begin
        function $(symbol("$func"))(n::$(n),i::Int64)
            for j = 1 : $a
                @sync for p in procs(getfield(n,j))
                    @async remotecall_wait((n,t)->($(symbol("$func"))(localpart(n),t)), p, getfield(n,j),i)
                end
                #map($(symbol("$func")),getfield(n,j))
            end
        end
    end   
    nothing    
end
