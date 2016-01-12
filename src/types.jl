
global num_neur = 0
global num_prop = 0

abstract Channel
abstract Prop
abstract Neuron
abstract NeuronPool
abstract Source

immutable Pt3d
    x::Float64
    y::Float64
    z::Float64
    d::Float64
    arc::Float64 #normalized distance from 0 to end
end

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

Section!(s::Array{Section,1},ind::Int64,p::UnitRange{Int64})=(s[ind]=Section(s[ind],p); nothing)

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
            secstack::Array{Section,1}
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

type Network{T <: NeuronPool}
    neur::T 
    t::FloatRange{Float64} 
    extra::Array{Extracellular,1} 
    #extracellular Stimulation
    intra::Array{Intracellular,1} 
    stim::Array{Stim,1}
end

function Network(neuron::Neuron,tstop::Float64; par=false)

    gen_neuronpool([neuron],par)

    pool=NeuronPool0([neuron],par)

    Network(pool,0.0:0.025:tstop,Array(Extracellular,0),Array(Intracellular,0),Array(Stim,0))

end

function Network{T<:Neuron}(neurons::Array{T,1},tstop::Float64; par=false)

    gen_neuronpool(neurons,par)

    pool=NeuronPool0(neurons,par)

    Network(pool,0.0:0.025:tstop,Array(Extracellular,0),Array(Intracellular,0),Array(Stim,0))  

end

#types of neurons found in network
function gen_neuronpool{T<:Neuron}(neur::Array{T,1}, par=false)

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
        myfields=[:($(symbol("N_$i"))::DArray{($(typeinds[i])),1}) for i=1:length(typeinds)]
    end
    
    @eval begin
        type $(symbol("NeuronPool0")) <: NeuronPool
            $(myfields...)
        end

        function NeuronPool0{T<:Neuron}(neur::AbstractArray{T,1},par=false)

            inds=Array(Array{Int64,1},0)
            mtypes=[typeof(neur[i]) for i=1:length(neur)]

            for i=1:length(fieldnames(NeuronPool0))
                push!(inds,find(mtypes.==eltype(fieldtype(NeuronPool0,i))))
            end

            if par==false
                NeuronPool0([typeof(neur[inds[i]][1])[neur[inds[i]]...] for i=1:length(inds)]...)
            else
                NeuronPool0([distribute(typeof(neur[inds[i]][1])[neur[inds[i]]...]) for i=1:length(inds)]...)
            end         
        end 
    end
end
