

global num_neur = 0
global num_prop = 0

abstract Channel

abstract Prop

type Prop0<:Prop
end

#associated 3d point
type Pt3d
    x::Float64
    y::Float64
    z::Float64
    d::Float64
    arc::Float64 #normalized distance from 0 to end
end

type Node
    ind::Int64
    area::Array{Float64,1} #surface area of left [1] and right[2] part of segment
    ri::Array{Float64,1}  #internal resistance of left[1] and right[2] part of segment
    parent::Int64 #index in node array of parent
    children::Array{Int64,1} #Node(s) from other sections attached to this one
    internal::Bool
    pt3d::UnitRange{Int64}
    prop::DataType
end

type Section
    refcount::Int64 #ID for section, also the place in secstack array
    mtype::Int64 #Cellbody=1,Axon=2,Dendrite=3,Apical=4
    pnode::Array{Int64,1} #one node at center of each segment
    child::Array{Section,1}
    pt3d::Array{Pt3d,1}
    Ra::Float64 #cytoplasmic resistivity
    length::Float64
end

abstract Neuron

type Neuron0 <: Neuron
    
    soma::Array{Prop0,1}
    axon::Array{Prop0,1}
    dendrite::Array{Prop0,1}
    apical::Array{Prop0,1}
    secstack::Array{Section,1}
    v::Array{Float64,1} #intracellular voltage
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

end

function Neuron0()
    Neuron0(Array(Prop,0),Array(Prop,0),Array(Prop,0),Array(Prop,0),Array(Section,0),zeros(Float64,0),zeros(Float64,0),zeros(Float64,0),zeros(Float64,0),zeros(Float64,0),0.0,0.0,0.025,Array(Node,0),zeros(Float64,0),zeros(Float64,0),zeros(Float64,0),[Array(Int64,0) for i=1:4],Array(Int64,0))
end

abstract NeuronPool

type Extra_coeffs
    c::Array{Float64,1}
    inds::Array{Int64,1}
end

abstract Source

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

#Default Line Sources for extracellular Electrodes
function Extracellular(xyz::Array{Float64,1})
    Extracellular{Line}(xyz,Array(Extra_coeffs,0),Array(Float64,0))
end

function Extracellular(source::Source,xyz::Array{Float64,1})
    Extracellular{typeof(source)}(xyz,Array(Extra_coeffs,0),Array(Float64,0))
end

type Stim
    Is::Array{Float64,1}
    mtype::Int64
    neur::Int64
    node::Int64
    tstart::Float64
    tstop::Float64
end

function Stim(Is::Float64,mtype::Int64,neur::Int64,node::Int64,tstart::Float64,tstop::Float64)
    Stim([Is],mtype,neur,node,tstart,tstop)
end

type Intracellular
    mtype::Int64
    neur::Int64
    node::Int64
    v::Array{Float64,1}
end

function Intracellular(mtype::Int64,neur::Int64,node::Int64)
    Intracellular(mtype,neur,node,Array(Float64,0))
end

type Network{T <: NeuronPool}
    neur::T #Array of Neurons for simulation
    t::FloatRange{Float64} #time range for simulation
    extra::Array{Extracellular,1} #Extracellular Recording
    #extracellular Stimulation
    intra::Array{Intracellular,1} #Intracellular Recording
    stim::Array{Stim,1} #Intracellular Stimulation
end

function Network(neuron::Neuron,tstop::Float64; par=false)

    #create neuron pool
    gen_neuronpool([neuron],par)

    pool=NeuronPool0([neuron],par)

    Network(pool,0.0:0.025:tstop,Array(Extracellular,0),Array(Intracellular,0),Array(Stim,0))

end

function Network{T<:Neuron}(neurons::Array{T,1},tstop::Float64; par=false)

    #create neuron pool
    gen_neuronpool(neurons,par)

    pool=NeuronPool0(neurons,par)

    Network(pool,0.0:0.025:tstop,Array(Extracellular,0),Array(Intracellular,0),Array(Stim,0))  

end


myconstants=Dict{ASCIIString, Float64}("ena"=>50.0, "ek"=>-77.0)

function make_prop()
end

function gen_prop(a::Channel,k::Int64)
    
    num_fields=sum(length(fieldnames(a)))

    @eval begin
        type $(symbol("Prop_$k")) <: Prop
            p::Array{Float64,2}
        end

        function make_prop(b::$(typeof(a)),n::Int64)
            $(symbol("Prop_$k"))(zeros(Float64,$(num_fields),n))
        end
      
    end

    c=make_prop(a,0)

    @eval begin
       
        function make_prop(a::$(typeof(c)),n::Int64)
            $(symbol("Prop_$k"))(zeros(Float64,$(num_fields),n))
        end

    end

    nothing
    
end

#Combination of channels found
function gen_prop{T<:Tuple}(a::T,k::Int64)

    num_fields=sum([length(fieldnames(a[i]))-1 for i=1:length(a)])

    @eval begin
        type $(symbol("Prop_$k")) <: Prop
            p::Array{Float64,2}
        end

        function make_prop(b::$(typeof(a)),n::Int64)
            $(symbol("Prop_$k"))(zeros(Float64,$(num_fields),n))
        end
          
    end

    c=make_prop(a,0)

    @eval begin
        
        function make_prop(a::$(typeof(c)),n::Int64)
            $(symbol("Prop_$k"))(zeros(Float64,$(num_fields),n))
        end
    end

    nothing

end

function make_neuron()
end

#Neuron with particular combination of channels
function gen_neuron(prop::Prop,k::Int64)

    @eval begin
        type $(symbol("Neuron_$k")) <: Neuron
            soma::($(typeof(prop)))
            axon::($(typeof(prop)))
            dendrite::($(typeof(prop)))
            apical::($(typeof(prop)))
            secstack::Array{Section,1}
            v::Array{Float64,1} #intracellular voltage
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

        function make_neuron(prop::$(typeof(prop)),n::Neuron,newnodes::Array{Node,1})
             $(symbol("Neuron_$k"))(prop,prop,prop,prop,n.secstack,n.v,n.a,n.b,n.d,n.rhs,n.Ra,n.Cm,n.dt,newnodes,n.i_vm,n.divm,n.diag_old,n.internal_nodes,n.par,zeros(Float64,length(n.v)),zeros(Float64,length(n.v)))
        end
        
    end   

end

function gen_neuron(prop::Tuple,k::Int64)

    @eval begin
        type $(symbol("Neuron_$k")) <: Neuron
            soma::($(typeof(prop[1])))
            axon::($(typeof(prop[2])))
            dendrite::($(typeof(prop[3])))
            apical::($(typeof(prop[4])))
            secstack::Array{Section,1}
            v::Array{Float64,1} #intracellular voltage
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

        function make_neuron(prop::$(typeof(prop)),n::Neuron,newnodes::Array{Node,1})
             $(symbol("Neuron_$k"))(prop[1],prop[2],prop[3],prop[4],n.secstack,n.v,n.a,n.b,n.d,n.rhs,n.Ra,n.Cm,n.dt,newnodes,n.i_vm,n.divm,n.diag_old,n.internal_nodes,n.par,zeros(Float64,length(n.v)),zeros(Float64,length(n.v)))
        end
        
    end   

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
