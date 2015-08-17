#=
Methods to construct and act on Section type, which mirrors the functions of Neuron's Section type
=#

type Node
    vars::Dict{ASCIIString,Float64}
    area
    a
    b
    rinv
    v_temp
    d
    rhs
    a_matelm
    b_matelm
    child::Section #section connected to this node
    sec::Section #this section
    prop::Array{Prop,1} #Array of abstract types (yucky!), each a subtype of 
end

function Node()
    Node()
end

function add_prop(node::Node,prop::Prop)
    for i=1:length(prop.nodevar)
        if !haskey(node.vars,prop.nodevar[i])
            node.vars[prop.nodevar[i]]=0.0
        end
    end
    nothing        
end

#associated 3d point
type Pt3d
    x
    y
    z
    d
    arc
end

abstract Prop #Property of section (HH, passive etc). contains all of the stuff you need to calc things

#=
In neuron, prop is a structure that is filled in by compiling the mod files for each channel (I think).

Everything in the mod files that is defined as a PARAMETER, STATE, or DERIVATIVE ends up being a part of param
Additionally, the last two entries in param seem to always be v, and _g (conductance)?
Things in param can be shared with other Props (like voltage is equal at every Prop for each node)

ppvar con contains other values, looks like usually related to ions

So as far as Julia goes, at each node there are state variables unique to one ion channel(e.g. sodium gate m), state variables shared among ion channels (e.g. voltage), constant paramters for each channel
=#

type Section
    refcount::Int64 #ID for section, also the place in secstack array
    nnode::Int64 #Number of nodes (nseg+1)
    pnode::Array{Node,1} #nseg+1
    parentsec::Section
    child::Section
    sibling::Section
    parentnode::Node
    order::Int64
    recalc_area::Int64
    volatile_mark::Int64
    npt3d::Int64
    pt3d::Array{Pt3d,1}
end

function Section(n::nseg)
end

function Section(section3d::Section3D) #like new_section
end

