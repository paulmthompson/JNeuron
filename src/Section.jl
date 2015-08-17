#=
Methods to construct and act on Section type, which mirrors the functions of Neuron's Section type
=#

abstract Prop #Property of section (HH, passive etc). contains all of the stuff you need to calc things

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

