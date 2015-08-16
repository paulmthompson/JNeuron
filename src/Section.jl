#=
Methods to construct and act on Section type, which mirrors the functions of Neuron's Section type
=#

type Node
    v
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
end

function Node()
    Node()
end

abstract Prop #Property of section (HH, passive etc). contains all of the stuff you need to calc things

type Section{T <: Prop}
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
    prop::T #property list
end

function Section{T<:Prop}(prop::T,n::nseg)
end

type Pt3d
    x
    y
    z
    d
    arc
end

function new_section(section3d::Section3D) # don't know if i need this rather than a Section constructor
    sec=Section()

    #add to secstack array
    
    return sec
end

