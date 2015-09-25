
abstract Prop #Property of section (HH, passive etc). contains all of the stuff you need to calc things

type Node{T}
    vars::Dict{ASCIIString,Float64}
    area::Array{Float64,1} #surface area of left [1] and right[2] part of segment
    ri::Array{Float64,1}  #internal resistance of left[1] and right[2] part of segment
    parent::T
    children::Array{T,1} #Node(s) from other sections attached to this one
    prop::Array{Prop,1} #Array of abstract types (yucky!), each a subtype of 
end

#associated 3d point
type Pt3d
    x::Float64
    y::Float64
    z::Float64
    d::Float64
    arc::Float64 #normalized distance from 0 to end
end

type Section
    refcount::Int64 #ID for section, also the place in secstack array
    mtype::Int64 #Cellbody=1,Axon=2,Dendrite=3,Apical=4
    pnode::Array{Node,1} #one node at center of each segment
    child::Array{Section,1}
    pt3d::Array{Pt3d,1}
    Ra::Float64 #cytoplasmic resistivity
    parentx::Float64
    length::Float64
end

type Neuron
    secstack::Array{Section,1}
    A::Array{Float64,2}
    V_new::Array{Float64,1}
    V_old::Array{Float64,1}
    delta_V::Array{Float64,1}
    rhs::Array{Float64,1}
    Ra::Float64
    Cm::Float64
    dt::Float64
end
