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

function add_prop!(node::Node,prop::Prop)
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
    mtype::Int64 #Cellbody=1,Axon=2,Dendrite=3,Apical=4
    nnode::Int64 #Number of nodes (nseg+1)
    pnode::Array{Node,1} #nseg+1
    child::Array{Section,0}
    npt3d::Int64
    pt3d::Array{Pt3d,1}
end

function Section(section3d::Section3D) #like new_section
    sec=Section(1,0,section3d.mytype,Array(Node,0),Array(Section,0),size(section3d.xyz,1),Array(Pt3d,size(section3d.xyz,1)))
    
    #add 3d points from 3d
    for i=1:length(section3d.d)
        sec.pt3d[i]=pt3d(section3d.xyz[i,:]...,section3d.d[i],i/length(section3d.d))
    end
    sec
end

function change_nseg!()

end

function define_shape!(neuron::Neuron)

end

