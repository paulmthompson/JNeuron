#=
Methods to construct and act on Section type, which mirrors the functions of Neuron's Section type
=#

abstract Prop #Property of section (HH, passive etc). contains all of the stuff you need to calc things

type Node
    vars::Dict{ASCIIString,Float64}
    area::Float64 #area of section at node
    a
    b
    ri #resistance between center of segment and parent segment
    v_temp
    d::Float64 #diameter of section at node
    rhs
    a_matelm
    b_matelm
    child::Section #section connected to this node?
    sec::Section #this section
    prop::Array{Prop,1} #Array of abstract types (yucky!), each a subtype of 
end

function Node()
    Node()
end

function add_prop!(node::Node,prop::Prop)

    if sum([typeof(node.prop[i])==typeof(prop) for i=1:length(node.prop)])==0

        push!(node.prop,prop)
        
        for i=1:length(prop.nodevar)
            if !haskey(node.vars,prop.nodevar[i])
                node.vars[prop.nodevar[i]]=0.0
            end
        end
    else
        println("Property also exists in node. Not inserting")
    end
           
    nothing        
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
    child::Array{Section,0}
    pt3d::Array{Pt3d,1}
    Ra::Float64 #cytoplasmic resistivity
end

function Section(section3d::Section3D) #like new_section
    sec=Section(1,section3d.mytype,Array(Node,0),Array(Section,0),Array(Pt3d,size(section3d.xyz,1)),)

    mylength=0.0
    #add 3d points from 3d
    for i=1:length(section3d.d)
        if i>0
            mylength+=sec.pt3d+
            sec.pt3d[i]=pt3d(section3d.xyz[i,:]...,section3d.d[i],norm(section3d.xyz[i,:],section3d.xyz[i-1,:]))
        else
            sec.pt3d[i]=pt3d(section3d.xyz[i,:]...,section3d.d[i],0.0)
        end
    end

    #normalize arc length
    for i=2:size(sec.pt3d,1)
        sec.pt3d[i].arc=sec.pt3d[i]/mylength
    end
    
    sec
end

function change_nseg!(sec::Section)
    #use triangular integration to find diameter, area, arc, (ri?) at each node point

end

function define_shape!(neuron::Neuron)

end

