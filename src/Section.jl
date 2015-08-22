#=
Methods to construct and act on Section type, which mirrors the functions of Neuron's Section type
=#

abstract Prop #Property of section (HH, passive etc). contains all of the stuff you need to calc things

type Node{T}
    vars::Dict{ASCIIString,Float64}
    area::Array{Float64,1} #surface area of left [1] and right[2] part of segment
    ri::Array{Float64,1}  #internal resistance of left[1] and right[2] part of segment
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
end

function Node(sec::Section,x::Int64,nseg::Int64)
    
    if length(sec.pnode)>0 #existing nodes
        myvars=sec.pnode[i].vars
        myprops=sec.pnode[i].prop
    else
        myvars=Dict{ASCIIString,Float64}("v"=>0.0)
        myprops=Array(Prop,0)
    end

    #find location of other nodes on section in normalized units
    
    cent=(2*x-1)/(2*nseg)
    interval=Float64[(x-1)/nseg, x/nseg]
    
    #calc area, ri

    #only if pt3d>2. If not, then there is going to need to be special case stuff
    area=zeros(Float64,2)
    ri=zeroes(Float64,2)
    volume=0.0
    
    arc3d=collect(Float64[sec.pt3d[i].arc for i=1:length(sec.pt3d)])

    first=findfirst(arc3d.>interval[1])-1
    mid=findnext(arc3d.>cent)-1
    last=findnext(arc3d.>=interval[2],first)-1

    frac=interp_area(sec.pt3d[first],interval[1],sec.pt3d[first+1])

    diam=frac[1]*sec.pt3d[first].d+frac[2]*sec.pt3d[first+1].d
    height=frac[2]*sec.pt3d[first+1].arc*sec.length
    
    for i=(first+1):mid
        area[1]+=frustrum_area(diam, sec.pt3d[i].d,height)
        ri[1]+=frustrum_resistance(diam,section.pt3d[i].d,height,sec.Ra)
        diam=sec.pt3d[i]
        height=sec.pt3d[i+1]*sec.length
    end

    frac=interp_area(sec.pt3d[mid], cent, sec.pt3d[mid+1])

    diam=frac[1]*sec.pt3d[mid].d+frac[2]*sec.pt3d[mid+1].d
    height=frac[1]*sec.pt3d[mid+1].arc*sec.length

    area[1]+=frustrum_area(sec.pt3d[mid].d,diam,height)
    ri[1]+=frustrum_resistance(sec.pt3d[mid].d,diam,height,sec.Ra)

    diam=sec.pt3d[mid+1].d
    height=frac[2]*sec.pt3d[mid+1].arc*sec.length
    
    for i=(mid+2):last
        area[2]+=frustrum_area(diam,sec.pt3d[i].d,height)
        ri[2]+=frustrum_resistance(diam,sec.pt3d[i].d,height,sec.Ra)
        diam=sec.pt3d[i]
        height=sec.pt3d[i+1]*sec.length
    end

    frac=interp_area(sec.pt3d[last], interval[2], sec.pt3d[last+1])

    diam=frac[1]*sec.pt3d[last].d+frac[2]*sec.pt3d[last+1].d
    height=frac[1]*sec.pt3d[last+1].arc*sec.length

    area[2]+=frustrum_area(sec.pt3d[last].d,diam,height)
    ri[2]+=frustrum_resistance(sec.pt3d[last].d,diam,height,sec.Ra)
     
    #Children are sections where this node is closest to child sections parentx
    child=Array(section,0)

    for i=1:length(sec.child)
        if (sec.child[i].parentx==1.0) & (x==nseg)
            push!(child,sec.child[i])
        elseif sec.child[i].parentx <=interval[2] & sec.child[i].parentx > interval[1] 
            push!(child,sec.child[i])
        end
    end
    
    Node(myvars, area, ri, child, myprops)
     
end

function frustrum_area(d1::Float64,d2::Float64,h::Float64)
    F=pi*(d1/2+d2/2)*sqrt((d1/2-d2/2)^2+h^2)
end

function frustrum_resistance(d1::Float64,d2::Float64,h::Float64,ri::Float64)
    R=p*h/(pi*(d1/2)*(d2/2))
end

function interp_area(x1::Float64, x2::Float64, x3::Float64)
    frac1=(x2-x1)/(x3-x1)
    frac2=(x3-x2)/(x3-x1)

    [frac1,frac2]
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

function Section(section3d::Section3D) #like new_section
    sec=Section(1,section3d.mytype,Array(Node,0),Array(Section,0),Array(Pt3d,size(section3d.xyz,1)),0.0,0.0,0.0)

    #add 3d points from 3d
    for i=1:length(section3d.d)
        if i>0
            sec.pt3d[i]=pt3d(section3d.xyz[i,:]...,section3d.d[i],norm(section3d.xyz[i,:],section3d.xyz[i-1,:]))
            sec.length+=sec.pt3d[i].arc
        else
            sec.pt3d[i]=pt3d(section3d.xyz[i,:]...,section3d.d[i],0.0)
        end
    end

    #normalize arc length
    for i=2:size(sec.pt3d,1)
        sec.pt3d[i].arc=sec.pt3d[i]/sec.length
    end

    sec.parentx=section3d.parentx
    
    sec
end

function change_nseg!(sec::Section,nseg::Int64)

    newnodes=Array(Node,0)
   
    for i=1:nseg
        push!(newnodes,Node(sec,i,nseg))
    end

    sec.pnode=newnodes

    nothing
    
end

function define_shape!(neuron::Neuron)

end

