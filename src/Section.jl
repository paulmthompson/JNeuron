#=
Methods to construct and act on Section type, which mirrors the functions of Neuron's Section type

Inputs:
sec - section containing node
x - ID of node to be created (range of 1 to nseg)
nseg - number of total segments in section

=#

function r_a_calc(sec::Section,x::Int64,nseg::Int64)
    
    #find location of other nodes on section in normalized units
    
    cent=(2*x-1)/(2*nseg)
    interval=Float64[(x-1)/nseg, x/nseg]
    
    #calc area, ri

    #only if pt3d>2. If not, then there is going to need to be special case stuff
    area=zeros(Float64,2)
    ri=zeros(Float64,2)

    #get all of the arc lengths into one data structures
    arc3d=collect(Float64[sec.pt3d[i].arc for i=1:length(sec.pt3d)])

    #Find the first 3d points 1) before the segment, 2) after the middle of the segment and 3) before the end of the segment
    first=findfirst(arc3d.>interval[1])-1
    mid=findfirst(arc3d.>cent)-1
    last=findnext(arc3d.>=interval[2],first)-1

    if length(arc3d)<4
        #If just a few points, whole thing is just modeled as a cylinder
        myarea=2*pi*sec.pt3d[1].d/2*sec.length*(interval[2]-interval[1])

        area[1]=myarea/2
        area[2]=myarea/2

        ri[1]=frustrum_resistance(sec.pt3d[1].d,sec.pt3d[1].d,sec.length*(interval[2]-interval[1])/2,sec.Ra)
        ri[2]=ri[1]
            
    else
        

    #=
        Area of a segment is modeled as a series of truncated cones between the 3d points. The diameter of the top and bottom of the frustrum is taken to be the average of the diameters of the nearest 3d points
    =#


    #left side to first point before center
    frac=interp_area(arc3d[first],interval[1],arc3d[first+1])

    diam=frac[1]*sec.pt3d[first].d+frac[2]*sec.pt3d[first+1].d
    height=(arc3d[first+1]-interval[1])*sec.length
    
    for i=(first+1):mid
        area[1]+=frustrum_area(diam, sec.pt3d[i].d,height)
        ri[1]+=frustrum_resistance(diam,sec.pt3d[i].d,height,sec.Ra)
        diam=sec.pt3d[i].d
        height=(arc3d[i+1]-arc3d[i])*sec.length
    end

    #First point before center to center
    frac=interp_area(arc3d[mid], cent, arc3d[mid+1])

    diam=frac[1]*sec.pt3d[mid].d+frac[2]*sec.pt3d[mid+1].d
    height=(cent-arc3d[mid])*sec.length

    area[1]+=frustrum_area(sec.pt3d[mid].d,diam,height)
    ri[1]+=frustrum_resistance(sec.pt3d[mid].d,diam,height,sec.Ra)

    #first point after center    
    height=(arc3d[mid+1]-cent)*sec.length

    area[2]+=frustrum_area(diam,sec.pt3d[mid+1].d,height)
    ri[2]+=frustrum_resistance(diam,sec.pt3d[mid+1].d,height,sec.Ra)

    diam=sec.pt3d[mid+1].d
    height=(arc3d[mid+2]-arc3d[mid+1])*sec.length

    for i=(mid+2):last
        area[2]+=frustrum_area(diam,sec.pt3d[i].d,height)
        ri[2]+=frustrum_resistance(diam,sec.pt3d[i].d,height,sec.Ra)
        diam=sec.pt3d[i].d
        height=(arc3d[i+1]-arc3d[i])*sec.length
    end

    #Final point
    frac=interp_area(arc3d[last], interval[2], arc3d[last+1])

    diam=frac[1]*sec.pt3d[last].d+frac[2]*sec.pt3d[last+1].d
    height=(interval[2]-arc3d[last])*sec.length

    area[2]+=frustrum_area(sec.pt3d[last].d,diam,height)
    ri[2]+=frustrum_resistance(sec.pt3d[last].d,diam,height,sec.Ra)

    end

    ri*=.01
    
    (area, ri)
     
end

function frustrum_area(d1::Float64,d2::Float64,h::Float64)
    F=pi*(d1/2+d2/2)*sqrt((d1/2-d2/2)^2+h^2)
end

function frustrum_resistance(d1::Float64,d2::Float64,h::Float64,ri::Float64)
    R=ri*h/(pi*(d1/2)*(d2/2))
end

function interp_area(x1::Float64, x2::Float64, x3::Float64)
    frac1=(x2-x1)/(x3-x1)
    frac2=(x3-x2)/(x3-x1)

    [frac1,frac2]
end

function add_prop!(node::Node,prop::Prop)

    if sum([typeof(node.prop[i])==typeof(prop) for i=1:length(node.prop)])==0

        push!(node.prop,typeof(prop)())
        
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

function add_prop!(neuron::Neuron,prop::Prop)

    for j=1:length(neuron.nodes)

        if neuron.nodes[j].internal==true
    
            if sum([typeof(neuron.nodes[j].prop[i])==typeof(prop) for i=1:length(neuron.nodes[j].prop)])==0

                push!(neuron.nodes[j].prop,typeof(prop)())
        
                for i=1:length(prop.nodevar)
                    if !haskey(neuron.nodes[j].vars,prop.nodevar[i])
                    neuron.nodes[j].vars[prop.nodevar[i]]=0.0
                    end
                end
            else
                println("Property also exists in node. Not inserting")
            end

        end
        
    end
             
    nothing
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

