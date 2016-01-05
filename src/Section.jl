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

        first=1
        last=length(arc3d)
            
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
    
    (area, ri,sec.pt3d[first:last])
     
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

function add(neuron::Neuron,prop_array::Channel)

    global num_neur::Int64
    global num_prop::Int64

    if method_exists(make_prop,(typeof(prop_array),Int))
    else
        num_prop+=1
        gen_prop(prop_array,num_prop)
    end
    
    myprop=make_prop(prop_array,0)

    gen_current(prop_array,myprop)

    newnodes=Array(Node,length(neuron.nodes))
    
    if method_exists(make_neuron,(typeof(myprop),typeof(neuron),typeof(newnodes)))
    else     
        num_neur+=1
        gen_neuron(myprop,num_neur)
    end

    for i=1:length(neuron.nodes)
        newnodes[i]=Node(neuron.nodes[i],make_prop(prop_array,0))

        for j=2:length(fieldnames(myprop))
            for k=1:length(getfield(myprop,j).nodevar)
                if !haskey(newnodes[i].vars,getfield(myprop,j).nodevar[k])
                    newnodes[i].vars[getfield(myprop,j).nodevar[k]]=0.0
                end
            end
        end
    end

    n=deepcopy(neuron)
    
    n1=make_neuron(myprop,n,newnodes)

    reset_pnode!(n1)

    n1
    
end

function add{T<:Tuple}(neuron::Neuron,prop_array::T)

    global num_neur::Int64
    global num_prop::Int64

    if method_exists(make_prop,(typeof(prop_array),Int))
    else
        num_prop+=1
        gen_prop(prop_array,num_prop)
    end
    
    myprop=make_prop(prop_array,0)

    gen_current(prop_array,myprop)

    newnodes=Array(Node,length(neuron.nodes))
    
    if method_exists(make_neuron,(typeof(myprop),typeof(neuron),typeof(newnodes)))
    else     
        num_neur+=1
        gen_neuron(myprop,num_neur)
    end

    for i=1:length(neuron.nodes)
        newnodes[i]=Node(neuron.nodes[i],make_prop(prop_array,0))

        for j=2:length(fieldnames(myprop))
            for k=1:length(getfield(myprop,j).nodevar)
                if !haskey(newnodes[i].vars,getfield(myprop,j).nodevar[k])
                    newnodes[i].vars[getfield(myprop,j).nodevar[k]]=0.0
                end
            end
        end
    end

    n=deepcopy(neuron)
    
    n1=make_neuron(myprop,n,newnodes)

    reset_pnode!(n1)

    n1
    
end

function add(neuron::Neuron,prop1,prop2,prop3,prop4)

    global num_neur::Int64
    global num_prop::Int64

    prop_array=(prop1,prop2,prop3,prop4)

    for i=1:4
        if method_exists(make_prop,(typeof(prop_array[i]),Int))
        else
            num_prop+=1
            gen_prop(prop_array[i],num_prop)
        end
    end
     
    new_prop_array=(make_prop(prop1,0),make_prop(prop2,0),make_prop(prop3,0),make_prop(prop4,0))

    gen_current(prop1,new_prop_array[1])
    gen_current(prop2,new_prop_array[2])
    gen_current(prop3,new_prop_array[3])
    gen_current(prop4,new_prop_array[4])

    newnodes=Array(Node,length(neuron.nodes))
    
    if method_exists(make_neuron,(typeof(new_prop_array),typeof(neuron),typeof(newnodes)))
    else     
        num_neur+=1
        gen_neuron(new_prop_array,num_neur)
    end

    for i=1:length(neuron.secstack)
        for j=1:length(neuron.secstack[i].pnode)
            ind=neuron.secstack[i].pnode[j].ind
            mtype=neuron.secstack[i].mtype

            newnodes[ind]=Node(neuron.secstack[i].pnode[j],make_prop(prop_array[mtype],0))

            for k=2:length(fieldnames(new_prop_array[mtype]))
                for h=1:length(getfield(new_prop_array[mtype],k).nodevar)
                    if !haskey(newnodes[ind].vars,getfield(new_prop_array[mtype],k).nodevar[h])
                        newnodes[ind].vars[getfield(new_prop_array[mtype],k).nodevar[h]]=0.0
                    end
                end
            end
        end
    end

    n=deepcopy(neuron)
    
    n1=make_neuron(new_prop_array,n,newnodes)

    reset_pnode!(n1)

    n1
    
end

function reset_pnode!(neuron::Neuron)

    for i=1:length(neuron.secstack)

        first=neuron.secstack[i].pnode[1].ind

        last=neuron.secstack[i].pnode[end].ind

        neuron.secstack[i].pnode=sub(neuron.nodes,first:last)
        
    end

    nothing
    
end


function Node{T<:Prop}(node::Node{Prop0},myprop::T)
    Node(node.ind,node.vars,node.area,node.ri,node.b,node.a,node.parent,node.children,node.internal,node.pt3d,deepcopy(myprop))
end


