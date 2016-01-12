
function set_nsegs!(neuron::Neuron,frequency=100.0,d_lambda=.1)

    neuron.internal_nodes=[Array(Int64,0) for i=1:4]
    neuron.nodes=Array(Node,0)
    
    nodesec=zeros(Int64,length(neuron.secstack)+1)
    nseglist=zeros(Int64,length(neuron.secstack))
    nodesec[1]=1
    first=0
    last=0

    #find total number of segments
    for i=1:length(neuron.secstack)
        
        lambda=lambda_f(frequency, neuron.secstack[i], neuron)
        nseg=floor(Int64, (neuron.secstack[i].length / (d_lambda * lambda) + .9) /2 )*2+1

        nodesec[i+1]=nodesec[i]+nseg
        
        nseglist[i]=nseg                        
    end
    
    for i=1:length(neuron.secstack)

        if neuron.secstack[i].mtype>1

            for j=1:nseglist[i]
           
                (area, ri,mypt3d) = r_a_calc(neuron.secstack[i],j,nseglist[i],neuron.Ra)
                
                if j==1
                    parent=nodesec[neuron.secstack[i].parent+1]+neuron.secstack[i].parent-1
                    if nseglist[i]==1
                        if length(neuron.secstack[i].child)==0
                            children=Array(Int64,0)
                        else
                            children=Int64[nodesec[neuron.secstack[i].child[k]]+neuron.secstack[i].child[k]-1 for k=1:length(neuron.secstack[i].child)]
                        end
                        
                    else
                        children=Int64[length(neuron.nodes)+2]
                    end
                elseif j==nseglist[i]
                    parent=length(neuron.nodes)
                    children=Int64[nodesec[neuron.secstack[i].child[k]]+neuron.secstack[i].child[k]-1 for k=1:length(neuron.secstack[i].child)]
                else
                    parent=length(neuron.nodes)
                    children=Int64[length(neuron.nodes)+2]
                end

                if (j==1)&&(j==nseglist[i])

                    first=length(neuron.nodes)+1
                    last=length(neuron.nodes)+2

                    internal_node(neuron,area,ri,parent,[length(neuron.nodes)+2],mypt3d,i)

                    edge_node(neuron,length(neuron.nodes),children,mypt3d[end])                
                elseif j==1
                    
                    first=length(neuron.nodes)+1

                    internal_node(neuron,area,ri,parent,children,mypt3d,i)                           
                elseif j==nseglist[i]

                    last=length(neuron.nodes)+2
                
                    internal_node(neuron,area,ri,parent,[length(neuron.nodes)+2],mypt3d,i)
                    
                    edge_node(neuron,length(neuron.nodes),children,mypt3d[end])
                else
                
                    internal_node(neuron,area,ri,parent,children,mypt3d,i)
                end
            end
        else

            #approximate soma as 1 segments with everything connected at middle
            first=length(neuron.nodes)+1
            last=length(neuron.nodes)+3
            for j=1:3

                if j ==1 
                    parent=length(neuron.nodes)+2
                    children=Array(Int64,0)

                    edge_node(neuron,parent,children,length(neuron.secstack[end].pt3d))
                    
                elseif j == 2

                    parent=0
                    children=Int64[nodesec[neuron.secstack[i].child[k]]+neuron.secstack[i].child[k]-1 for k=1:length(neuron.secstack[i].child)]
                    
                    (area, ri, mypt3d) = r_a_calc(neuron.secstack[i],1,1,neuron.Ra)
                    push!(children,length(neuron.nodes))
                    push!(children,length(neuron.nodes)+2)

                    internal_node(neuron,area,ri,parent,children,mypt3d,i)
                    
                else
                    parent=length(neuron.nodes)
                    children=Array(Int64,0)

                    edge_node(neuron,parent,children,1)                    
                end              
            end
        end                   
        Section!(neuron.secstack,i,first:last)             
    end
    nothing   
end

function lambda_f(frequency::Float64,sec::Section,neuron::Neuron)
    
    x1 = 0.0
    d1 = sec.pt3d[1].d
    lam = 0.0

    for i=2:length(sec.pt3d)
        x2=sec.pt3d[i].arc*sec.length
        d2=sec.pt3d[i].d
        lam+= (x2-x1)/sqrt(d1+d2)
        x1=x2
        d1=d2
    end

    lam *= sqrt(2) * 1e-5 * sqrt(4*pi*frequency*neuron.Ra*neuron.Cm)

    sec.length/lam
    
end

internal_node(n,a,ri,p,c,pt,i)=(nd=Node(a,ri,p,c,pt); add_node(n,nd); add_inode(n,i))

edge_node(n,p,c,pt)=(nd=Node(p,c,pt); add_node(n,nd))

function r_a_calc(sec::Section,x::Int64,nseg::Int64,Ra::Float64)
    
    #find location of other nodes on section in normalized units
    
    cent=(2*x-1)/(2*nseg)
    interval=Float64[(x-1)/nseg, x/nseg]
    
    #get all of the arc lengths into one data structures
    arc3d=collect(Float64[sec.pt3d[i].arc for i=1:length(sec.pt3d)])

    #Find the first 3d points 1) before the segment, 2) after the middle of the segment and 3) before the end of the segment
    first=findfirst(arc3d.>interval[1])-1
    mid=findfirst(arc3d.>cent)-1
    last=findnext(arc3d.>=interval[2],first)-1

    if length(arc3d)<4
        
        #If just a few points, whole thing is just modeled as a cylinder
        myarea=2*pi*sec.pt3d[1].d/2*sec.length*(interval[2]-interval[1])

        area=SegArea(myarea/2,myarea/2,myarea)

        r=frustrum_resistance(sec.pt3d[1].d,sec.pt3d[1].d,sec.length*(interval[2]-interval[1])/2,Ra)
        ri=SegRi(.01*r,.01*r)

        first=1
        last=length(arc3d)
            
    else

        al=0.0
        ar=0.0
        rl=0.0
        rr=0.0    
        
    #=
        Area of a segment is modeled as a series of truncated cones between the 3d points. The diameter of the top and bottom of the frustrum is taken to be the average of the diameters of the nearest 3d points
    =#

        #left side to first point before center
        frac=interp_area(arc3d[first],interval[1],arc3d[first+1])

        diam=frac[1]*sec.pt3d[first].d+frac[2]*sec.pt3d[first+1].d
        height=(arc3d[first+1]-interval[1])*sec.length
    
        for i=(first+1):mid
            al+=frustrum_area(diam, sec.pt3d[i].d,height)
            rl+=frustrum_resistance(diam,sec.pt3d[i].d,height,Ra)
            diam=sec.pt3d[i].d
            height=(arc3d[i+1]-arc3d[i])*sec.length
        end

        #First point before center to center
        frac=interp_area(arc3d[mid], cent, arc3d[mid+1])

        diam=frac[1]*sec.pt3d[mid].d+frac[2]*sec.pt3d[mid+1].d
        height=(cent-arc3d[mid])*sec.length

        al+=frustrum_area(sec.pt3d[mid].d,diam,height)
        rl+=frustrum_resistance(sec.pt3d[mid].d,diam,height,Ra)

        #first point after center    
        height=(arc3d[mid+1]-cent)*sec.length

        ar+=frustrum_area(diam,sec.pt3d[mid+1].d,height)
        rr+=frustrum_resistance(diam,sec.pt3d[mid+1].d,height,Ra)

        diam=sec.pt3d[mid+1].d
        height=(arc3d[mid+2]-arc3d[mid+1])*sec.length

        for i=(mid+2):last
            ar+=frustrum_area(diam,sec.pt3d[i].d,height)
            rr+=frustrum_resistance(diam,sec.pt3d[i].d,height,Ra)
            diam=sec.pt3d[i].d
            height=(arc3d[i+1]-arc3d[i])*sec.length
        end

        #Final point
        frac=interp_area(arc3d[last], interval[2], arc3d[last+1])

        diam=frac[1]*sec.pt3d[last].d+frac[2]*sec.pt3d[last+1].d
        height=(interval[2]-arc3d[last])*sec.length

        ar+=frustrum_area(sec.pt3d[last].d,diam,height)
        rr+=frustrum_resistance(sec.pt3d[last].d,diam,height,Ra)

        area=SegArea(al,ar,al+ar)
        ri=SegRi(.01*rl,.01*rr)
    end
    
    (area, ri, first:last)
     
end

frustrum_area(d1::Float64,d2::Float64,h::Float64)=pi*(d1/2+d2/2)*sqrt((d1/2-d2/2)^2+h^2)

frustrum_resistance(d1::Float64,d2::Float64,h::Float64,ri::Float64)=ri*h/(pi*(d1/2)*(d2/2))

interp_area(x1::Float64, x2::Float64, x3::Float64)=((x2-x1)/(x3-x1),(x3-x2)/(x3-x1))
