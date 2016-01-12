
function set_nsegs!(n::Neuron,frequency=100.0,d_lambda=.1)

    n.internal_nodes=[Array(Int64,0) for i=1:4]
    n.nodes=Array(Node,0)
    
    nodesec=zeros(Int64,length(n.secs)+1)
    nseglist=zeros(Int64,length(n.secs))
    nodesec[1]=1
    
    #find total number of segments
    for i=1:length(n.secs)
        
        lambda=lambda_f(frequency, n.secs[i], n)
        nseg=floor(Int64, (n.secs[i].length / (d_lambda * lambda) + .9) /2 )*2+1

        nodesec[i+1]=nodesec[i]+nseg
        
        nseglist[i]=nseg                        
    end

    first=0
    last=0
    for i=1:length(n.secs)

        if n.secs[i].mtype>1

            for j=1:nseglist[i]
           
                (area, ri,mypt3d) = r_a_calc(n.secs[i],j,nseglist[i],n.Ra)

                if (j==1)&&(j==nseglist[i])

                    first=length(n.nodes)+1
                    last=length(n.nodes)+2
                    single_node(n,nodesec,area,ri,mypt3d,i)

                elseif j==1
                    
                    first=length(n.nodes)+1
                    first_node(n,nodesec,area,ri,mypt3d,i)
                    
                elseif j==nseglist[i]

                    last=length(n.nodes)+2
                    last_node(n,area,ri,mypt3d,nodesec,i)
                else
                    middle_node(n,area,ri,mypt3d,i)
                end
            end
        else
            #approximate soma as 1 segment with everything connected at middle
            first=length(n.nodes)+1
            last=length(n.nodes)+3

            parent=length(n.nodes)+2
            children=Array(Int64,0)
            edge_node(n,parent,children,length(n.secs[end].pt3d))
                    
            parent=0
            children=edge_children(n,nodesec,i)            
            (area, ri, mypt3d) = r_a_calc(n.secs[i],1,1,n.Ra)
            push!(children,length(n.nodes))
            push!(children,length(n.nodes)+2)
            internal_node(n,area,ri,parent,children,mypt3d,i)             

            parent=length(n.nodes)
            children=Array(Int64,0)
            edge_node(n,parent,children,1)                    
        end                   
        Section!(n.secs,i,first:last)             
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

edge_children(n,ns,i)=Int64[ns[n.secs[i].child[k]]+n.secs[i].child[k]-1 for k=1:length(n.secs[i].child)]

first_par(n,ns,i)=ns[n.secs[i].parent+1]+n.secs[i].parent-1

first_childs(n)=Int64[length(n.nodes)+2]

first_node(n,ns,area,ri,pt,i)=(p=first_par(n,ns,i); c=first_childs(n); internal_node(n,area,ri,p,c,pt,i))

middle_node(n,area,ri,pt,i)=(p=length(n.nodes);c=first_childs(n);internal_node(n,area,ri,p,c,pt,i))

function last_node(n,area,ri,pt,ns,i)
    p=length(n.nodes)
    c=edge_children(n,ns,i)               
    internal_node(n,area,ri,p,[length(n.nodes)+2],pt,i)           
    edge_node(n,length(n.nodes),c,pt[end])
end

function single_node(n,ns,area,ri,pt,i) 
    p=first_par(n,ns,i)
    if length(n.secs[i].child)==0
        c=Array(Int64,0)
    else
        c=edge_children(n,ns,i)
    end
    internal_node(n,area,ri,p,[length(n.nodes)+2],pt,i)    
    edge_node(n,length(n.nodes),c,pt[end])
end

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
