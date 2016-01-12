

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

function add(neuron::Neuron,prop_array...)

    num_arg=length(prop_array)

    if num_arg==1

        prop_array=prop_array[1]
        
        gen_prop_check(prop_array)
        
        myprop=make_prop(prop_array,0)
        gen_current(prop_array,myprop)
        
    elseif num_arg==4
        
        for i=1:4
            gen_prop_check(prop_array[i])
        end

        myprop=(make_prop(prop_array[1],0),make_prop(prop_array[2],0),make_prop(prop_array[3],0),make_prop(prop_array[4],0))

        for i=1:4
            gen_current(prop_array[i],myprop[i])
        end

    else
        #error checking
    end

    gen_neur_check(myprop)
     
    n1=make_neuron(myprop,neuron.nodes,neuron.secstack,neuron.internal_nodes)
    
end



