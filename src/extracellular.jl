
function extracellular{T}(e::Extracellular{T},n::Neuron,sigma::Float64)

    coeffs=zeros(Float64,0)
    inds=zeros(Int64,0)

    for i=1:length(n.secs)
        for j in n.secs[i].pnode
            if n.nodes[j].internal == true

                push!(inds,j)
                mypt3d=n.secs[i].pt3d[n.nodes[j].pt3d]
                myc=calc_coeffs(e,mypt3d,n.nodes[j].parent)
                push!(coeffs,(1/(4*pi*sigma))*myc)
                
            end
        end
    end
      
    Extra_coeffs(coeffs,inds)      
end

calc_coeffs{T<:Point}(e::Extracellular{T},p::Array{Pt3d,1},par::Int64)=point_coeffs(p,e.xyz)

calc_coeffs{T<:Line}(e::Extracellular{T},p::Array{Pt3d,1},par::Int64)=line_coeffs(p,e.xyz)

function calc_coeffs{T<:Mixed}(e::Extracellular{T},p::Array{Pt3d,1},par::Int64)
    if par==0
        c=point_coeffs(p,e.xyz)
    else
        c=line_coeffs(p,e.xyz)
    end
    c
end

function extracellular{T<:Neuron}(n::Array{T,1},e::Extracellular)

    mycoeffs=Array(Extra_coeffs,length(n))
    
    for i=1:length(n)
        mycoeffs[i]=extracellular(e,n[i],0.3)
    end

    mycoeffs   
end

function extracellular{T<:Line}(extra::Extracellular{T},xyz::Array{Float64,2},minds::Array{Array{Int64,1},1},sigma::Float64,cinds::Array{Int64,1})

    coeffs=zeros(Float64,div(size(xyz,1),2))

    i=1
    count=1
    while i<size(xyz,1)
        @inbounds coeffs[count]=(1/(4*pi*sigma))*line_coeffs(xyz,extra.xyz,i)
        i+=2
        count+=1
    end 
    Extra_coeffs(coeffs,cinds)  
end

function line_coeffs(pt3d::Array{Pt3d,1},xyz::Array{Float64,1})
    
    (unit_ds, delta_s)=pt3d_vec(pt3d[1],pt3d[end])

    d1=vect(pt3d[1],xyz)
    d2=vect(pt3d[end],xyz)

    #ln=dot(unit_ds,d1)
    hn=dot(unit_ds,d2)
    ln=delta_s+hn

    if (ln>0) && (hn<0)

        a=sqrt(sum(d2.^2))+abs(hn)
        b=(ln+sqrt(sum(d1.^2))) / abs(sum(d2.^2) - hn*hn)

        d=log(a*b)
        
    else
        
        a=sqrt(sum(d2.^2))+abs(hn)
        b=sqrt(sum(d1.^2))+abs(ln)

        d=log(a/b)

    end

    1/delta_s*d
end

function line_coeffs(xyz::Array{Float64,2},xyze::Array{Float64,1},i::Int64)

    @inbounds (ds_x,ds_y,ds_z,delta_s)=pt3d_vec(xyz[i,1],xyz[i,2],xyz[i,3],xyz[i+1,1],xyz[i+1,2],xyz[i+1,3])

    @inbounds (d1_x,d1_y,d1_z)=vect(xyz[i,1],xyz[i,2],xyz[i,3],xyze[1],xyze[2],xyze[3])
    @inbounds (d2_x,d2_y,d2_z)=vect(xyz[i+1,1],xyz[i+1,2],xyz[i+1,3],xyze[1],xyze[2],xyze[3])

    #ln=pdot(ds_x,ds_y,ds_z,d1_x,d1_y,d1_z)
    hn=pdot(ds_x,ds_y,ds_z,d2_x,d2_y,d2_z)
    ln=delta_s+hn

    if (ln>0) && (hn<0)

        a=sqrt(d2_x*d2_x+d2_y*d2_y+d2_z*d2_z)+hn
        b=(ln+sqrt(d1_x*d1_x+d1_y*d1_y+d1_z*d1_z)) / abs(d2_x*d2_x+d2_y*d2_y+d2_z*d2_z - hn*hn)

        d=log(a*b)

    else
        
        a=sqrt(d2_x*d2_x+d2_y*d2_y+d2_z*d2_z)+abs(hn)
        b=sqrt(d1_x*d1_x+d1_y*d1_y+d1_z*d1_z)+abs(ln)

        d=log(a/b)

    end

    1/delta_s*d
end

function point_coeffs(p::Array{Pt3d,1},xyz::Array{Float64,1})

    middle=round(Int64, length(p)/2)
            
    dist1=dist(p[middle],xyz)

    1/(dist1)          
end

dist(p1::Pt3d, p2::Pt3d)=sqrt((p1.x-p2.x)^2+(p1.y-p2.y)^2+(p1.z-p2.z)^2)

dist(p1::Pt3d,xyz::Array{Float64,1})=sqrt((p1.x-xyz[1])^2+(p1.y-xyz[2])^2+(p1.z-xyz[3])^2)

dist(x1::Float64,y1::Float64,z1::Float64,x2::Float64,y2::Float64,z2::Float64)=sqrt((x1-x2)^2+(y1-y2)^2+(z1-z2)^2)

vect(p::Pt3d,xyz::Array{Float64,1})=[xyz[1]-p.x,xyz[2]-p.y,xyz[3]-p.z]

vect(p1::Pt3d,p2::Pt3d)=[p2.x-p1.x,p2.y-p1.y,p2.z-p1.z]

vect(x1::Float64,y1::Float64,z1::Float64,x2::Float64,y2::Float64,z2::Float64)=(x2-x1,y2-y1,z2-z1)

pt3d_vec(p1::Pt3d,p2::Pt3d)=(n=vect(p1,p2); mag=dist(p1,p2); n /= mag; (n, mag))

function pt3d_vec(x1::Float64,y1::Float64,z1::Float64,x2::Float64,y2::Float64,z2::Float64)
    (vx,vy,vz)=vect(x1,y1,z1,x2,y2,z2)
    mag=dist(x1,y1,z1,x2,y2,z2)
    (vx/mag,vy/mag,vz/mag,mag)   
end

pdot(x1,y1,z1,x2,y2,z2)=x1*x2+y1*y2+z1*z2

function getv(im::Array{Float64,2},e::Extra_coeffs,v::Array{Float64,2},k::Int64)
    
    mystart=20
    
    for i=1:size(v,1)      
        for j=1:length(e.c)
            @inbounds v[i,k]+=e.c[j]*im[e.inds[j],mystart+i-1] 
        end     
    end
    nothing
end

function getv(im::Array{Float64,2},e::Extra_coeffs,v::Array{Float64,3},k::Int64,l::Int64)
    
    mystart=20
    
    for i=1:size(v,1)      
        for j=1:length(e.c)
            @inbounds v[i,k,l]+=e.c[j]*im[e.inds[j],mystart+i-1] 
        end     
    end
    nothing
end

function nete(n::Neuron,im::Array{Float64,2},ex::Extracellular,num::Int64)
     
    mystart=20
    v=zeros(Float64,length(mystart:size(im,2)),num)    
    (xyz,inds)=getxyz(n);   
    cinds=findcind(n)   
    xyz1=copyn(xyz)
    
    for i=1:num
        resetn!(xyz,xyz1)
        randomize_shape!(n,xyz1,inds)
        translate3d!(xyz1,rand(-500.0:500.0),rand(-500.0:500.0),rand(-500.0:500.0))        
        c=extracellular(ex,xyz1,inds,.3,cinds)
        getv(im,c,v,i)      
    end   
    v
end

function nete{T}(n::Neuron,im::Array{Float64,2},ex::Array{T,1},num::Int64)

    mystart=20
    v=zeros(Float64,length(mystart:size(im,2)),num,length(ex))    
    (xyz,inds)=getxyz(n);   
    cinds=findcind(n)   
    xyz1=copyn(xyz)
    
    for i=1:num
        resetn!(xyz,xyz1)
        randomize_shape!(n,xyz1,inds)
        translate3d!(xyz1,rand(-500.0:500.0),rand(-500.0:500.0),rand(-500.0:500.0))
        for j=1:length(ex)
            c=extracellular(ex[j],xyz1,inds,.3,cinds)
            getv(im,c,v,i,j)  
        end       
    end   
    v
end

function extrap(v::Array{Float64,2},t::Float64)

    v2=zeros(Float64,length(0.0:.025:t))

    spike=find_cspikes(v)
    if length(spike)>0
        cur_spike=spike[1]
        spike_ind=1
    else
        cur_spike=0
    end

    st=[Array(Int64,0) for i=1:length(spike)]

    firing=zeros(Int64,length(0.0:.025:t)-size(v,1))

    @inbounds for i=1:size(v,2)

        d=Poisson(100*rand()/40000)

        rand!(d,firing)

        j=1
        if i==cur_spike
            while j<length(firing)
                if firing[j]>0
                    push!(st[spike_ind],j)
                    for k=1:size(v,1)
                        if abs(v[k,i])>abs(v2[j])
                            v2[j]=v[k,i]+randn()*.01
                        end
                        j+=1
                    end
                end
                j+=1
            end
            spike_ind+=1
            if spike_ind<=length(spike)
                cur_spike=spike[spike_ind]
            end
        else
            while j<length(firing)
                if firing[j]>0
                    for k=1:size(v,1)
                        if abs(v[k,i])>abs(v2[j])
                            v2[j]=v[k,i]+randn()*.01
                        end
                        j+=1
                    end
                end
                j+=1
            end
        end
    end
    (v2,st)
end

function extrap(v::Array{Float64,3},t::Float64)

    v2=zeros(Float64,length(0.0:.025:t),size(v,3))

    spike=find_cspikes(mean(v,3))
    if length(spike)>0
        cur_spike=spike[1]
        spike_ind=1
    else
        cur_spike=0
    end

    st=[Array(Int64,0) for i=1:length(spike)]

    firing=zeros(Int64,length(0.0:.025:t)-size(v,1))

    @inbounds for i=1:size(v,2)

        d=Poisson(100*rand()/40000)

        rand!(d,firing)

        j=1
        if i==cur_spike
            while j<length(firing)
                if firing[j]>0
                    push!(st[spike_ind],j)
                    for k=1:size(v,1)
                        for l=1:size(v2,2)
                            if abs(v[k,i,l])>abs(v2[j,l])
                                v2[j,l]=v[k,i,l]+randn()*.01
                            end
                        end
                        j+=1
                    end
                end
                j+=1
            end
            spike_ind+=1
            if spike_ind<=length(spike)
                cur_spike=spike[spike_ind]
            end
        else
            while j<length(firing)
                if firing[j]>0
                    for k=1:size(v,1)
                        for l=1:size(v2,2)
                            if abs(v[k,i,l])>abs(v2[j,l])
                                v2[j,l]=v[k,i,l]+randn()*.01
                            end
                        end
                        j+=1
                    end
                end
                j+=1
            end
        end
    end
    (v2,st)
end

function find_cspikes(v)

    spike_num=Array(Int64,0)  
    l=size(v,1)
    for i=1:size(v,2)
        for j=1:l
            @inbounds if v[j,i]<-0.12
                push!(spike_num,i)
                break
            end
        end 
    end
    spike_num
end

function runc(n::NetworkS,init=false)

    myi=zeros(Float64,length(n.neur.N_1[1].i_vm),length(n.t))    
    if init==true
        init!(n)
    end  
    for i=1:5
	n.neur.N_1[1].v[i]=30.0
    end
    for i=0:2
	n.neur.N_1[1].v[end-i]=30
    end                       
    for i=1:length(n.t)
	
        @inbounds main(n.neur)
        
        myi[:,i]=n.neur.N_1[1].i_vm[:]          
    end
    myi       
end
