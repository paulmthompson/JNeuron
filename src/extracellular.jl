#=

Extracellular 
(x,y,z) extracellular position
sigma: extracellular conductivity
=#

#Point Source

function extracellular{T}(e::Extracellular{T},n::Neuron,sigma::Float64)

    coeffs=zeros(Float64,0)
    inds=zeros(Int64,0)

    for i=1:length(n.secs)
        for j in n.secs[i].pnode
            if n.nodes[j].internal == true

                push!(inds,j)
                mypt3d=n.secs[i].pt3d[n.nodes[j].pt3d]
                myc=calc_coeffs(e,mypt3d,n)
                push!(coeffs,(1/(4*pi*sigma))*myc)
                
            end
        end
    end
      
    Extra_coeffs(coeffs,inds)      
end

calc_coeffs{T<:Point}(e::Extracellular{T},p::Array{Pt3d,1},n::Neuron)=point_coeffs(p,e.xyz)

calc_coeffs{T<:Line}(e::Extracellular{T},p::Array{Pt3d,1},n::Neuron)=line_coeffs(p,e.xyz)

function calc_coeffs{T<:Mixed}(e::Extracellular{T},p::Array{Pt3d,1},n::Neuron)
    if n.nodes[j].parent==0
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

    ln=dot(unit_ds,d1)
    hn=dot(unit_ds,d2)

    a=sqrt(sum(d2.^2))-hn
    b=sqrt(sum(d1.^2))-ln

    1/delta_s*log(abs(a/b))  
end

function line_coeffs(xyz::Array{Float64,2},xyze::Array{Float64,1},i::Int64)

    @inbounds (ds_x,ds_y,ds_z,delta_s)=pt3d_vec(xyz[i,1],xyz[i,2],xyz[i,3],xyz[i+1,1],xyz[i+1,2],xyz[i+1,3])

    @inbounds (d1_x,d1_y,d1_z)=vect(xyz[i,1],xyz[i,2],xyz[i,3],xyze[1],xyze[2],xyze[3])
    @inbounds (d2_x,d2_y,d2_z)=vect(xyz[i+1,1],xyz[i+1,2],xyz[i+1,3],xyze[1],xyze[2],xyze[3])

    ln=pdot(ds_x,ds_y,ds_z,d1_x,d1_y,d1_z)
    hn=pdot(ds_x,ds_y,ds_z,d2_x,d2_y,d2_z)

     a=sqrt(d2_x*d2_x+d2_y*d2_y+d2_z*d2_z)-hn
     b=sqrt(d1_x*d1_x+d1_y*d1_y+d1_z*d1_z)-ln

    1/delta_s*log(abs(a/b))
    
end

function point_coeffs(p::Array{Pt3d,1},xyz::Array{Float64,1})

    middle=round(length(p)/2)
            
    dist1=vect(p[middle],xyz)

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

function getv(i,e,v,l)
    
    mystart=round(Int,size(i,2)/2)
    
    for k=1:size(v,1)      
        for j=1:length(e.c)
            @inbounds v[k,l]+=e.c[j]*i[e.inds[j],mystart+k-1] 
        end
        
    end
    nothing
end

function mn(n,i,ex,cinds,xyz,inds,xyz1,v,k)
    
    resetn!(xyz,xyz1)
    
    randomize_shape!(n,xyz1,inds)
    
    translate3d!(xyz1,rand(-500.0:500.0),rand(-500.0:500.0),rand(-500.0:500.0))
    
    c=extracellular(ex,xyz1,inds,.3,cinds)
    
    getv(i,c,v,k)   
end


function nete(n,i,ex,num)
     
    mystart=round(Int,size(i,2)/2)
    v=zeros(Float64,length(mystart:size(i,2)),num)
    
    (xyz,inds)=getxyz(n);
    
    cinds=findcind(n)
    
    xyz1=copyn(xyz)
    
    for k=1:num
        mn(n,i,ex,cinds,xyz,inds,xyz1,v,k)
    end
    
    v
end

function extrap(v,t)
   
    v2=zeros(Float64,length(0.0:.025:t))

    spike=find_cspikes(v)

    st=[Array(Int64,0) for i=1:length(spike)]
    
    sigl=size(v,1)
    
    act=zeros(Int64,size(v,2))
    
    @inbounds for i=1:length(v2)       
        for j=1:size(v,2)      
            if act[j]==0
                if  rand()>.9995
                    if sum(j.==spike)!=0
                        push!(st[find(j.==spike)[1]],i)
                    end
                    act[j]+=1
                    v2[i]+=v[1,j]
                end
            else
                v2[i]+=v[act[j],j]
                act[j]+=1
                    if act[j]>=sigl
                        act[j]=0
                    end
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
            @inbounds if v[j,i]<-.03
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
    for i=1:length(n.t)

        for j=1:length(n.stim)
            @inbounds getfield(n.neur,n.stim[j].mtype)[n.stim[j].neur].rhs[n.stim[j].node]+=n.stim[j].Is[i]
        end
        @inbounds main(n.neur)
        
        myi[:,i]=n.neur.N_1[1].i_vm[:]          
    end
    myi       
end
