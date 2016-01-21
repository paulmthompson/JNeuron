#=

Extracellular 
(x,y,z) extracellular position
sigma: extracellular conductivity
=#

#Point Source

function extracellular{T<:Point}(extra::Extracellular{T},neuron::Neuron,sigma::Float64)

    coeffs=zeros(Float64,0)
    inds=zeros(Int64,0)

    for i=1:length(neuron.secs)
        for j in neuron.secs[i].pnode
            if neuron.nodes[j].internal == true

                push!(inds,j)
                mypt3d=neuron.secs[i].pt3d[neuron.nodes[j].pt3d]
                push!(coeffs,(1/(4*pi*sigma))*point_coeffs(mypt3d,extra.xyz))
                
            end
        end
    end
      
    Extra_coeffs(coeffs,inds)      
end

#Line Source

function extracellular{T<:Line}(extra::Extracellular{T},neuron::Neuron,sigma::Float64)

    coeffs=zeros(Float64,0)
    inds=zeros(Int64,0)

    for i=1:length(neuron.secs)
        for j in neuron.secs[i].pnode
            if neuron.nodes[j].internal == true

                push!(inds,j)
                mypt3d=neuron.secs[i].pt3d[neuron.nodes[j].pt3d]
                push!(coeffs,(1/(4*pi*sigma))*line_coeffs(mypt3d,extra.xyz))
                
            end
        end
    end 
    Extra_coeffs(coeffs,inds)  
end

#Mixed Source (Soma as a point, everything else as line)

function extracellular{T<:Mixed}(extra::Extracellular{T},neuron::Neuron,sigma::Float64)

    coeffs=zeros(Float64,0)
    inds=zeros(Int64,0)

    for i=1:length(neuron.secs)
        for j in neuron.secs[i].pnode
            if neuron.nodes[j].internal == true

                push!(inds,j)
                mypt3d=neuron.secs[i].pt3d[neuron.nodes[j].pt3d]

                if neuron.nodes[j].parent==0
                    push!(coeffs,(1/(4*pi*sigma))*point_coeffs(mypt3d,extra.xyz))
                else
                    push!(coeffs,(1/(4*pi*sigma))*line_coeffs(mypt3d,extra.xyz))

                end            
            end
        end
    end  
    Extra_coeffs(coeffs,inds)  
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
        coeffs[count]=(1/(4*pi*sigma))*line_coeffs(xyz,extra.xyz,i)
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

    a=sqrt(sum(d1.^2))-hn
    b=sqrt(sum(d2.^2))-ln

    1/delta_s*log(abs(a/b))  
end

function line_coeffs(xyz::Array{Float64,2},xyze::Array{Float64,1},i::Int64)

    (ds_x,ds_y,ds_z,delta_s)=pt3d_vec(xyz[i,1],xyz[i,2],xyz[i,3],xyz[i+1,1],xyz[i+1,2],xyz[i+1,3])

    (d1_x,d1_y,d1_z)=vect(xyz[i,1],xyz[i,2],xyz[i,3],xyze[1],xyze[2],xyze[3])
    (d2_x,d2_y,d2_z)=vect(xyz[i+1,1],xyz[i+1,2],xyz[i+1,3],xyze[1],xyze[2],xyze[3])

    ln=dot(ds_x,ds_y,ds_z,d1_x,d1_y,d1_z)
    hn=dot(ds_x,ds_y,ds_z,d2_x,d2_y,d2_z)

     a=sqrt(d1_x*d1_x+d1_y*d1_y+d1_z*d1_z)-hn
     b=sqrt(d2_x*d2_x+d2_y*d2_y+d2_z*d2_z)-ln

    1/delta_s*log(abs(a/b))
    
end

function point_coeffs(p::Array{Pt3d,1},xyz::Array{Float64,1})

    middle=round(length(p)/2)
            
    dist1=vect(p[middle],xyz)

    1/(dist1)
            
end

dist(p1::Pt3d, p2::Pt3d)=sqrt((p1.x-p2.x)^2+(p1.y-p2.y)^2+(p1.z-p2.z)^2)

dist(p1::Pt3d,xyz::Array{Float64,1})=sqrt((p1.x-xyz[1])^2+(p1.y-xyz[2])^2+(p1.z-xyz[3])^2)

function dist(x1,y1,z1,x2,y2,z2)
    sqrt((x1-x2)^2+(y1-y2)^2+(z1-z2)^2)
end

vect(p::Pt3d,xyz::Array{Float64,1})=[xyz[1]-p.x,xyz[2]-p.y,xyz[3]-p.z]

vect(p1::Pt3d,p2::Pt3d)=[p2.x-p1.x,p2.y-p1.y,p2.z-p1.z]

vect(x1,y1,z1,x2,y2,z2)=(x2-x1,y2-y1,z2-z1)

pt3d_vec(p1::Pt3d,p2::Pt3d)=(n=vect(p1,p2); mag=dist(p1,p2); n /= mag; (n, mag))

function pt3d_vec(x1::Float64,y1::Float64,z1::Float64,x2::Float64,y2::Float64,z2::Float64)

    (vx,vy,vz)=vect(x1,y1,z1,x2,y2,z2)

    mag=dist(x1,y1,z1,x2,y2,z2)

    vx /= mag
    vy /= mag
    vz /= mag

    (vx,vy,vz,mag)
    
end

function dot(x1,y1,z1,x2,y2,z2)
    sum=x1*x2
    sum+=y1*y2
    sum+=z1*z2   
end
