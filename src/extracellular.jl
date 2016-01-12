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
      
    (coeffs,inds)      
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
    (coeffs,inds)  
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
    (coeffs,inds)  
end

function line_coeffs(pt3d::Array{Pt3d,1},xyz::Array{Float64,1})
    
    (unit_ds, delta_s)=pt3d_vec(pt3d[end],pt3d[1])

    d1=vect(pt3d[1],xyz)
    d2=vect(pt3d[end],xyz)

    ln=dot(unit_ds,d1)
    hn=dot(unit_ds,d2)

    a=sqrt(sum(d1.^2))-hn
    b=sqrt(sum(d2.^2))-ln

    1/delta_s*log(abs(a/b))  
end

function point_coeffs(p::Array{Pt3d,1},xyz::Array{Float64,1})

    middle=round(length(p)/2)
            
    dist1=vect(p[middle],xyz)

    1/(dist1)
            
end

dist(p1::Pt3d, p2::Pt3d)=sqrt((p1.x-p2.x)^2+(p1.y-p2.y)^2+(p1.z-p2.z)^2)

dist(p1::Pt3d,xyz::Array{Float64,1})=sqrt((p1.x-xyz[1])^2+(p1.y-xyz[2])^2+(p1.z-xyz[3])^2)

vect(p::Pt3d,xyz::Array{Float64,1})=[xyz[1]-p.x,xyz[2]-p.y,xyz[3]-p.z]

vect(p1::Pt3d,p2::Pt3d)=[p2.x-p1.x,p2.y-p1.y,p2.z-p1.z]

pt3d_vec(p1::Pt3d,p2::Pt3d)=(n=vect(p1,p2); mag=dist(p1,p2); n /= mag; (n, mag))
