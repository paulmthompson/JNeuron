#=

Extracellular 
(x,y,z) extracellular position
sigma: extracellular conductivity
=#

#Point Source

function extracellular{T<:Point}(extra::Extracellular{T},neuron::Neuron,sigma::Float64)

    coeffs=zeros(Float64,0)

    for i=1:length(neuron.nodes)

        if neuron.nodes[i].internal == true

            push!(coeffs,(1/(4*pi*sigma))*point_coeffs(neuron.nodes[i].pt3d,extra.xyz))
            
        else
        end

    end
    
    coeffs  
    
end

#Line Source

function extracellular{T<:Line}(extra::Extracellular{T},neuron::Neuron,sigma::Float64)

    coeffs=zeros(Float64,0)
    
    for i=1:length(neuron.nodes)

        if neuron.nodes[i].internal == true
            
            push!(coeffs,(1/(4*pi*sigma))*line_coeffs(neuron.nodes[i].pt3d,extra.xyz))
            
        else
        end

    end
    
    coeffs  

end

#Mixed Source (Soma as a point, everything else as line)

function extracellular{T<:Mixed}(extra::Extracellular{T},neuron::Neuron,sigma::Float64)

    coeffs=zeros(Float64,0)
    
    for i=1:length(neuron.nodes)

        if neuron.nodes[i].internal == true

            if neuron.nodes[i].parent==0
                
                push!(coeffs,(1/(4*pi*sigma))*point_coeffs(neuron.nodes[i].pt3d,extra.xyz))
            else
                push!(coeffs,(1/(4*pi*sigma))*line_coeffs(neuron.nodes[i].pt3d,extra.xyz))
            end
 
        else
        end

    end
    
    coeffs  

end

function line_coeffs(pt3d::Array{Pt3d,1},xyz::Array{Float64,1})
    
    (unit_ds, delta_s)=pt3d_vec(pt3d[end],pt3d[1])

    dist1=pt3d_xyz_vec(pt3d[1],xyz)
    dist2=pt3d_xyz_vec(pt3d[end],xyz)

    ln=dot(unit_ds,dist1)
    hn=dot(unit_ds,dist2)

    a=sqrt(sum(dist1.^2))-hn
    b=sqrt(sum(dist2.^2))-ln

    1/delta_s*log(abs(a/b))
    
end

function point_coeffs(pt3d::Array{Pt3d,1},xyz::Array{Float64,1})

    middle=round(length(pt3d)/2)
            
    dist1=pt3d_xyz_vec(pt3d[middle],xyz)

    1/(dist1)
            
    
end

#distance between two 3d points
function dist(pt3d1::Pt3d, pt3d2::Pt3d)
    dist=sqrt((pt3d1.x-pt3d2.x)^2+(pt3d1.y-pt3d2.y)^2+(pt3d1.z-pt3d2.z)^2)
end

function dist(pt3d1::Pt3d,xyz::Array{Float64,1})
    dist=sqrt((pt3d1.x-xyz[1])^2+(pt3d1.y-xyz[2])^2+(pt3d1.z-xyz[3])^2)
end

function pt3d_xyz_vec(pt3d1::Pt3d,xyz::Array{Float64,1})
    mynorm=zeros(Float64,3)
    mynorm[1]=(xyz[1]-pt3d1.x)
    mynorm[2]=(xyz[2]-pt3d1.y)
    mynorm[3]=(xyz[3]-pt3d1.z)

    mynorm
end

function pt3d_vec(pt3d1::Pt3d,pt3d2::Pt3d)

    mynorm=zeros(Float64,3)
    mynorm[1]=(pt3d2.x-pt3d1.x)
    mynorm[2]=(pt3d2.y-pt3d1.y)
    mynorm[3]=(pt3d2.z-pt3d1.z)

    mag=dist(pt3d1,pt3d2)

    mynorm /= mag

    (mynorm, mag)
end
