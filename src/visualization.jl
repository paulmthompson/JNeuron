
#=

Basic 3D reconstruction

=#

function translate3d!(neuron::Neuron,x::Float64,y::Float64,z::Float64)

    for i=1:length(neuron.secstack)

        for j=1:length(neuron.secstack[i].pt3d)

            neuron.secstack[i].pt3d[j].x+=x
            neuron.secstack[i].pt3d[j].y+=y
            neuron.secstack[i].pt3d[j].z+=z
            
        end
        
    end
    
end

function rotate3d!(neuron::Neuron,theta::Float64,myaxis::Int64)

    rotmat=zeros(Float64,3,3)
    
    if myaxis==1

        rotmat[1,1]=1.0
        rotmat[2,2]=cos(theta)
        rotmat[2,3]=-sin(theta)
        rotmat[3,2]=sin(theta)
        rotmat[3,3]=cos(theta)
        
    elseif myaxis==2

        rotmat[1,1]=cos(theta)
        rotmat[1,3]=sin(theta)
        rotmat[2,2]=1.0
        rotmat[3,1]=-sin(theta)
        rotmat[3,3]=cos(theta)

    else

        rotmat[1,1]=cos(theta)
        rotmat[1,2]=-sin(theta)
        rotmat[2,1]=sin(theta)
        rotmat[2,2]=cos(theta)
        rotmat[3,3]=1.0
        
    end

    for i=1:length(neuron.secstack)

        for j=1:length(neuron.secstack[i].pt3d)

            xyz=rotmat*[neuron.secstack[i].pt3d[j].x,neuron.secstack[i].pt3d[j].y,neuron.secstack[i].pt3d[j].z]
                  
            neuron.secstack[i].pt3d[j].x=xyz[1]
            neuron.secstack[i].pt3d[j].y=xyz[2]
            neuron.secstack[i].pt3d[j].z=xyz[3]
            
        end
        
    end
    
end

