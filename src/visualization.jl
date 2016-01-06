
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

function randomize_shape!(neuron::Neuron)

    for i=1:length(neuron.secstack)

        if (neuron.secstack[i].parent!=length(neuron.secstack))&&(neuron.secstack[i].parent!=0)

            rot_seg(neuron,neuron.secstack[i].parent,i)
            
        end

    end
    
    nothing
    
end

function rot_3d(a,b,c,u,v,w,x,y,z,theta)
    
    xnew=(a*(v*v+w*w)-u*(b*v + c*w - u*x - v*y - w*z))*(1-cosd(theta))+x*cosd(theta)+(-c*v + b*w - w*y + v*z)*sind(theta)
    ynew=(b*(u*u+w*w)-v*(a*u + c*w - u*x - v*y - w*z))*(1-cosd(theta))+y*cosd(theta)+(c*u - a*w + w*x - u*z)*sind(theta)
    znew=(c*(u*u+v*v)-w*(a*u + b*v - u*x - v*y - w*z))*(1-cosd(theta))+z*cosd(theta)+(-b*u + a*v - v*x + u*y)*sind(theta)
    
    (xnew,ynew,znew)
end

function rot_seg(neuron::Neuron,pind::Int64,cind::Int64)

    #last point of previous segment
    a=neuron.secstack[pind].pt3d[end].x
    b=neuron.secstack[pind].pt3d[end].y
    c=neuron.secstack[pind].pt3d[end].z

    (uvw,mag)=pt3d_vec(neuron.secstack[pind].pt3d[end],neuron.secstack[pind].pt3d[end-1])

    theta=randn()*10

    rot_child(neuron,cind,a,b,c,uvw,theta)

    nothing
end

function rot_child(neuron::Neuron,cind::Int64,a::Float64,b::Float64,c::Float64,uvw::Array{Float64,1},theta::Float64)
    
    for sec in neuron.secstack[cind].child

        i=sec.refcount
        rot_child(neuron,i,a,b,c,uvw,theta)
        
    end

    @inbounds for j=1:length(neuron.secstack[cind].pt3d)

        x=neuron.secstack[cind].pt3d[j].x
        y=neuron.secstack[cind].pt3d[j].y
        z=neuron.secstack[cind].pt3d[j].z
            
        xyz=rot_3d(a,b,c,uvw[1],uvw[2],uvw[3],x,y,z,theta)

        neuron.secstack[cind].pt3d[j].x=xyz[1]
        neuron.secstack[cind].pt3d[j].y=xyz[2]
        neuron.secstack[cind].pt3d[j].z=xyz[3]
    end

    nothing
    
end

