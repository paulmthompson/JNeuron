
#=

Basic 3D reconstruction

=#

function plot_arrays(neuron::Neuron)
    
    x=[zeros(Float64,0) for i=1:length(neuron.secs)]
    y=[zeros(Float64,0) for i=1:length(neuron.secs)]
    
    for i=1:length(neuron.secs), j=1:length(neuron.secs[i].pt3d)
        push!(x[i],neuron.secs[i].pt3d[j].x)
        push!(y[i],neuron.secs[i].pt3d[j].y)
    end

    (x,y)
end

function translate3d!(neuron::Neuron,x::Float64,y::Float64,z::Float64)

    for i=1:length(neuron.secs), j=1:length(neuron.secs[i].pt3d)
            neuron.secs[i].pt3d[j]=add(neuron.secs[i].pt3d[j],x,y,z)     
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

    for i=1:length(neuron.secs)
        for j=1:length(neuron.secs[i].pt3d)

            xyz=rotmat*[neuron.secs[i].pt3d[j].x,neuron.secs[i].pt3d[j].y,neuron.secs[i].pt3d[j].z]

            neuron.secs[i].pt3d[j]=Pt3d(neuron.secs[i].pt3d[j],xyz)          
        end     
    end    
end

function randomize_shape!(neuron::Neuron)

    for i=1:length(neuron.secs)

        if (neuron.secs[i].parent!=length(neuron.secs))&&(neuron.secs[i].parent!=0)

            rot_seg(neuron,neuron.secs[i].parent,i)         
        end
    end   
    nothing   
end

function rot_3d(a,b,c,u,v,w,x,y,z,dcos,dsin)
  
    xnew=(a*(v*v+w*w)-u*(b*v + c*w - u*x - v*y - w*z))*(1-dcos)+x*dcos+(-c*v + b*w - w*y + v*z)*dsin
    ynew=(b*(u*u+w*w)-v*(a*u + c*w - u*x - v*y - w*z))*(1-dcos)+y*dcos+(c*u - a*w + w*x - u*z)*dsin
    znew=(c*(u*u+v*v)-w*(a*u + b*v - u*x - v*y - w*z))*(1-dcos)+z*dcos+(-b*u + a*v - v*x + u*y)*dsin

    (xnew,ynew,znew)
end

function rot_seg(neuron::Neuron,pind::Int64,cind::Int64)

    #last point of previous segment
    a=neuron.secs[pind].pt3d[end].x
    b=neuron.secs[pind].pt3d[end].y
    c=neuron.secs[pind].pt3d[end].z

    (uvw,mag)=pt3d_vec(neuron.secs[pind].pt3d[end],neuron.secs[pind].pt3d[end-1])

    theta=randn()*10

    dcos=cosd(theta)
    dsin=sind(theta)

    rot_child(neuron,cind,a,b,c,uvw,dcos,dsin)
    nothing
end

function rot_child(neuron::Neuron,cind::Int64,a::Float64,b::Float64,c::Float64,uvw::Array{Float64,1},dcos::Float64,dsin::Float64)
    
    for i in neuron.secs[cind].child

        rot_child(neuron,i,a,b,c,uvw,dcos,dsin)      
    end

    @inbounds for j=1:length(neuron.secs[cind].pt3d)

        x=neuron.secs[cind].pt3d[j].x
        y=neuron.secs[cind].pt3d[j].y
        z=neuron.secs[cind].pt3d[j].z
            
        xyz=rot_3d(a,b,c,uvw[1],uvw[2],uvw[3],x,y,z,dcos,dsin)

        neuron.secs[cind].pt3d[j]=Pt3d(xyz[1],xyz[2],xyz[3],neuron.secs[cind].pt3d[j].d,neuron.secs[cind].pt3d[j].arc)

    end
    nothing   
end

