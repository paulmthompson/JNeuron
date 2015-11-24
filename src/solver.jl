
#=

We start with the cable equation, spatially discretize it, then temporal discretize it. We have two methods of temporally discretizing the cable equation, backward Euler and Crank-Nicolson. Luckily, these are easily related, so we only need to derive the backward Euler cable equation.

The main matrix equation to go from one time step to another for an array of voltage values for all nodes, V_n

 A \ (V_n - sum(i)) = delta_V

Where A is basically the matrix of resistances between nodes, and i is the current from each channel

Only the diagonals of A change from iteration to iteration, by a factor proportional to the conductances of the active components (di/dv).

=#

function fillA!(neuron::Neuron)

    ext=false

    nodelen=length(neuron.nodes)
    
    neuron.v=zeros(Float64,nodelen)
    neuron.a=zeros(Float64,nodelen)
    neuron.b=zeros(Float64,nodelen)
    neuron.d=zeros(Float64,nodelen)
    neuron.par=zeros(Int64,nodelen)
    neuron.rhs=zeros(Float64,nodelen)
    neuron.i_vm=zeros(Float64,nodelen)
    neuron.divm=zeros(Float64,nodelen)
    neuron.diag_old=zeros(Float64,nodelen)
    neuron.v1=zeros(Float64,nodelen)
    neuron.i2=zeros(Float64,nodelen)
    
    for i=1:length(neuron.nodes)

        if neuron.nodes[i].parent != 0
            
            r_p=neuron.nodes[neuron.nodes[i].parent].ri[2]+neuron.nodes[i].ri[1]
            area_p=sum(neuron.nodes[neuron.nodes[i].parent].area)
            area_n=sum(neuron.nodes[i].area)

            neuron.nodes[i].b=100/(r_p*area_n)
            neuron.nodes[i].a=100/(r_p*area_p)

            neuron.diag_old[i]=100/(r_p*area_n)+.001*neuron.Cm/neuron.dt
        
            for j in neuron.nodes[i].children
                r_c=neuron.nodes[i].ri[2]+neuron.nodes[j].ri[1]
                neuron.diag_old[i] += 100/(r_c*area_n)
                
            end
            
        else

            neuron.diag_old[i]=.001*neuron.Cm/neuron.dt

            area_n=sum(neuron.nodes[i].area)
            
            for j in neuron.nodes[i].children
                r_c=neuron.nodes[i].ri[2]+neuron.nodes[j].ri[1]
                neuron.diag_old[i] += 100/(r_c*area_n)
                
            end

        end

        neuron.a[i]=neuron.nodes[i].a
        neuron.b[i]=neuron.nodes[i].b
        neuron.par[i]=neuron.nodes[i].parent
        
    end


    for i=1:length(neuron.secstack)
        for j=1:length(neuron.secstack[i].pnode)
            if neuron.secstack[i].pnode[j].internal==true
                push!(neuron.internal_nodes[neuron.secstack[i].mtype],neuron.secstack[i].pnode[j].ind)
            end
        end
    end

    neuron.soma=typeof(neuron.soma)(length(neuron.internal_nodes[1]))
    neuron.axon=typeof(neuron.axon)(length(neuron.internal_nodes[2]))
    neuron.dendrite=typeof(neuron.dendrite)(length(neuron.internal_nodes[3]))
    neuron.apical=typeof(neuron.apical)(length(neuron.internal_nodes[4]))
 
    if ext==true
        neuron.diag_ext=diagview(neuron.A_ext)
        neuron.diag_ext_old=diag(neuron.A_ext)
    end

    nothing
       
end

function initialcon!(neuron::Neuron, vi=-65.0,dt=.025)

    neuron.dt=dt #default
    
    #fill matrix
    fillA!(neuron)

    #initial V
    neuron.v[:]=vi

    for i=1:4
        init!(getfield(neuron,i),neuron.v,neuron.internal_nodes[i])
    end
    
    #initial states of channels at nodes
    for i=1:length(neuron.nodes)
        mykeys=keys(neuron.nodes[i].vars)
        for mykey in mykeys
            neuron.nodes[i].vars[mykey]=myconstants[mykey]
        end
        
    end

    neuron
    
end

function main(neuron::Neuron)

    #t=tentry+dt for euler, t=tentry+dt/2 for CN
    ext=false
    
    for ind=1:4

        cur!(getfield(neuron,ind),neuron.v,neuron.i_vm,neuron.internal_nodes[ind])
        cur!(getfield(neuron,ind),neuron.v1,neuron.i2,neuron.internal_nodes[ind])
        
    end

    @inbounds @simd for i=1:length(neuron.v)
        neuron.divm[i] = (neuron.i2[i]-neuron.i_vm[i])/.001
    end

    rhs_diag!(neuron)
    
    #solve A \ rhs to get delta_v
    hines_solve!(neuron)

    #update voltages v_new = delta_v + v_old for euler, v_new  = 2*delta_v + v_old for CN
    if ext==true
        neuron.vext[:] += neuron.delta_vext
        neuron.v[:] += neuron.delta_v + neuron.delta_vext
    else   
        add_delta!(neuron)
    end

    #find non voltage states (like gate variables for conductances)

    for ind=1:4
        con!(getfield(neuron,ind),neuron.v,neuron.internal_nodes[ind])     
    end
 
    neuron
    
end

function main_ext(neuron::Neuron)
    
     for i=1:length(neuron.nodes)
       
         #calculate current
         iext=1/neuron.xg*(neuron.vext[i]-neuron.ex)
                
         #add iext to rhs
         neuron.rhs_ext[i]+=iext

         #add cm/dt*delta_v + di/vm*delta_v + i(v_m) to rhs (membrane current entering ext space)
         neuron.rhs_ext+=neuron.Cm/neuron.dt*neuron.delta_v[i] + neuron.divm[i]*neuron.delta_v[i] + neuron.i_vm[i]
                
         #calculate current entering node from parent and exiting to children
         #add to rhs for that node
        
         #parent current
         i_p=(neuron.vext[neuron.nodes[i].parent]-neuron.vext[i])/neuron.enodes[i].parent_r

         #children current
         i_c=0.0
         for j=1:length(neuron.nodes[i].children)
             i_c+=(neuron.vext[i]-neuron.vext[neuron.nodes[i].children[j]])/neuron.enodes[i].children_r[j]
         end

         neuron.rhs_ext[i]+=i_p
         neuron.rhs_ext[i]+=i_c
            
     end
end

function add_delta!(neuron::Neuron)
    @fastmath @inbounds @simd for i=1:length(neuron.v)
        neuron.v[i] += neuron.rhs[i]
        neuron.rhs[i] = 0.0
        neuron.v1[i]=neuron.v[i]+.001
        neuron.i2[i]=0.0
        neuron.i_vm[i]=0.0
    end
    nothing
end

function rhs_diag!(neuron::Neuron)

    dv=0.0
    
    @fastmath @inbounds @simd for i=1:(length(neuron.v)-2)
        dv=neuron.v[neuron.par[i]]-neuron.v[i]
        neuron.rhs[i] += neuron.b[i]*dv
        neuron.rhs[neuron.par[i]] -= neuron.a[i]*dv
        neuron.d[i] = neuron.diag_old[i] + neuron.divm[i]
        neuron.rhs[i] -= neuron.i_vm[i]
    end

    i=length(neuron.v)-1
    
    neuron.d[i] = neuron.diag_old[i] + neuron.divm[i]
    neuron.rhs[i] -= neuron.i_vm[i]

    i=length(neuron.v)

    dv=neuron.v[neuron.par[i]]-neuron.v[i]
    neuron.rhs[i] += neuron.b[i]*dv
    neuron.rhs[neuron.par[i]] -= neuron.a[i]*dv
    neuron.d[i] = neuron.diag_old[i] + neuron.divm[i]
    neuron.rhs[i] -= neuron.i_vm[i]
    
    nothing
    
end

function hines_solve!(neuron::Neuron)

    p=0.0

    @fastmath @inbounds @simd for i=(length(neuron.d)-3):-1:1

	p = -neuron.a[i] / neuron.d[i]
	neuron.d[neuron.par[i]] -= p * -neuron.b[i]
	neuron.rhs[neuron.par[i]] -= p * neuron.rhs[i]

    end

    i=length(neuron.d)-2

    p = -neuron.a[i] / neuron.d[i]
    neuron.d[neuron.par[i]] -= p * -neuron.b[i]
    neuron.rhs[neuron.par[i]] -= p * neuron.rhs[i]

    i=length(neuron.d)-2

    p = -neuron.a[i] / neuron.d[i]
    neuron.d[neuron.par[i]] -= p * -neuron.b[i]
    neuron.rhs[neuron.par[i]] -= p * neuron.rhs[i]

    neuron.rhs[length(neuron.d)-1] /= neuron.d[length(neuron.d)-1]
	
    i=length(neuron.d)

    neuron.rhs[i] -= -neuron.b[i] * neuron.rhs[neuron.par[i]]
    neuron.rhs[i] /= neuron.d[i]

    @fastmath @inbounds @simd for i=1:(length(neuron.d)-2)
	neuron.rhs[i] -= -neuron.b[i] * neuron.rhs[neuron.par[i]]
	neuron.rhs[i] /= neuron.d[i]
    end

nothing

end

