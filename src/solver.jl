
#=

We start with the cable equation, spatially discretize it, then temporal discretize it. We have two methods of temporally discretizing the cable equation, backward Euler and Crank-Nicolson. Luckily, these are easily related, so we only need to derive the backward Euler cable equation.

The main matrix equation to go from one time step to another for an array of voltage values for all nodes, V_n

 A \ (V_n - sum(i)) = delta_V

Where A is basically the matrix of resistances between nodes, and i is the current from each channel

Only the diagonals of A change from iteration to iteration, by a factor proportional to the conductances of the active components (di/dv).

=#

function fillA!(neuron::Neuron)

    for i=1:length(neuron.nodes)

        if neuron.nodes[i].parent != 0
            
            r_p=neuron.nodes[neuron.nodes[i].parent].ri.r+neuron.nodes[i].ri.l
            area_p=neuron.nodes[neuron.nodes[i].parent].area.t
            area_n=neuron.nodes[i].area.t

            neuron.b[i]=100/(r_p*area_n)
            neuron.a[i]=100/(r_p*area_p)

            neuron.diag_old[i]=100/(r_p*area_n)+.001*neuron.Cm/neuron.dt
        
            for j in neuron.nodes[i].children
                r_c=neuron.nodes[i].ri.r+neuron.nodes[j].ri.l
                neuron.diag_old[i] += 100/(r_c*area_n)              
            end         
        else
            
            neuron.diag_old[i]=.001*neuron.Cm/neuron.dt
            area_n=neuron.nodes[i].area.t
            
            for j in neuron.nodes[i].children
                r_c=neuron.nodes[i].ri.r+neuron.nodes[j].ri.l
                neuron.diag_old[i] += 100/(r_c*area_n)             
            end
        end

        neuron.par[i]=neuron.nodes[i].parent      
    end
    nothing       
end

function initialcon!(neuron::Neuron, vi=-65.0,dt=.025)

    neuron.dt=dt

    neuron.v[:]=vi

    for i=1:4
        init!(getfield(neuron,i),neuron.v,neuron.internal_nodes[i])
    end   
end

function main(neuron::Neuron)

    #reset membrane current
    @inbounds @simd for i=1:length(neuron.v)
        neuron.i_vm[i] = 0.0
    end
    
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
    add_delta!(neuron)

    #find non voltage states (like gate variables for conductances)
    for ind=1:4
        con!(getfield(neuron,ind),neuron.v,neuron.internal_nodes[ind])     
    end
    nothing
end

function add_delta!(neuron::Neuron)
    @fastmath @inbounds @simd for i=1:length(neuron.v)
        neuron.v[i] += neuron.rhs[i]
        neuron.rhs[i] = 0.0
        neuron.v1[i]=neuron.v[i]+.001
        neuron.i2[i]=0.0
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

    @fastmath @inbounds for i=(length(neuron.d)-3):-1:1

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

    @fastmath @inbounds for i=1:(length(neuron.d)-2)
	neuron.rhs[i] -= -neuron.b[i] * neuron.rhs[neuron.par[i]]
	neuron.rhs[i] /= neuron.d[i]
    end
    nothing
end

