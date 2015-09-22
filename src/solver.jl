
#=

We start with the cable equation, spatially discretize it, then temporal discretize it. We have two methods of temporally discretizing the cable equation, backward Euler and Crank-Nicolson. Luckily, these are easily related, so we only need to derive the backward Euler cable equation.

The main matrix equation to go from one time step to another for an array of voltage values for all nodes, V_n

 A \ (V_n + w_(n+1)) = V_n+1

Where A is basically the matrix of resistances between nodes, and w_n+1 is proportional to the conductance at the next time step.

Only the diagonals of A change from iteration to iteration, by a factor proportional to the conductances of the active components.

=#


function fillA!(neuron::Neuron)

    neuron.A=zeros(Float64,length(neuron.secstack),length(neuron.secstack))

    #find parent resistance parent_r=delta_t/(Rp_node*Cp_node)

    #find children resistance children_r=delta_t/(Rc_node*Cc_node)

    
    for i=1:length(neuron.secstack)
        #populate diagonal
        1+parent_r+sum(children_r)

        #populate parent
        -parent_r

        #populate children (for each child)
        for child_r in children_r
            -child_r
        end
          
    end
end

