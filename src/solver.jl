
#=

We start with the cable equation, spatially discretize it, then temporal discretize it. We have two methods of temporally discretizing the cable equation, backward Euler and Crank-Nicolson. Luckily, these are easily related, so we only need to derive the backward Euler cable equation.

The main matrix equation to go from one time step to another for an array of voltage values for all nodes, V_n

 A \ (V_n - sum(i)) = delta_V

Where A is basically the matrix of resistances between nodes, and i is the current from each channel

Only the diagonals of A change from iteration to iteration, by a factor proportional to the conductances of the active components (di/dv).

=#

function fillA!(neuron::Neuron)

    neuron.A=zeros(Float64,length(neuron.secstack),length(neuron.secstack))
    
    for i=1:length(neuron.secstack)

        for j=1:length(neuron.secstack[i].pnode)

            #find parent resistance parent_r=1/(Rp_node)
            parent_r=1/(neuron.secstack[i].pnode[j].ri[1]+neuron.secstack[i].pnode[j].parent.ri[2])
            
            #find children resistance children_r=1/(Rc_node)
            children_r=zeros(Float64,length(neuron.secstack[i].pnode[j].children))

            for k=1:length(children_r)
                children_r[k]=1/(neuron.secstack[i].pnode[j].ri[2]+neuron.secstack[i].pnode[j].children[k].ri[1])
            end
           
            #populate diagonal
            C/neuron.dt+parent_r+sum(children_r)

            #populate parent
            -parent_r

            #populate children (for each child)
            for child_r in children_r
                -child_r
            end

        end
        
          
    end
end

function main()

    #t=tentry+dt/2
    
    #calculate di/dv and i
    
    #add di/dv to A and -i to rhs

    #solve A \ rhs to get delta_v

    #if second order correct, currents are updated?

    #update voltages

    #t=tentry+dt

    #find non voltage states (like gate variables for conductances)

end
