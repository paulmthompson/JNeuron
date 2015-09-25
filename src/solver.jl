
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

function initialcon!(neuron::Neuron)

    #fill matrix

    #initial states of channels at nodes

    
end

function main(neuron::Neuron)

    #t=tentry+dt for euler, t=tentry+dt/2 for CN
    
    #calculate di/dv and i

    for i=1:length(neuron.secstack)

        for j=1:length(neuron.secstack[i].pnode)

            for k=1:length(neuron.secstack[i].pnode[j].prop)

                #calculate current
                prop_calc(neuron.secstack[i].pnode[j].prop[k],neuron.secstack[i].pnode[j])

                #add to A diagonal

                #calculate dv/di

                #add dv/di to rhs
                
            end
        end
    end

    #solve A \ rhs to get delta_v

    #if second order correct, currents are updated?

    #update voltages v_new = delta_v + v_old for euler, v_new  = 2*delta_v + v_old for CN

    #t=tentry+dt for CN

    #find non voltage states (like gate variables for conductances)

end
