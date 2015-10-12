
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


#=
Adapted from matlab script by Matt Fig
=#
function tridiagonalize(A::Array{Float64,2})
    
    lngth=size(A,2)
    v=zeros(Float64,lngth)
    I=eye(lngth)
    Anew=zeros(Float64,size(A))
    Aold=zeros(Float64,size(A))
    
    Aold[:]=A[:]

    for jj=1:lngth-2
        v[1:jj] = 0.0
        S=ss(Aold,jj)
        v[jj+1] = sqrt(.5.*(1+abs(Aold[jj+1,jj])/(S+2*eps()) ))
        v[jj+2:lngth] = Aold[jj+2:lngth,jj]*sign(Aold[jj+1,jj]) / (2*v[jj+1]*S+2*eps())
        P = I - 2*v*v'
        Anew = P*Aold*P
        Aold[:]=Anew[:]
        
    end

    Anew[abs(Anew).<5e-14]=0.0

    Anew

end

function ss(A::Array{Float64,2},jj::Int64)
    sqrt(sum(A[(jj+1):end,jj].^2))
end

function initialcon!(neuron::Neuron)

    #fill matrix

    #initial states of channels at nodes

    
end

function main(neuron::Neuron)

    #t=tentry+dt for euler, t=tentry+dt/2 for CN
    
    for i=1:length(neuron.nodes)

        i1=0.0
        i2=0.0
        i=0.0
            
        for k=1:length(neuron.secstack[i].pnode[j].prop)

            #calculate current
            i1=cur_calc(neuron.node[i].prop[k],neuron.node[i],neuron.v_new[i])

            i1=cur_calc(neuron.node[i].prop[k],neuron.node[i],neuron.v_new[i]+.001)
                
            #add -i to rhs
            neuron.rhs[i]-=i1

            #add dv/di to diagonal of A
            neuron.diag[i]+=(i2-i1)
                
        end

            #calculate current entering node from parent and exiting to children
            #add to rhs for that node
    end

    #solve A \ rhs to get delta_v

    #if second order correct, currents are updated?

    #update voltages v_new = delta_v + v_old for euler, v_new  = 2*delta_v + v_old for CN

    #t=tentry+dt for CN

    #find non voltage states (like gate variables for conductances)

    

end
