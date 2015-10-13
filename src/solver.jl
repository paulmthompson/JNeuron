
#=

We start with the cable equation, spatially discretize it, then temporal discretize it. We have two methods of temporally discretizing the cable equation, backward Euler and Crank-Nicolson. Luckily, these are easily related, so we only need to derive the backward Euler cable equation.

The main matrix equation to go from one time step to another for an array of voltage values for all nodes, V_n

 A \ (V_n - sum(i)) = delta_V

Where A is basically the matrix of resistances between nodes, and i is the current from each channel

Only the diagonals of A change from iteration to iteration, by a factor proportional to the conductances of the active components (di/dv).

=#

function fillA!(neuron::Neuron)

    neuron.A=zeros(Float64,length(neuron.secstack),length(neuron.secstack))
    
    for i=1:length(neuron.nodes)

            #find parent resistance parent_r=1/(Rp_node)
            neuron.nodes[i].parent_r=1/(neuron.nodes[i].ri[1]+neuron.nodes[i].parent.ri[2])
            
            #find children resistance children_r=1/(Rc_node)
            neuron.nodes[i].children_r=zeros(Float64,length(neuron.nodes[i].children))

            for k=1:length(neuron.nodes[i].children_r)
                neuron.nodes[i].children_r[k]=1/(neuron.nodes[i].ri[2]+neuron.nodes[i].children[k].ri[1])
            end
           
            #populate diagonal
            A[i,i]=neuron.C/neuron.dt+neuron.nodes[i].parent_r+sum(neuron.nodes[i].children_r)

            #populate parent
            neuron.A[i,neuron.nodes[i].parent.ind]=-neuron.nodes[i].parent_r

            #populate children (for each child)
            for j=1:length(neuron.nodes[i].children_r)
                neuron.A[i,neuron.nodes[i].children[j].ind]=-neuron.nodes[i].children_r[j]
            end
               
    end

    neuron.mydiag=diagview(A)
    neuron.diag_old=diag(A)
    
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
            
        for k=1:length(neuron.nodes[i].prop)

            #calculate current
            i1=cur_calc(neuron.nodes[i].prop[k],neuron.nodes[i],neuron.v_new[i])

            i1=cur_calc(neuron.nodes[i].prop[k],neuron.nodes[i],neuron.v_new[i]+.001)
                
            #add -i to rhs
            neuron.rhs[i]-=i1

            #add dv/di to diagonal of A
            neuron.diag[i]+=(i2-i1)
                
        end

        #calculate current entering node from parent and exiting to children
        #add to rhs for that node
        
        #parent current
        i_p=(neuron.vold[neuron.nodes[i].parent.ind]-neuron.vold[i])/neuron.nodes[i].parent_r

        #children current
        i_c=0.0
        for j=1:length(neuron.nodes[i].children)
            i_c+=(neuron.v_old[i]-neuron.v_old[neuron.nodes[i].children[j].ind])/neuron.nodes[i].children_r[j]
        end

        neuron.rhs[i]+=i_p
        neuron.rhs[i]+=i_c
            
    end

    #solve A \ rhs to get delta_v
    neuron.delta_v = neuron.A \ neuron.rhs

    #if second order correct, currents are updated?

    #update voltages v_new = delta_v + v_old for euler, v_new  = 2*delta_v + v_old for CN
    neuron.v_old=neuron.delta_v + neuron.v_old

    #t=tentry+dt for CN

    #find non voltage states (like gate variables for conductances)
    
    for i=1:length(neuron.nodes)
        for j=1:length(neuron.nodes[i].prop)
            con_calc(neuron.nodes[i].prop[j],neuron.nodes[i],v[i])
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
