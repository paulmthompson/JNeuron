
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
    
    neuron.A=zeros(Float64,nodelen,nodelen)
    neuron.v=zeros(Float64,nodelen)
    neuron.delta_v=zeros(Float64,nodelen)
    neuron.rhs=zeros(Float64,nodelen)
    neuron.i_vm=zeros(Float64,nodelen)
    neuron.divm=zeros(Float64,nodelen)

    #find resistances between nodes

    for i=1:length(neuron.secstack)
        for j=1:length(neuron.secstack[i].pnode)
            #boundary conditions
            if j==1
                neuron.secstack[i].pnode[1].a=1/neuron.secstack[i].pnode[1].ri[1]
	        neuron.secstack[i].pnode[1].b=100/(2*neuron.secstack[i].pnode[1].ri[1]*neuron.secstack[i].pnode[1].area[1])
            else
                if neuron.secstack[i].pnode[j].parent==0
                    neuron.secstack[i].pnode[j].b=1/neuron.secstack[i].pnode[j].ri[1]
                else

                    #parent
                    
                    mya=sum(neuron.secstack[i].pnode[j].area)
                    myr=neuron.secstack[i].pnode[j].ri[1]+neuron.nodes[neuron.secstack[i].pnode[j].parent].ri[2]
                    neuron.secstack[i].pnode[j].b=100/(mya*myr)

                    #child
                    
                    mya=sum(neuron.nodes[neuron.secstack[i].pnode[j].parent].area)
            
                    neuron.secstack[i].pnode[j].a=100/(mya*myr)
                    
                end
                
            end
        end   
    end
        
    for i=1:length(neuron.nodes)

        #populate parent
        if neuron.nodes[i].parent!=0
            neuron.A[i,neuron.nodes[i].parent]=-neuron.nodes[i].b
        end

        #populate children (for each child)
        for j=1:length(neuron.nodes[i].children)
            neuron.A[i,neuron.nodes[i].children[j]]=-neuron.nodes[neuron.nodes[i].children[j]].a
        end

        #populate diagonal
        #don't really understand this
        #Not right for no child case
        if length(neuron.nodes[i].children)==0
            if length(neuron.nodes[neuron.nodes[i].parent].children)>1
                neuron.A[i,i]=.001*neuron.Cm/neuron.dt+2*neuron.nodes[i].b
            else
                neuron.A[i,i]=.001*neuron.Cm/neuron.dt+2*neuron.nodes[i].b+neuron.nodes[neuron.nodes[i].parent].a
            end
            
        elseif (length(neuron.nodes[i].children)>1)&&(neuron.nodes[i].parent!=0)
            if length(neuron.nodes[neuron.nodes[i].parent].children)>1
                neuron.A[i,i]=.001*neuron.Cm/neuron.dt+2*neuron.nodes[i].b
            else
                neuron.A[i,i]=.001*neuron.Cm/neuron.dt+2*neuron.nodes[i].b+neuron.nodes[neuron.nodes[i].parent].a
            end

        else
                
            neuron.A[i,i]=.001*neuron.Cm/neuron.dt+sum(abs(neuron.A[i,:]))
        end
        
        if ext==true
            #find parent resistance parent_r=1/(Rp_node)
            neuron.enodes[i].parent_r=1/(neuron.enodes[i].ri[1]+neuron.enodes[neuron.nodes[i].parent].ri[2])

            for k=1:length(neuron.nodes[i].children_r)
                neuron.enodes[i].children_r[k]=1/(neuron.enodes[i].ri[2]+neuron.enodes[neuron.nodes[i].children[k]].ri[1])
            end

            #populate diagonal
            neuron.A_ext[i,i]=neuron.Cx/neuron.dt+neuron.enodes[i].parent_r+sum(neuron.enodes[i].children_r)

            #because the extracellular node has passive membrane dynamics, di/dvx term is constant
            neuron.A_ext[i,i]+=1/neuron.xg

            #populate parent
            neuron.A_ext[i,neuron.nodes[i].parent]=-neuron.enodes[i].parent_r

            #populate children (for each child)
            for j=1:length(neuron.enodes[i].children_r)
                neuron.A_ext[i,neuron.nodes[i].children[j]]=-neuron.enodes[i].children_r[j]
            end
        
        end
        
    end
  
    neuron.diag=diagview(neuron.A)
    neuron.diag_old=diag(neuron.A)

    if ext==true
        neuron.diag_ext=diagview(neuron.A_ext)
        neuron.diag_ext_old=diag(neuron.A_ext)
    end
       
end

function initialcon!(neuron::Neuron)

    #fill matrix
    fillA!(neuron)

    #initial V?
    neuron.v=neuron.vi.*ones(Float64,length(neuron.v))
    
    #initial states of channels at nodes
    for i=1:length(neuron.nodes)
        for j=1:length(neuron.nodes[i].prop)
            prop_init(neuron.nodes[i].prop[j],neuron.nodes[i],neuron.v[i])
        end

        mykeys=keys(neuron.nodes[i].vars)
        for mykey in mykeys
            neuron.nodes[i].vars[mykey]=myconstants[mykey]
        end
        
    end
    
end

function main(neuron::Neuron)

    #t=tentry+dt for euler, t=tentry+dt/2 for CN
    ext=false

    #reset diagonal
    neuron.diag[:]=neuron.diag_old[:]

    #reset rhs
    neuron.rhs[:]=0.0

    for i=1:length(neuron.nodes)

        neuron.i_vm[i] = 0.0
        neuron.divm[i] = 0.0
            
        for k=1:length(neuron.nodes[i].prop)

            #calculate current
            i1=cur_calc(neuron.nodes[i].prop[k],neuron.nodes[i],neuron.v[i])

            i2=cur_calc(neuron.nodes[i].prop[k],neuron.nodes[i],neuron.v[i]+.001)
            
            neuron.i_vm[i] += i1
            neuron.divm[i] += (i2-i1)/.001
                               
        end

        #add -i to rhs
        neuron.rhs[i] -= neuron.i_vm[i]

        #add dv/di to diagonal
        neuron.diag[i] += neuron.divm[i]

        #calculate current entering node from parent and exiting to children
        #add to rhs for that node
        
        #parent current
        #=
        if neuron.nodes[i].parent!=0
            i_p=(neuron.v[neuron.nodes[i].parent]-neuron.v[i])*neuron.nodes[i].b
        else
            i_p=0.0
        end

        #children current
        i_c=0.0
        for j=1:length(neuron.nodes[i].children)
            i_c+=(neuron.v[i]-neuron.v[neuron.nodes[i].children[j]])*neuron.nodes[neuron.nodes[i].children[j]].a
        end

        neuron.rhs[i]+=i_p
        neuron.rhs[i]+=i_c
        =#

        if neuron.nodes[i].parent != 0
            dv=neuron.v[neuron.nodes[i].parent]-neuron.v[i]
            neuron.rhs[i] += neuron.nodes[i].b*dv
            neuron.rhs[neuron.nodes[i].parent] -= neuron.nodes[i].a*dv
        end
                                 
    end

    #solve A \ rhs to get delta_v
    neuron.delta_v = neuron.A \ neuron.rhs

    if ext==true
        main_ext(neuron)
        neuron.delta_vext = neuron.A_ext \ neuron.rhs_ext
    end
 
    #if second order correct, currents are updated?

    #update voltages v_new = delta_v + v_old for euler, v_new  = 2*delta_v + v_old for CN
    if ext==true
        neuron.vext += neuron.delta_vext
        neuron.v += neuron.delta_v + neuron.delta_vext
    else   
        neuron.v += neuron.delta_v
    end
    
    #t=tentry+dt for CN

    #find non voltage states (like gate variables for conductances)
    
    for i=1:length(neuron.nodes)
        for j=1:length(neuron.nodes[i].prop)
            con_calc(neuron.nodes[i].prop[j],neuron.nodes[i],neuron.v[i],neuron.dt)
        end
    end

    nothing
    
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
