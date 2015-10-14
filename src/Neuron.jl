
#=
Main data container for JNeuron
=#

function add_sec(neuron::Neuron, sec::Section) #Add section to workspace
    push!(neuron.secstack,sec)
end

function set_nsegs(neuron::Neuron,frequency::Float64,d_lambda::Float64)

    for i=1:length(neuron.secstack)
        lambda=lambda_f(frequency, neuron.secstack[i], neuron)
        nseg=round(Int64, neuron.secstack[i].length / (d_lambda * lambda) + .9)

        newstart=length(neuron.nodes)+1
        #create the appropriate number of nodes
        for j=1:nseg

            #add to main node matrix
            push!(neuron.nodes,Node(neuron.secstack[i],j,nseg))

            neuron.nodes[end].ind=length(neuron.nodes)
            
        end

        #map pnode of this section to newly created nodes            
        neuron.secstack[i].pnode=view(neuron.nodes,newstart:(newstart+nseg))
        
    end
    
end

function lambda_f(frequency::Float64,sec::Section,neuron::Neuron)
    
    x1 = 0.0
    d1 = sec.pt3d[1].d
    lam = 0.0

    for i=2:length(sec.pt3d)
        x2=sec.pt3d[i].arc*sec.length
        d2=sec.pt3d[i].d
        lam+= (x2-x1)/sqrt(d1+d2)
        x1=x2
        d1=d2
    end

    lam *= sqrt(2) * 1e-5 * sqrt(4*pi*frequency*sec.Ra*neuron.Cm)

    sec.length/lam
    
end

function connect_nodes(neuron::Neuron)

end
