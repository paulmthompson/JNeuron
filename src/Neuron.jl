
#=
Main data container for JNeuron
=#

function add_sec(neuron::Neuron, sec::Section) #Add section to workspace
    push!(neuron.secstack,sec)
end

function set_nsegs(neuron::Neuron,frequency::Float64,d_lambda::Float64)

    nodesec=zeros(Int64,length(neuron.secstack))
    nseglist=zeros(Int64,length(neuron.secstack)+1)
    nodesec[1]=1

    #find total number of segments
    for i=1:length(neuron.secstack)
        
        lambda=lambda_f(frequency, neuron.secstack[i], neuron)
        nseg=round(Int64, neuron.secstack[i].length / (d_lambda * lambda) + .9)

        nodesec[i+1]=nodesec[i]+nseg
        
        nseglist[i]=nseg
        
    end
    
    for i=1:length(neuron.secstack)

        if neuron.secstack[i].mtype>1

            for j=1:nseglist[i]

                if j==1
                    parent=
                    children=Int64[length(neuron.nodes)+2]
                elseif j==nseglist[i]
                    parent=length(neuron.nodes)
                    children=Int64[nodesec[neuron.secstack[i].child[j].refcount] for j=1:length(neuron.secstack[i].child)]
                else
                    parent=length(neuron.nodes)
                    children=Int64[length(neuron.nodes)+2]
                end

                (area, ri) = r_a_calc(neuron.secstack[i],j,nseglist[i])
                myvars=Dict{ASCIIString,Float64}()
                push!(neuron.nodes,Node(length(neuron.nodes)+1,myvars,area,ri,0.0,zeros(Float64,length(children)),parent,children,Array(Prop,0))) 
                
            end

        else

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
