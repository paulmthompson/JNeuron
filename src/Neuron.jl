
#=
Main data container for JNeuron
=#

function add_sec(neuron::Neuron, sec::Section) #Add section to workspace
    push!(neuron.secstack,sec)
end

function set_nsegs!(neuron::Neuron,frequency=100.0,d_lambda=.1)

    nodesec=zeros(Int64,length(neuron.secstack)+1)
    nseglist=zeros(Int64,length(neuron.secstack))
    nodesec[1]=1
    parents=zeros(Int64,length(neuron.secstack))
    internodes=0
    first=0
    last=0
    

    #find total number of segments
    for i=1:length(neuron.secstack)
        
        lambda=lambda_f(frequency, neuron.secstack[i], neuron)
        nseg=floor(Int64, (neuron.secstack[i].length / (d_lambda * lambda) + .9) /2 )*2+1

        nodesec[i+1]=nodesec[i]+nseg
        
        nseglist[i]=nseg

        for j=1:length(neuron.secstack[i].child)
            parents[neuron.secstack[i].child[j].refcount]=i
        end
                        
    end
    
    for i=1:length(neuron.secstack)

        for j=1:nseglist[i]

            if neuron.secstack[i].mtype>1
                
                if j==1
                    parent=nodesec[parents[i]+1]+parents[i]-1
                    if nseglist[i]==1
                        if length(neuron.secstack[i].child)==0
                            children=Array(Int64,0)
                        else
                            children=Int64[nodesec[neuron.secstack[i].child[k].refcount]+neuron.secstack[i].child[k].refcount-1 for k=1:length(neuron.secstack[i].child)]
                        end
                        
                    else
                        children=Int64[length(neuron.nodes)+2]
                    end
                elseif j==nseglist[i]
                    parent=length(neuron.nodes)
                    children=Int64[nodesec[neuron.secstack[i].child[k].refcount]+neuron.secstack[i].child[k].refcount-1 for k=1:length(neuron.secstack[i].child)]
                else
                    parent=length(neuron.nodes)
                    children=Int64[length(neuron.nodes)+2]
                end

            else
                
                middle=round(Int64,(1+nseglist[i])/2)

                if j<middle
                    parent=length(neuron.nodes)+1
                    children=Array(Int64,0)
                elseif length(neuron.nodes)==(nodesec[end]-1)
                    parent=length(neuron.nodes)
                    children=Array(Int64,0)
                elseif j>middle
                    parent=length(neuron.nodes)
                    children=Int64[length(neuron.nodes)+2]
                else
                    parent=0
                    children=Int64[nodesec[neuron.secstack[i].child[j].refcount] for j=1:length(neuron.secstack[i].child)]
                end
            end         
            
            (area, ri,mypt3d) = r_a_calc(neuron.secstack[i],j,nseglist[i])
            myvars=Dict{ASCIIString,Float64}()

            if (j==1)&&(j==nseglist[i])

                first=length(neuron.nodes)+1
                last=length(neuron.nodes)+2
                neuron.internal_nodes+=1
                internodes+=1

                #create regular node
                push!(neuron.nodes,Node(length(neuron.nodes)+1,myvars,area,ri,0.0,0.0,parent,[length(neuron.nodes)+2],true,mypt3d,Prop()))

                #create internode at end of section
                push!(neuron.nodes,Node(length(neuron.nodes)+1,myvars,[100.0],[0.0,0.0],0.0,0.0,length(neuron.nodes),children,false,[mypt3d[end]],Prop()))
                
            elseif j==1

                first=length(neuron.nodes)+1
                neuron.internal_nodes+=1
                
                #create regular node
                push!(neuron.nodes,Node(length(neuron.nodes)+1,myvars,area,ri,0.0,0.0,parent,children,true,mypt3d,Prop()))
                            
            elseif j==nseglist[i]

                last=length(neuron.nodes)+2
                neuron.internal_nodes+=1
                internodes+=1
                
                #create regular node
                push!(neuron.nodes,Node(length(neuron.nodes)+1,myvars,area,ri,0.0,0.0,parent,[length(neuron.nodes)+2],true,mypt3d,Prop()))

                #create internode at end of section
                push!(neuron.nodes,Node(length(neuron.nodes)+1,myvars,[100.0],[0.0,0.0],0.0,0.0,length(neuron.nodes),children,false,[mypt3d[end]],Prop()))

            else

                neuron.internal_nodes+=1
                
                push!(neuron.nodes,Node(length(neuron.nodes)+1,myvars,area,ri,0.0,0.0,parent,children,true,mypt3d,Prop()))

            end
                    
        end

        #map pnode of this section to newly created nodes            
        neuron.secstack[i].pnode=sub(neuron.nodes,first:last)
               
    end

    nothing
      
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

