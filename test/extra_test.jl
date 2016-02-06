
module Extra_Test

using FactCheck, JNeuron

myimport=input(string(dirname(Base.source_path()),"/../examples/data/cell2.asc"));
blank_neuron=instantiate(myimport);
set_nsegs!(blank_neuron);

myneuron=add(blank_neuron,(HH(),Passive()));
myneuron2=add(blank_neuron,HH());

mynetwork=Network(myneuron,15.0);

electrode1_1=Extracellular([0.0,0.0,0.0]);
add!(mynetwork,electrode1_1);

electrode2_1=Extracellular(JNeuron.Point(),[0.0,0.0,0.0])
add!(mynetwork,electrode2_1);

electrode3_1=Extracellular(JNeuron.Mixed(),[0.0,0.0,0.0])
add!(mynetwork,electrode3_1);


mynetwork2=Network((deepcopy(myneuron),myneuron2),15.0)

electrode1_2=Extracellular([0.0,0.0,0.0]);
add!(mynetwork2,electrode1_2);

electrode2_2=Extracellular(JNeuron.Point(),[0.0,0.0,0.0])
add!(mynetwork2,electrode2_2);

electrode3_2=Extracellular(JNeuron.Mixed(),[0.0,0.0,0.0])
add!(mynetwork2,electrode3_2);

facts() do

    for i=1:3
        @fact length(mynetwork.extra[i].v) --> 601
        @fact length(mynetwork.extra[i].coeffs) --> 1
        @fact length(mynetwork.extra[i].coeffs[1].c) --> 726
        @fact length(mynetwork.extra[i].coeffs[1].inds) --> 726
    end

    
    for i=1:3
        @fact length(mynetwork2.extra[i].coeffs) --> 2
        @fact length(mynetwork2.extra[i].v) --> 601
        for j=1:2   
            @fact length(mynetwork2.extra[i].coeffs[j].c) --> 726
            @fact length(mynetwork2.extra[i].coeffs[j].inds) --> 726
        end
    end

end

end
