
module Neuron_Test

using FactCheck, JNeuron

facts() do

    myimport=input(string(dirname(Base.source_path()),"/../examples/data/cell2.asc"));
    blank_neuron=instantiate(myimport);
    @fact typeof(blank_neuron) --> JNeuron.Neuron_0

    set_nsegs!(blank_neuron);

    @fact length(blank_neuron.nodes) --> 925
    
    @fact Int64[length(blank_neuron.internal_nodes[i]) for i=1:4] --> [1,179,193,353]

    #Soma
    @fact blank_neuron.nodes[924].parent --> 0
    @fact blank_neuron.nodes[924].children --> [1,217,229,317,321,397,417,453,473,923,925]
    @fact blank_neuron.secs[198].pnode --> 923:925

    myneuron1=add(blank_neuron,HH());
    @fact typeof(myneuron1) --> JNeuron.Neuron_1
    myneuron2=add(blank_neuron,(HH(),Passive()));
    @fact typeof(myneuron2) --> JNeuron.Neuron_2
    myneuron3=add(blank_neuron,(HH(),Passive()),(HH(),Passive()),Passive(),Passive());
    @fact typeof(myneuron3) --> JNeuron.Neuron_3

    @fact size(myneuron1.soma.p) --> (16,1)
    @fact size(myneuron1.axon.p) --> (16,179)
    @fact size(myneuron1.dendrite.p) --> (16,193)
    @fact size(myneuron1.apical.p) --> (16,353)

    @fact size(myneuron2.soma.p) --> (17,1)
    @fact size(myneuron2.axon.p) --> (17,179)
    @fact size(myneuron2.dendrite.p) --> (17,193)
    @fact size(myneuron2.apical.p) --> (17,353)

    @fact size(myneuron3.soma.p) --> (17,1)
    @fact size(myneuron3.axon.p) --> (17,179)
    @fact size(myneuron3.dendrite.p) --> (3,193)
    @fact size(myneuron3.apical.p) --> (3,353)


    #mynetwork=Network((myneuron1,myneuron2,myneuron3),100.0);

    #run!(mynetwork,true);
    #run!(mynetwork,false);

    #@fact mynetwork.neur.N_1[1].v[1] --> roughly(-64.9604)
    #@fact mynetwork.neur.N_2[1].v[1] --> roughly(-67.5302)
    
end

end
