
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

run!(mynetwork,true);
run!(mynetwork2,true);

facts() do

    for i=1:3
        @fact length(mynetwork.extra[i].v) --> 601
        @fact length(mynetwork.extra[i].coeffs) --> 1
        @fact length(mynetwork.extra[i].coeffs[1].c) --> 726
        @fact length(mynetwork.extra[i].coeffs[1].inds) --> 726
        @fact mynetwork.extra[i].v[end] --> roughly(0.0, atol=1e-4)
    end

    for i=1:3
        @fact length(mynetwork2.extra[i].coeffs) --> 2
        @fact length(mynetwork2.extra[i].v) --> 601
        @fact mynetwork2.extra[i].v[end] --> roughly(0.0, atol=1e-4)
        for j=1:2   
            @fact length(mynetwork2.extra[i].coeffs[j].c) --> 726
            @fact length(mynetwork2.extra[i].coeffs[j].inds) --> 726
        end
    end

end

mynetwork3=Network(deepcopy(myneuron),15.0)
mystim=Stim(5.0,1,1,924,5.0,5.125)
add!(mynetwork3,mystim);
i3=JNeuron.runc(mynetwork3,true);
mye=Extracellular([500.0,125.0,0.0])
mye2=Extracellular([500.0,75.0,0.0])
myv=JNeuron.nete(myneuron,i3,[mye,mye2],100);
(myv3, mys)=JNeuron.extrap(myv,50.0);

facts() do

    @fact size(i3) --> (925,601)
    @fact size(myv) --> (302,100,2)
    @fact size(myv3) --> (2001,2)

end


end
