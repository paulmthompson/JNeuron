
module Import_Test

using FactCheck, JNeuron

facts() do

    myimport=input(string(dirname(Base.source_path()),"/../examples/data/cell2.asc"));
    @fact length(myimport.sections) --> 198

    blank_neuron=instantiate(myimport);
    @fact length(blank_neuron.secs) --> 198

    blank_neuron2=instantiate(string(dirname(Base.source_path()),"/../examples/data/cell2.asc"));
    @fact length(blank_neuron2.secs) --> 198

    a=zeros(Int64,4)
    for i=1:length(blank_neuron.secs)
        a[blank_neuron.secs[i].mtype]+=1
    end

    @fact a --> [1, 37, 63, 97]
    
    @fact blank_neuron.secs[198].child --> [1, 38, 43, 62, 63, 80, 87, 96, 101]

end


end
