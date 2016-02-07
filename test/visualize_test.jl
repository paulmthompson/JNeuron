
module Visualize_Test

using FactCheck, JNeuron

myimport=input(string(dirname(Base.source_path()),"/../examples/data/cell2.asc"));
blank_neuron=instantiate(myimport);
set_nsegs!(blank_neuron);

xy=JNeuron.plot_arrays(blank_neuron)

facts() do

    @fact length(xy) --> 198
    @fact size(xy[1],1) --> size(blank_neuron.secs[1].pt3d,1)
    @fact xy[1][1,1] --> blank_neuron.secs[1].pt3d[1].x
    @fact xy[1][1,2] --> blank_neuron.secs[1].pt3d[1].y
end

blank_neuron2=deepcopy(blank_neuron)

JNeuron.translate3d!(blank_neuron2,1.0,1.0,1.0)

facts() do

    @fact blank_neuron.secs[1].pt3d[1].x-->roughly(blank_neuron2.secs[1].pt3d[1].x-1.0)
    @fact blank_neuron.secs[1].pt3d[1].y-->roughly(blank_neuron2.secs[1].pt3d[1].y-1.0)
    @fact blank_neuron.secs[1].pt3d[1].z-->roughly(blank_neuron2.secs[1].pt3d[1].z-1.0)
end

blank_neuronx=deepcopy(blank_neuron)

JNeuron.rotate3d!(blank_neuronx,pi/2,1)

facts() do

    @fact blank_neuronx.secs[1].pt3d[1].x-->roughly(blank_neuron.secs[1].pt3d[1].x, atol=1e-2)
    @fact blank_neuronx.secs[1].pt3d[1].y-->roughly(-blank_neuron.secs[1].pt3d[1].z, atol=1e-2)
    @fact blank_neuronx.secs[1].pt3d[1].z-->roughly(blank_neuron.secs[1].pt3d[1].y, atol=1e-2)

end

blank_neurony=deepcopy(blank_neuron)

JNeuron.rotate3d!(blank_neurony,pi/2,2)

facts() do

    @fact blank_neurony.secs[1].pt3d[1].x-->roughly(blank_neuron.secs[1].pt3d[1].z, atol=1e-2)
    @fact blank_neurony.secs[1].pt3d[1].y-->roughly(blank_neuron.secs[1].pt3d[1].y, atol=1e-2)
    @fact blank_neurony.secs[1].pt3d[1].z-->roughly(-blank_neuron.secs[1].pt3d[1].x, atol=1e-2)

end

blank_neuronz=deepcopy(blank_neuron)

JNeuron.rotate3d!(blank_neuronz,pi/2,3)

facts() do

    @fact blank_neuronz.secs[1].pt3d[1].x-->roughly(-blank_neuron.secs[1].pt3d[1].y, atol=1e-2)
    @fact blank_neuronz.secs[1].pt3d[1].y-->roughly(blank_neuron.secs[1].pt3d[1].x, atol=1e-2)
    @fact blank_neuronz.secs[1].pt3d[1].z-->roughly(blank_neuron.secs[1].pt3d[1].z, atol=1e-2)

end

blank_neuronr=deepcopy(blank_neuron)

JNeuron.randomize_shape!(blank_neuronr)

facts() do
    
    @fact blank_neuron.secs[10].pt3d[1].x-->not(roughly(blank_neuronr.secs[1].pt3d[1].x))
    @fact blank_neuron.secs[10].pt3d[1].y-->not(roughly(blank_neuronr.secs[1].pt3d[1].y))
    @fact blank_neuron.secs[10].pt3d[1].z-->not(roughly(blank_neuronr.secs[1].pt3d[1].z))
end

end
