
#################
Getting Started
#################

*************
Installation
*************

*************
Foundations
*************

Julia is based on types, or data containers for specific purposes. Consequently, JNeuron defines types for a network (or everything you are modeling) and everything that you can put into a network for simulation.

==========
Networks
==========

The main type of JNeuron is the Network type, which is actually a relatively simple data container:

.. code-block:: julia

	type Network{T <: NeuronPool}
    		neur::T #Array of Neurons for simulation
    		t::FloatRange{Float64} #time range for simulation
    		extra::Array{Extracellular,1} #Extracellular Recording
    		#extracellular Stimulation
    		intra::Array{Intracellular,1} #Intracellular Recording
    		stim::Array{Stim,1} #Intracellular Stimulation
	end

As you can see, the Network type contains fields for any Neurons, extracellular electrodes or intracellular electrodes that you made add, as well as the current time current to all of these components. 

Because extracelluar and intracellular electrodes depend on the neurons that are present, Networks should first be initialized by 1) the neurons that are to be present, and 2) the simulation time (this can be chagned later. 

.. code-block:: julia

	runtime=15.0 #in ms
	mynetwork=Network(myneuron,runtime)

Networks can also be initalized by groups of neurons

.. code-block:: julia

	mynetwork=Network([myneuron1,myneuron2],runtime)

This is also the step where you specify whether to take advantage of parallel processing.

.. code-block:: julia

	mynetwork=Network([myneuron1,myneuron2],runtime, par=true)


------------------
Adding Components
------------------

Anything that can be in a network is placed there with the add! function. 

.. code-block:: julia
	
	add!(mynetwork,myelectrode)

--------------------
Running Simulations
--------------------

Once you have specified the components of the network, it has everything it needs to carry out the simulation.

.. code-block:: julia
	
	run!(mynetwork)

Any components that you specified to calculate (like extracellular voltage above) will be stored in the data container at the end of the run for you to access:

.. code-block:: julia
	
	mynetwork.extra[1].v #voltage calculated on extracellular electrode


******************
Network Components
******************

The components available to place into a network for simulation are 1) neurons, 2) intracellular stimulation electrodes, 3) intracellular recording, 4) extracellular recording, and 5) extracellular stimulation.

===========
NeuronPool
===========

Subtypes of the Neuron abstract type are data containers with all of the variables necessary to solve the cable equation at each iteration. It can be created from 3D morphology data, and contain a variety of ion channels. Neurons can not only have different shapes and sizes, but also different channels types, distributed differently among the soma, axon and dendrites. To maintain performance, under the hood JNeuron will be creating concrete types for each combinations of ion channels that you use, called Neuron_1, Neuron_2, Neuron_3 etc. The user doesn't have to worry about this, and can create neurons as follows:

---------------------
Loading 3D Structure
---------------------

Neurons can be created from 3D reconstructions from imaging. Right now only Neurolucida file formats (.asc) are supported, but with others soon to come! The 3D data is first parsed by JNeuron to fill a  Import3D type that contains information about the 3D structure. Then this type can be used to generate the Neuron type that holds all of the data necessary to solve the cable equation.

.. code-block:: julia

	filepath="/path/to/file/cell.asc"
	myimport=input(filepath)
	blank_neuron=instantiate(myimport)

If you don't need the 3D information in Import3D, you can just call instantiate with a filepath ending in .asc

.. code-block:: julia

	filepath="/path/to/file/cell.asc"
	blank_neuron=instantiate(filepath)

----------------
Discretization
----------------

In JNeuron, the cable equation is discretized so that voltages are solved at certain points, or nodes (see Theory section). We have to decide for a given neuron shape, what is the appropriate number of nodes? By default, JNeuron uses the lambda-rule (see Theory). 

.. code-block:: julia

	set_nsegs!(blank_neuron)

----------------
Adding Channels
----------------

Neurons can have a variety of ion channels. Many from the literature are already defined in JNeuron (see Channels). Neurons are considered to have 4 types of sections: 1) cell body, 2) axon, 3) basal dendrites and 4) apical dendrites. Different types of ion channels can be present in each of these 4 sections. If the same channels are present everywhere in the neuron, they can easily be added my calling the add method with an instance of the channel type. Multiple channels are added as a tuple of channel types:

.. code-block:: julia

	myneuron1=add(blank_neuron,HH()) #add hodgkin huxley channels

	myneuron2=add(blank_neuron,(HH(),Passive())) #add hodgkin huxley and Passive channels

If we instead want Hodgkin Huxley and Passive channels in the soma and axon, but only passive channels in the basal and apical dendrites, we can call the add function with a channel/tuple argument for each location:

.. code-block:: julia
	
	# 1) cell body, 2) axon 3) basal 4) apical
	myneuron3=add(blank_neuron,(HH(),Passive()),(HH(),Passive()),Passive(),Passive());

Notice that in this step the input to the add method is the "blank_neuron" which has no channel types. When we call the add function, JNeuron is actually performing several metaprogramming steps under the hood, by putting together all of the methods that get called for each ion channel, and generating a new type unique for that neuron.

---------------------------------------------------
Constructing A Network with different neuron types
---------------------------------------------------

Above, we have 3 different neurons, which all have different combinations of neuron channels in different places. We would add these to a network as a tuple:

.. code-block:: julia

	simulation_time=100.0 #ms
	mynetwork1=Network((myneuron1,myneuron2,myneuron3),simulation_time)

The first field of the network, neur, is a NeuronPool type, which will have an array field for each of the different types of neurons. In the example above, the pool would have 3 fields each with one entry. Constrast this to the example below, where all of the neurons have the same channel types, and therefore there would be a neuron pool with 1 field with three entries:

.. code-block:: julia

	simulation_time=100.0 #ms
	myneuron1_1=deepcopy(myneuron1)
	myneuron1_2=deepcopy(myneuron1)
	mynetwork2=Network((myneuron1,myneuron1_1,myneuron1_2),simulation_time)
		

========================
Extracellular Recording
========================

JNeuron supports detecting the extracellular potential at 3D locations in a network of neurons. Electrodes are initialized by their 3D position and then can be added to the network:

.. code-block:: julia

	myelectrode=Extracellular([400.0,200.0,0.0])
	add!(mynetwork,myelectrode)

All of the necessary relationships between the collection of neurons in the network and the 3D position are calculated when the extracellular potential is added to the network. Multiple methods of calculating the extracellular potential, as well as different electrode shape approximates are supported (see Recording).

=========================
Extracellular Stimulation
=========================

Coming soon

========================
Intracellular Recording
========================

The intracellular potentials over the course of the simulation at particular locations in the neuron can be saved. To place an intracellular recording, you must specify 1) position of the neuron in the neuron pool, 2) the index of the neuron in your network, as well as 3) the particular node of interest.

.. code-block:: julia

	myintra=Intracellular(1,1,100)
	add!(mynetwork,myintra)

=========================
Intracellular Stimulation
=========================

A period of intracellular stimulation is defined by its 1) magnitude 2) time window of activity 3) and location. Once these are specified, a stimulation instance can be added to a network

.. code-block:: julia

	amp=2.0 #stimulation amplitude in nA
	neuron_type=1 #position of neuron type in neuron pool
	neuron_num=1 #index of neuron receiving stimulation
	node=40 #node of above neuron to input current
	tstart=1.0 #in ms
	tstop=2.0 #in ms
	mystim=Stim(amp,neuron_type,neuron_num,node,tstart,tstop)
	add!(mynetwork,mystim)





