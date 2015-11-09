
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

	type Network{T <: AbstractArray{Neuron,1}}
    		neur::T #Array of Neurons for simulation
    		t::FloatRange{Float64} #time range for simulation
    		extra::Array{Extracellular,1} #Extracellular Recording
    		intra::Array{Intracellular,1} #Intracellular Recording
    		stim::Array{Stim,1} #Intracellular Stimulation
	end

As you can see, the Network type contains spaces for any Neurons, extracellular electrodes or intracellular electrodes that you made add, as well as the current time current to all of these components. 

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
	
	run(mynetwork)

Any components that you specified to calculate (like extracellular voltage above) will be stored in the data container at the end of the run for you to access:

.. code-block:: julia
	
	mynetwork.extra[1].v #voltage calculated on extracellular electrode


******************
Network Components
******************

The components available to place into a network for simulation are 1) neurons, 2) intracellular stimulation electrodes, 3) intracellular recording, 4) extracellular recording, and 5) extracellular stimulation.

=======
Neuron
=======

The Neuron type is the data container with all of the variables necessary to solve the cable equation at each iteration. It can be created from 3D morphology data, and contain a variable of ion channels.

---------------------
Loading 3D Structure
---------------------

Neurons can be created from 3D reconstructions from imaging. Right now only Neurolucida file formats (.asc) are supported, but with others soon to come! The 3D data is first parsed by JNeuron to fill a  Import3D type that contains information about the 3D structure. Then this type can be used to generate the Neuron type that holds all of the data necessary to solve the cable equation.

.. code-block:: julia

	filepath="/path/to/file/cell.asc"
	myimport=input(filepath)
	myneuron=instantiate(myimport)

If you don't need the 3D information in Import3D, you can just call instantiate with a filepath ending in .asc

.. code-block:: julia

	filepath="/path/to/file/cell.asc"
	myneuron=instantiate(filepath)

----------------
Discretization
----------------

In JNeuron, the cable equation is discretized so that voltages are solved at certain points, or nodes (see Theory section). We have to decide for a given neuron shape, what is the appropriate number of nodes? By default, JNeuron uses the lambda-rule (see Theory). 

.. code-block:: julia

	set_nsegs!(myneuron)

----------------
Adding Channels
----------------

Neurons can have a variety of ion channels. Many from the literature are already defined in JNeuron (see Channels). If channels are present everywhere in the neuron, they can easily be added as follows:

.. code-block:: julia

	add!(myneuron,HH()) #add hodgkin huxley channels

For channels distributed non-uniforming:

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

The intracellular potentials over the course of the simulation at particular locations in the neuron can be saved. To place an intracellular recording, you must specify the index of the neuron in your network, as well as the particular node of interest.

.. code-block:: julia

	myintra=Intracellular(1,100)
	add!(mynetwork,myintra)

=========================
Intracellular Stimulation
=========================

A period of intracellular stimulation is defined by its 1) magnitude 2) time window of activity 3) and location. Once these are specified, a stimulation instance can be added to a network

.. code-block:: julia

	amp=2.0 #stimulation amplitude in nA
	neuron_num=1 #index of neuron receiving stimulation
	node=40 #node of above neuron to input current
	tstart=1.0 #in ms
	tstop=2.0 #in ms
	mystim=Stim(amp,neuron_num,node,tstart,tstop)
	add!(mynetwork,mystim)





