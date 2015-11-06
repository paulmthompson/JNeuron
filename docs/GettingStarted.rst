
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

In JNeuron, the cable equation is discretized so that voltages are solved at certain points, or nodes (see Theory section). 

----------------
Adding Channels
----------------

========================
Extracellular Recording
========================

=========================
Extracellular Stimulation
=========================

========================
Intracellular Recording
========================

=========================
Intracellular Stimulation
=========================

