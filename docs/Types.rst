
######
Types
######

***************
Neuron Related
***************

===========
Neuron Pool
===========

Imagine you have a similar with three kinds of neurons with different ion channels present in different parts of the cell:

.. code-block:: julia

	filepath="/path/to/file/cell.asc"
	blank_neuron=instantiate(filepath)
	set_nsegs!(blank_neuron)
	myneuron1=add(blank_neuron,[HH()])
	myneuron2=add(blank_neuron,[HH(),Passive()])
	myneuron3=add(blank_neuron,Array[[HH(),Passive()],[HH(),Passive()],[Passive()],[Passive()]]);

You could easily construct an array of the neurons above to be accessed during simulations:

.. code-block:: julia
	
	my_neur_array=Array(Neuron,3)
	my_neur_array[1]=myneuron1
	my_neur_array[2]=myneuron2
	my_neur_array[3]=myneuron3

But this is actually inefficient, because the methods used for solving the equation in the first one are different than the second, which are different than the third becuase each has different channel types: consequently, the compiler has to figure out which ones to use.

Consequently, you see above that the Network type contains a "Neuron Pool" which for the case above would look like this:

.. code-block:: julia

	type NeuronPool0
		N_1::Array{Neuron_1,1}
		N_2::Array{Neuron_2,1}
		N_3::Array{Neuron_3,1}
	end
