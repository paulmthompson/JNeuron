
# Overview

# Derivation

The generic cable equation for a passive piece of membrane is:

where C is Capacitance, V is membrane voltage, R is the resistance and g is the passive membrane conductance. 

To solve this equation for voltage for space and time, we must spatially discretize the equation, and then discretize the time. Our spatial discretization procedure must be mindful of 1) changing resistances between nodes and 2) possible branch points where the cable can split and merge. The temporal discretization needs to be stable. Additionally, the above cable equation needs to be modified to account for active components, which will also have their own differential equations. These couple differential equations must be solved on time scales that maintain their accuracy levels over time. Finally, we want to do all of this as efficiently as possible.


