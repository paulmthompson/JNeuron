
##########
Theory
##########

**************
Cable Equation
**************

The spatial and temporal distribution of voltage along a neuron can be represented by the following partial differential equation:

.. math:: C \frac{\partial{V}}{\partial{t}} = \frac{a}{2R_a}\frac{\partial^2{V}}{\partial{x^2}} - i_m(V)

Where :math:`C` is the specific membrane capacitance :math:`(F/cm^2)`, :math:`V` is the membrane voltage, :math:`a` is the radius of the neuron (cm), :math:`R_a` is the cytoplasmic resistivity :math:`(\Omega cm)`, and :math:`i_m(V)` is the membrane current from either active (ion channels) or passive processes.

To simulate a neuron, we will need to discretize the above equation in both the spatial and temporal domain, while maintaining a desired level of accuracy. Numerous publications have described how to attain a second order correct solution; Mascagni 1989 is particularly informative.

=======================
Spatial Discretization
=======================

----------
Conceptual
----------

A real neuron will have a membrane potential that is continuously distributed along the length of its soma, axon and dendrites. We can't calculate this for arbitrary shapes, so we'll have to settle for determing voltages at discrete locations. What should these locations be?

If we wanted to divide a neuron into pieces or <b>segments</b> as we shall call them, a natural location to cut a neuron into pieces would be at its own divisions in structure. NEURON uses the term, and the term we will also adopt here, of <b>section</b> to describe a continuous piece of an unbranching cable. We will therebefore want to calculate voltages at the intersections of all of these sections. How about along a section? We can divide a section into smaller pieces (segments) as well, and calculate voltage at these discrete points. 

-----
Math
-----

Rather than considering the continuous cable equation above, we can consider a neuron that is divided into segments, each of length :math:`\Delta x`. We can then replace the second order patial derivative with a second order correct finite difference:

.. math:: C \frac{\partial{V}}{\partial{t}} = \frac{a}{2R_a} \frac{V_{i+1}-2V_i+V_{i-1}}{(\Delta x)^2} - i_m(V_i), \ \ i=0, ..., N

Now we will have a system of equations, one for each of the N segments.

=======================
Temporal Discretization
=======================

-----------
Conceptual
-----------

Most of us have probably at some time or another solved an ordinary differential equation, and used one of the many easy to use computational tools, such as matlabs ode45 or ode23s to do so. Unforunately, the cable equation is not an ODE, but a parabolic partial differential equation, and the solution is going to be more tedious and idiosyncratic to the problem.

The main points are that we need to maintain some desired level of accuracy, be it first order or second order; and we need to consider the stability of our solution. You can use a so called forward difference method on the cable equation above, and this would be similar to some of the under the hood methods in solvers like ode45, but it would NOT be stable. Instead we need to use so called implicit methods, where the solution for time t=t1 depends on some values at t=t2. As you would expect, we will need to employ some trickery to approximate the values we need at t=t2, since we haven't solved for that time yet.

------
Math
------

The simpliest implicit method is called the Implicit Euler method, and will approximate the above equation at finite time steps as

.. math:: C \frac{V_i^{n+1}-V_i^n}{\Delta t} = \frac{a}{2R_a} \frac{V_{i+1}^{n+1}-2V_i^{n+1}+V_{i-1}^{n+1}}{(\Delta x)^2} - i_m(V_i^{n+1}), \ \ i=0, ..., N

As you can see above, this solution depends on the current at a future time, which we will have to approximate. We will need to do this and some reordering to make this equation something easier to solve. First, for reasons that will become clear later, we should arrange our voltages in terms of voltage at the current time, and the change in voltage, or :math:`\Delta V_i`. Applying this change of variables gives:

.. math:: C \frac{\Delta V_i}{\Delta t} = \frac{a}{2R_a} \frac{(\Delta V_{i+1} + V_{i+1})-2(\Delta V_i + V_i)+(\Delta V_{i-1} + V_{i-1})}{(\Delta x)^2} - i_m(V_i^{n+1})

.. math:: C \frac{\Delta V_i}{\Delta t} - \frac{a}{2R_a} \frac{\Delta V_{i+1} - \Delta V_i - \Delta V_i + \Delta V_{i-1}}{(\Delta x)^2}= \frac{a}{2R_a} \frac{V_{i+1} - V_i - V_i + V_{i-1}}{(\Delta x)^2} - i_m(V_i^{n+1})

The term :math:`\frac{a}{2 R_a (\Delta x)^2}` can be expressed in more meaningful terms as :math:`\frac{1}{a_i R_{ij}}` where :math:`a_i` is surface area and :math:`R_{ij}` is the resistance between the sections. We will make this substitution and change our notation slightly to view the voltage at a node not dependent on the nodes on either side but on its "parent" node and "child" nodes. By definition, segment will have a parent (except the first segment), and can have 0 to M children.

.. math:: C \frac{\Delta V_i}{\Delta t} + \frac{\Delta V_p - \Delta V_i}{a_i R_{pi}} - \sum \frac{\Delta V_c - \Delta V_i}{a_i R_{ci}} = \frac{V_p - V_i}{a_i R_{pi}} - \sum \frac{V_c - V_i}{a_i R_{ci}} - i_m(V_i^{n+1}) - i_m(V_i^{n+1})

Where :math:`R_{pi}` is the resistance between a given node and its parent, :math:`R_{ci}` is the resistance between a given node and one of its children, and :math:`a_i` is the surface area of the given node.

The equation is almost solvable: we can calculate all of the parameters by the geometry of the neuron, but we do not know the membrane current at the next time step. We can approximate this current as the following:

.. math:: i_i(V_i(t + \Delta t)) = i(V_i(t)) + \Delta V_i \frac{d i_i}{d V_i}

Where we can approximate the instanenous conductance in our simulation by evaluating the current at the known voltage and then at v=v+.001. This gives us the final equation:

.. math:: C \frac{\Delta V_i}{\Delta t} +\frac{d i_i}{d V_i} \Delta V_i + \frac{\Delta V_p - \Delta V_i}{a_i R_{pi}} - \sum \frac{\Delta V_c - \Delta V_i}{a_i R_{ci}} = \frac{V_p - V_i}{a_i R_{pi}} - \sum \frac{V_c - V_i}{a_i R_{ci}} - i_m(V_i)

=============
Matrix Math
=============

-----------
Conceptual
-----------

Everything our the right hand side of the above equation: 1) the current density entering a given segment from the parent segment; 2) the current density(s) leaving the given segment to its child segments, and 3) the membrane current due to activate and passive components at the current time step can all be calculated during a given iteration. 

On the left hand side, every coefficient of the :math:`\Delta V` terms are constant with the exception of the di/dv term which changes from iteration to iteration. Therefore, the most efficient way to simultaneously solve for the change in voltage at all of the segments would be to calculate the constant coefficients for the left hand side, which will complete a matrix, A, relating the change in voltage at one node to the change in voltage at another.

.. math:: A \Delta V = \frac{V_p - V_i}{a_i R_{pi}} - \sum \frac{V_c - V_i}{a_i R_{ci}} - i_m(V_i)

At every iteration then, we calculate that right hand term for each segment, calculate di/dv for each segment and add it to the diagonal term, and then solve for :math:`\Delta V` by simple matrix math.

------
Math
------



***********
References
***********
