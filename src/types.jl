
#=
Data structures need to be organized to either emphasize topology, or computational efficiency. The computational data structures are used to set most things up, so we can get contiguous blocks of floats for matrix math. 

Each neuron is divided into sections (think branches or compartments) and each section will have multiple segments, with each segment having one node. The node is where the action is located: it is where we calculate the voltage for a segment. Every time step, we calculate a change in voltage, and add it to the old voltages at the nodes to get new voltages. Therefore we need arrays of floats to represent:

1) delta_v=Array(Float64,nodenum)
2) v_old=Array(Float64,nodenum)
3) v_new=Array(Float64,nodenum)

Unfortunately because of a lot of complexities in discretizing the cable equation, maintaining stability and low error, we don't end up with a nice equation with a linear relationship between the old voltage and new voltage:

A * v_old = V_new

We can however create a matrix problem of the following form by discretizing the cable equation:

A(t) * delta_v = current_into_nodes(t)

Note that since we are iterating, the current_into_nodes term will be time dependent (or need to be updated every iteration) and the matrix relating the delta_v terms to one another also need to be updated, but luckily only along the diagonal.

So we need to construct the data structure for the relationships between the change in voltages between nodes:

4) A=Array(Float64,nodenum,nodenum)

I mentioned that the diagonal elements of A are the only terms that change, so we might want to keep those handy

5) diag=diagview(A)

Additionally, we should keep the original diagonal values to allow for resetting after each iteration

6) diag_old=diag(A)

We will also want a vector of the currents into each node at the previous time step. because current_into_nodes is cumbersome, we will adopt the neuron notation of rhs

6) rhs=Array(Float64,nodenum)

We will then solve for the change in voltages with the following:

A \ rhs = delta_v

In addition to the voltages at nodes, there are other meaningful quantities that we want to keep track of at each node, like what channel equations we will use, what the resistances are between other nodes, what nodes it is connected to etc. We should therefore keep a vector of nodes:

7) nodes = Array(Node,nodenum)

Where each node is its own type, which keeps track of that information

These data containers don't help us to easily see the big picture of how the neuron is structured. This falls under the concept of a Section, which can be thought of as a continous piece of a dendrite, soma or axon. Dendritic branches would be connected together to form the dentritic tree. Each Section would keep track of what it is connected to, and also what nodes are included in it. So while for the actual simulation a "Section" is not an incredibly important concept, it is very important for the abstract from a neuron to the math to solve for voltage at different locations. The main data structure should contain a list of these sections:

8) secstack = Array(Section, sectionnum)


=#


abstract Prop #Property of section (HH, passive etc). contains all of the stuff you need to calc things

#associated 3d point
type Pt3d
    x::Float64
    y::Float64
    z::Float64
    d::Float64
    arc::Float64 #normalized distance from 0 to end
end

type Node
    ind::Int64
    vars::Dict{ASCIIString,Float64}
    area::Array{Float64,1} #surface area of left [1] and right[2] part of segment
    ri::Array{Float64,1}  #internal resistance of left[1] and right[2] part of segment
    b::Float64 #resistance between node and parent node
    a::Float64 #resistance between node and each child node
    parent::Int64 #index in node array of parent
    children::Array{Int64,1} #Node(s) from other sections attached to this one
    internal::Bool
    pt3d::Array{Pt3d,1}
    prop::Array{Prop,1} #Array of abstract types (yucky!), each a subtype of 
end

type Node_ext
    ind::Int64
    ri::Array{Float64,1}
    parent_r::Float64
    children_r::Array{Float64,1}
end

type Section
    refcount::Int64 #ID for section, also the place in secstack array
    mtype::Int64 #Cellbody=1,Axon=2,Dendrite=3,Apical=4
    pnode::Array{Node,1} #one node at center of each segment
    child::Array{Section,1}
    pt3d::Array{Pt3d,1}
    Ra::Float64 #cytoplasmic resistivity
    parentx::Float64
    length::Float64
end

type Neuron
    secstack::Array{Section,1}
    A::SparseMatrixCSC{Float64,Int64}
    v::Array{Float64,1} #intracellular voltage
    delta_v::Array{Float64,1} #change in membrane voltage
    rhs::Array{Float64,1}
    Ra::Float64
    Cm::Float64
    dt::Float64
    nodes::Array{Node,1}
    i_vm::Array{Float64,1}
    divm::Array{Float64,1}
    diag_old::Array{Float64,1}

    #extracellular data structures (where to put these? different type probably)
    #enodes::Array{Node_ext,1}
    #vext::Array{Float64,1}
    #delta_vext::Array{Float64,1}
    #A_ext::Array{Float64,2}
    #rhs_ext::Array{Float64,1}
    #diag_ext::Array{Float64,1}
    #diag_ext_old::Array{Float64,1}
end

function Neuron()
    Neuron(Array(Section,0),spzeros(Float64,0,0),zeros(Float64,0),zeros(Float64,0),zeros(Float64,0),0.0,0.0,0.0,Array(Node,0),zeros(Float64,0),zeros(Float64,0),zeros(Float64,0))
end

type Network
    neur::Array{Neuron,1}
    t::FloatRange{Float64}
end

function Network(neuron::Neuron,tstop::Float64)
    Network([neuron],0.0:neuron.dt:tstop)
end

function run(network::Network)

end


myconstants=Dict{ASCIIString, Float64}("ena"=>50.0, "ek"=>-77.0)
