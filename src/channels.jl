
#=
Passive
still needs range variables
=#

function prop_init(prop::Prop,node::Node)
end

type Passive <: Prop
    nodevar::Array{ASCIIString,1}
    v::Float64
    i::Float64
    g::Float64
    e::Float64
end

function Passive()
    myvars=["v"]
    Passive(myvars,0.0,.001,.07)
end

function prop_calc(prop::Passive,node::Node)
    prop.i=prop.g*(node.vars["v"]-prop.e)
end
    
#=
HH
=#

type HH <: Prop
    nodevar::Array{ASCIIString,1}
    gnabar::Float64
    gkbar::Float64
    gl::Float64
    el::Float64
    gna::Float64
    gk::Float64
    il::Float64
    m::Float64
    h::Float64
    n::Float64
    minf::Float64
    hinf::Float64
    ninf::Float64
    ina::Float64
    ik::Float64
    mtau::Float64
    ntau::Float64
    htau::Float64
end

function HH()
    myvars=["v", "ena", "ek"]
    HH(.12,.036,.0003,-54.3,zeros(Float64,14)...)
end

function prop_calc(prop::HH,node::Node)

    rates(prop,node)

    prop.m = prop.m + (1. - exp(node.dt*(( ( ( - 1.0 ) ) ) / prop.mtau)))*(- ( ( ( prop.minf ) ) / prop.mtau ) / ( ( ( ( - 1.0) ) ) / prop.mtau ) - prop.m)
    prop.n = prop.n + (1. - exp(node.dt*(( ( ( - 1.0 ) ) ) / prop.ntau)))*(- ( ( ( prop.ninf ) ) / prop.ntau ) / ( ( ( ( - 1.0) ) ) / prop.ntau ) - prop.n)
    prop.h = prop.h + (1. - exp(node.dt*(( ( ( - 1.0 ) ) ) / prop.htau)))*(- ( ( ( prop.hinf ) ) / prop.htau ) / ( ( ( ( - 1.0) ) ) / prop.htau ) - prop.h)
        
    prop.gna = prop.gnabar * prop.m^3 * prop.h
    prop.ina = prop.gna * (node.vars["v"] - node.vars["ena"])
    prop.gk = prop.gkbar * prop.n^4
    prop.ik = prop.gk * (node.vars["v"] - node.vars["ek"])
    prop.il = prop.gl * (node.vars["v"] - prop.el)

    return prop.ina + prop.ik + prop.il
    
end

function prop_init(prop::HH,node::Node)
    rates(prop,node)
    prop.m=prop.minf
    prop.h=prop.hinf
    prop.n=prop.ninf
end

function rates(prop::HH,node::Node)

    q10 = 3^((node.temp - 6.3)/10)
    
    #"m" sodium activation system
    alpha = .1 * vtrap(-(v+40),10)
    beta =  4 * exp(-(v+65)/18)
    mysum = alpha + beta
    prop.mtau = 1/(q10*mysum)
    prop.minf = alpha/mysum
    
    #"h" sodium inactivation system
    alpha = .07 * exp(-(v+65)/20)
    beta = 1 / (exp(-(v+35)/10) + 1)
    mysum = alpha + beta
    prop.htau = 1/(q10*mysum)
    prop.hinf = alpha/mysum
    
    #"n" potassium activation system
    alpha = .01*vtrap(-(v+55),10) 
    beta = .125*exp(-(v+65)/80)
    mysum = alpha + beta
    prop.ntau = 1/(q10*mysum)
    prop.ninf = alpha/mysum

end

function vtrap(x::Float64,y::Float64)
    if abs(x/y) < 1e-6
        trap = y*(1 - x/y/2)
    else
        trap = x/(exp(x/y) - 1)
    end
    return trap
end





