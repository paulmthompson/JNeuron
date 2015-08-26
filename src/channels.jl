
#=
Passive
still needs range variables
adapt for Julia by PMT
=#

function prop_init(prop::Prop,node::Node)
end

function implicit_euler(x::Float64, dt::Float64,xtau::Float64,xinf::Float64)
    x = x + (1.0 - exp(dt*(( ( ( - 1.0 ) ) ) / xtau)))*(- ( ( ( xinf ) ) / xtau ) / ( ( ( ( - 1.0) ) ) / xtau ) - x)
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
 This is the original Hodgkin-Huxley treatment for the set of sodium, 
  potassium, and leakage channels found in the squid giant axon membrane.
  ("A quantitative description of membrane current and its application 
  conduction and excitation in nerve" J.Physiol. (Lond.) 117:500-544 (1952).)
 Membrane voltage is in absolute mV and has been reversed in polarity
  from the original HH convention and shifted to reflect a resting potential
  of -65 mV.
SW Jaslove  6 March, 1992
adapted for Julia by PMT
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
    HH(myvars,.12,.036,.0003,-54.3,zeros(Float64,14)...)
end

function prop_calc(prop::HH,node::Node)

    rates(prop,node)

    prop.m=implicit_euler(prop.m,node.dt,prop.mtau,prop.minf)
    prop.n=implicit_euler(prop.n,node.dt,prop.ntau,prop.ninf)
    prop.h=implicit_euler(prop.h,node.dt,prop.htau,prop.hinf)
  
    prop.gna = prop.gnabar * prop.m^3 * prop.h
    prop.ina = prop.gna * (node.vars["v"] - node.vars["ena"])
    prop.gk = prop.gkbar * prop.n^4
    prop.ik = prop.gk * (node.vars["v"] - node.vars["ek"])
    prop.il = prop.gl * (node.vars["v"] - prop.el)

    prop.ina + prop.ik + prop.il
    
end

function prop_init(prop::HH,node::Node)
    rates(prop,node)
    prop.m=prop.minf
    prop.h=prop.hinf
    prop.n=prop.ninf
    nothing
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

    nothing
end

function vtrap(x::Float64,y::Float64)
    if abs(x/y) < 1e-6
        trap = y*(1 - x/y/2)
    else
        trap = x/(exp(x/y) - 1)
    end
    return trap
end

#=
Ca_HVA

Reuveni, Friedman, Amitai, and Gutnick, J.Neurosci. 1993
=#

type Ca_HVA <: Prop
    nodevar::Array{ASCIIString,1}
    gCa_HVAbar::Float64
    gCA_HVA::Float64
    m::Float64
    h::Float64
    minf::Float64
    hinf::Float64
    ica::Float64
    mtau::Float64
    htau::Float64
end

function Ca_HVA()
    myvars=["v", "eca"]
    Ca_HVA(myvars,.00001,zeros(Float64,7)...)
end

function prop_calc(prop::Ca_HVA, node::Node)
    rates(prop,node)

    prop.m=implicit_euler(prop.m,node.dt,prop.mtau,prop.minf)
    prop.h=implicit_euler(prop.h,node.dt,prop.htau,prop.hinf)
    
    prop.gCa_HVA = prop.gCa_HVAbar*prop.m^2*prop.h
    prop.ica = prop.gCA_HVA*(node.vars["v"]-node.vars["eca"])
end

function prop_init(prop::Ca_HVA, node::Node)
    rates(prop,node)
    prop.m=prop.mInf
    prop.h=prop.hInf
    nothing
end

function rates(prop::Ca_HVA, node::Node)
    if node.vars["v"] == -27.0
        v=v+.0001
    end

    mAlpha =  (0.055*(-27-node.vars["v"]))/(exp((-27-node.vars["v"])/3.8) - 1)        
    mBeta  =  (0.94*exp((-75-node.vars["v"])/17))

    prop.mInf = mAlpha/(mAlpha + mBeta)
    prop.mTau = 1/(mAlpha + mBeta)

    hAlpha =  (0.000457*exp((-13-node.vars["v"])/50))
    hBeta  =  (0.0065/(exp((-node.vars["v"]-15)/28)+1))

    prop.hInf = hAlpha/(hAlpha + hBeta)
    prop.hTau = 1/(hAlpha + hBeta)
    nothing
end

