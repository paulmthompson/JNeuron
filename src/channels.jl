
export Passive,HH


#=
Passive
still needs range variables
adapted for Julia by PMT
=#

function implicit_euler(x::Float64, dt::Float64,xtau::Float64,xinf::Float64)
    x = x + (1.0 - exp(dt*(( ( ( - 1.0 ) ) ) / xtau)))*(- ( ( ( xinf ) ) / xtau ) / ( ( ( ( - 1.0) ) ) / xtau ) - x)
end

function prop_init(prop::Channel,v::Float64)
end

type Passive <: Channel
    nodevar::Array{ASCIIString,1}
    g::Float64
    e::Float64
end

function Passive()
    myvars=Array(ASCIIString,0)
    Passive(myvars,.001,.07)
end

function con_calc(prop::Passive,v::Float64,dt::Float64)

end

function cur_calc(prop::Passive,vars::Dict{ASCIIString,Float64},v::Float64)
    prop.g*(v-prop.e)
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

const q10 = 3^((6.3 - 6.3)/10)

type HH <: Channel
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
    myvars=["ena", "ek"]
    HH(myvars,.12,.036,.0003,-54.3,zeros(Float64,14)...)
end

function con_calc(prop::HH,v::Float64,dt::Float64)
    rates(prop,v)

    prop.m=implicit_euler(prop.m,dt,prop.mtau,prop.minf)
    prop.n=implicit_euler(prop.n,dt,prop.ntau,prop.ninf)
    prop.h=implicit_euler(prop.h,dt,prop.htau,prop.hinf)
  
    prop.gna = prop.gnabar * prop.m^3 * prop.h
    prop.gk = prop.gkbar * prop.n^4

    nothing
    
end

function cur_calc(prop::HH,vars::Dict{ASCIIString,Float64},v::Float64)
    
    prop.ina = prop.gna * (v - vars["ena"])
    prop.ik = prop.gk * (v - vars["ek"])
    prop.il = prop.gl * (v - prop.el)

    prop.ina + prop.ik + prop.il
end


function prop_init(prop::HH,v::Float64)
    rates(prop,v)
    prop.m=prop.minf
    prop.h=prop.hinf
    prop.n=prop.ninf

    prop.gna = prop.gnabar * prop.m^3 * prop.h
    prop.gk = prop.gkbar * prop.n^4
    
    nothing
end

function rates(prop::HH,v::Float64)

    #q10 = 3^((node.temp - 6.3)/10)
    
    
    #"m" sodium activation system
    alpha = .1 * vtrap(-(v+40),10.0)
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
    alpha = .01*vtrap(-(v+55),10.0) 
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
    trap
end

#=
Ca_HVA

Reuveni, Friedman, Amitai, and Gutnick, J.Neurosci. 1993
=#

type Ca_HVA <: Channel
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
    myvars=["eca"]
    Ca_HVA(myvars,.00001,zeros(Float64,7)...)
end

function con_calc(prop::Ca_HVA,v::Float64,dt::Float64)

    rates(prop,v)

    prop.m=implicit_euler(prop.m,dt,prop.mtau,prop.minf)
    prop.h=implicit_euler(prop.h,dt,prop.htau,prop.hinf)

    prop.gCa_HVA = prop.gCa_HVAbar*prop.m^2*prop.h

    nothing
end

function cur_calc(prop::Ca_HVA,vars::Dict{ASCIIString,Float64},v::Float64)

    prop.ica = prop.gCA_HVA*(v-vars["eca"])
    
end

function prop_init(prop::Ca_HVA,v::Float64)
    rates(prop,v)
    prop.m=prop.mInf
    prop.h=prop.hInf
    nothing
end

function rates(prop::Ca_HVA,v::Float64)
    if v == -27.0
        v+=.0001
    else
    end

    mAlpha =  (0.055*(-27-v))/(exp((-27-v)/3.8) - 1)        
    mBeta  =  (0.94*exp((-75-v)/17))

    prop.mInf = mAlpha/(mAlpha + mBeta)
    prop.mTau = 1/(mAlpha + mBeta)

    hAlpha =  (0.000457*exp((-13-v)/50))
    hBeta  =  (0.0065/(exp((-v-15)/28)+1))

    prop.hInf = hAlpha/(hAlpha + hBeta)
    prop.hTau = 1/(hAlpha + hBeta)
    nothing
end

#=
LVA ca channel

Note: mtau is an approximation from the plots
Avery and Johnston 1996, tau from Randall 1997
shifted by -10 mv to correct for junction potential
corrected rates using q10 = 2.3, target temperature 34, orginal 21
=#

type Ca_LVAst <: Channel
    nodevar::Array{ASCIIString,1}
    gCa_LVAstbar::Float64
    gCA_LVAst::Float64
    m::Float64
    h::Float64
    minf::Float64
    hinf::Float64
    ica::Float64
    mtau::Float64
    htau::Float64
    qt::Float64
end

function Ca_LVAst()
    myvars=["eca"]
    Ca_LVAst(myvars,.00001,zeros(Float64,7)...,2.3^((34-21)/10))
end

function con_calc(prop::Ca_LVAst,v::Float64,dt::Float64)

    rates(prop,v)

    prop.m=implicit_euler(prop.m,dt,prop.mtau,prop.minf)
    prop.h=implicit_euler(prop.h,dt,prop.htau,prop.hinf)
    
    prop.gCa_LVAst = prop.gCa_LVAstbar*prop.m^2*prop.h

    nothing
end

function cur_calc(prop::Ca_LVAst,vars::Dict{ASCIIString,Float64},v::Float64)

    prop.ica = prop.gCA_LVAst*(v-vars["eca"])
    
end

function prop_init(prop::Ca_LVAst,v::Float64)
    rates(prop,v)
    prop.m=prop.mInf
    prop.h=prop.hInf
    nothing
end

function rates(prop::Ca_LVAst,v::Float64)

    prop.mInf = 1.0000/(1+ exp((v - -30.000)/-6))
    prop.mTau = (5.0000 + 20.0000/(1+exp((v - -25.000)/5)))/prop.qt
    prop.hInf = 1.0000/(1+ exp((v - -80.000)/6.4))
    prop.hTau = (20.0000 + 50.0000/(1+exp((v - -40.000)/7)))/prop.qt
    
    nothing
end

type Ih <: Channel
    nodevar::Array{ASCIIString,1}
    gIhbar::Float64
    ehcn::Float64
    gIh::Float64
    m::Float64
    minf::Float64
    ihcn::Float64
    mtau::Float64
end

function Ih()
    myvars=Array(ASCIIString,0)
    Ih(myvars,.00001, -45, zeros(Float64, 5)...)
end

function con_calc(prop::Ih,v::Float64,dt::Float64)

    rates(prop,v)

    prop.m=implicit_euler(prop.m,dt,prop.mtau,prop.minf)
    
    prop.gIh = prop.gIhbar*prop.m

    nothing
end

function cur_calc(prop::Ih,vars::Dict{ASCIIString,Float64},v::Float64)

    prop.ica = prop.gIh*(v-prop.ehcn)
    
end

function prop_init(prop::Ih,v::Float64)
    rates(prop,v)
    prop.m=prop.mInf
    nothing
end

function rates(prop::Ih,v::Float64)

    if v == -154.9
        v=+ .0001
    else
    end
    
    mAlpha =  0.001*6.43*(v+154.9)/(exp((v+154.9)/11.9)-1)
    mBeta  =  0.001*193*exp(v/33.1)
    prop.mInf = mAlpha/(mAlpha + mBeta)
    prop.mTau = 1/(mAlpha + mBeta)
    
    nothing
end

