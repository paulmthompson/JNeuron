

function gen_current{T<:Channel}(props::Array{T,1},myprop::Prop)

    confields=Array(Expr,length(props))
    curfields=Array(Expr,length(props))
    initfields=Array(Expr,length(props))
    mylengths=0

    for i=1:length(props)
        confields[i]=speed_con(props[i],mylengths)
        curfields[i]=speed_cur(props[i],mylengths)
        initfields[i]=speed_init(props[i],mylengths)
        mylengths=length(fieldnames(props[i]))-1
    end
    
    @eval begin

        function con!(p::$(typeof(myprop)),v::Array{Float64,1},internal::Array{Int64,1},dt=.025)

            @fastmath @simd for i=1:length(internal)
                k=internal[i]
                $([confields[j] for j=1:length(confields)]...)
            end
            
            nothing
        end

        function cur!(p::$(typeof(myprop)),v::Array{Float64,1},im::Array{Float64,1},internal::Array{Int64,1})

            @fastmath @simd for i=1:length(internal)
                k=internal[i]
                $([curfields[j] for j=1:length(curfields)]...)
            end

            nothing
        end

        function init!(p::$(typeof(myprop)),v::Array{Float64,1},internal::Array{Int64,1})
            @fastmath @simd for i=1:length(internal)
                k=internal[i]
                $([initfields[j] for j=1:length(initfields)]...)
            end
        end
        

            
    end    
end

function speed_cur(myprop::HH,l::Int64)

    :(begin
      
      @inbounds ina = p.p[$(5+l),i] * (v[k] - ena)
      @inbounds ik  = p.p[$(6+l),i] * (v[k] - ek)
      @inbounds il  = p.p[$(3+l),i] * (v[k] - p.p[$(4+l),i])

      @inbounds im[k]+=ina+ik+il

      end)
    
end

function speed_con(myprop::HH,l::Int64)

    :(begin

        #"m" sodium activation system
        @inbounds alpha = .1 * vtrap(-(v[k]+40.0),10.0)
        @inbounds beta =  4.0 * exp(-(v[k]+65.0)/18.0)
        mysum = alpha + beta
        @inbounds p.p[$(13+l),i] = 1.0/(q10*mysum)
        @inbounds p.p[$(10+l),i] = alpha/mysum
 
        #"h" sodium inactivation system
        @inbounds alpha = .07 * exp(-(v[k]+65.0)/20.0)
        @inbounds beta = 1.0 / (exp(-(v[k]+35.0)/10.0)+1)
        mysum = alpha + beta
        @inbounds p.p[$(14+l),i] = 1.0/(q10*mysum)
        @inbounds p.p[$(11+l),i] = alpha/mysum
    
        #"n" potassium activation system
        @inbounds alpha = .01*vtrap(-(v[k]+55.0),10.0) 
        @inbounds beta = .125*exp(-(v[k]+65.0)/80.0)
        mysum = alpha + beta
        @inbounds p.p[$(15+l),i] = 1/(q10*mysum)
        @inbounds p.p[$(12+l),i] = alpha/mysum

        @inbounds p.p[$(7+l),i]=implicit_euler(p.p[$(7+l),i],dt,p.p[$(13+l),i],p.p[$(10+l),i])
        @inbounds p.p[$(8+l),i]=implicit_euler(p.p[$(8+l),i],dt,p.p[$(14+l),i],p.p[$(11+l),i])
        @inbounds p.p[$(9+l),i]=implicit_euler(p.p[$(9+l),i],dt,p.p[$(15+l),i],p.p[$(12+l),i])

        @inbounds p.p[$(5+l),i] = p.p[$(1+l),i] * p.p[$(7+l),i] * p.p[$(7+l),i] * p.p[$(7+l),i] * p.p[$(8+l),i]
        @inbounds p.p[$(6+l),i] = p.p[$(2+l),i] * p.p[$(9+l),i] * p.p[$(9+l),i] * p.p[$(9+l),i] * p.p[$(9+l),i]


    end)
           
end

function speed_init(myprop::HH,l::Int64)
    
    :(begin

      @inbounds p.p[$(1+l),i]=.12
      @inbounds p.p[$(2+l),i]=.036
      @inbounds p.p[$(3+l),i]=.0003
      @inbounds p.p[$(4+l),i]=-54.3

        #"m" sodium activation system
        @inbounds alpha = .1 * vtrap(-(v[k]+40.0),10.0)
        @inbounds beta =  4.0 * exp(-(v[k]+65.0)/18.0)
        mysum = alpha + beta
        @inbounds p.p[$(13+l),i] = 1.0/(q10*mysum)
        @inbounds p.p[$(10+l),i] = alpha/mysum
 
        #"h" sodium inactivation system
        @inbounds alpha = .07 * exp(-(v[k]+65.0)/20.0)
        @inbounds beta = 1.0 / (exp(-(v[k]+35.0)/10.0)+1)
        mysum = alpha + beta
        @inbounds p.p[$(14+l),i] = 1.0/(q10*mysum)
        @inbounds p.p[$(11+l),i] = alpha/mysum
    
        #"n" potassium activation system
        @inbounds alpha = .01*vtrap(-(v[k]+55.0),10.0) 
        @inbounds beta = .125*exp(-(v[k]+65.0)/80.0)
        mysum = alpha + beta
        @inbounds p.p[$(15+l),i] = 1/(q10*mysum)
        @inbounds p.p[$(12+l),i] = alpha/mysum

        @inbounds p.p[$(7+l),i]=p.p[$(10+l),i]
        @inbounds p.p[$(8+l),i]=p.p[$(11+l),i]
        @inbounds p.p[$(9+l),i]=p.p[$(12+l),i]

        @inbounds p.p[$(5+l),i] = p.p[$(1+l),i] * p.p[$(7+l),i] * p.p[$(7+l),i] * p.p[$(7+l),i] * p.p[$(8+l),i]
        @inbounds p.p[$(6+l),i] = p.p[$(2+l),i] * p.p[$(9+l),i] * p.p[$(9+l),i] * p.p[$(9+l),i] * p.p[$(9+l),i]
    end)

end


const ena=50.0
const ek=-77.0
