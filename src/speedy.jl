

function gen_current{T<:Channel}(props::Array{T,1})

    confields=Array(Expr,length(props))
    curfields=Array(Expr,length(props))
    mylengths=0

    for i=1:length(props)
        confields[i]=speed_con(props[i],mylengths)
        curfields[i]=speed_cur(props[i],mylengths)
        mylengths=length(fieldnames(props[i]))-1
    end
    
    @eval begin

        function con(a::Array{Float64,2},v::Array{Float64,1},dt=.025)

            @fastmath @inbounds @simd for i=1:size(a,2)

                $(confields...)
                
            end
            
            nothing
        end

    end
    
end

function speed_con(myprop::HH,l::Int64)

    :(begin
    
        #"m" sodium activation system
        alpha = .1 * vtrap(-(v[i]+40.0),10.0)
        beta =  4.0 * exp(-(v[i]+65.0)/18.0)
        mysum = alpha + beta
        a[$(16+l),i] = 1.0/(q10*mysum)
        a[$(11+l),i] = alpha/mysum
    
        #"h" sodium inactivation system
        alpha = .07 * exp(-(v[i]+65.0)/20.0)
        beta = 1.0 / (exp(-(v[i]+35.0)/10.0) + 1.0)
        mysum = alpha + beta
        a[$(18+l),i] = 1/(q10*mysum)
        a[$(13+l),i] = alpha/mysum
    
        #"n" potassium activation system
        alpha = .01*vtrap(-(v[i]+55.0),10.0) 
        beta = .125*exp(-(v[i]+65.0)/80.0)
        mysum = alpha + beta
        a[$(17+l),i] = 1/(q10*mysum)
        a[$(12+l),i] = alpha/mysum
        
        a[$(8+l),i]=implicit_euler(a[$(8+l),i],dt,a[$(16+l),i],a[$(11+l),i])
        a[$(9+l),i]=implicit_euler(a[$(9+l),i],dt,a[$(17+l),i],a[$(12+l),i])
        a[$(10+l),i]=implicit_euler(a[$(10+l),i],dt,a[$(18+l),i],a[$(13+l),i])

        a[$(5+l),i] = a[$(1+l),i] * a[$(8+l),i] * a[$(8+l),i] * a[$(8+l),i] * a[$(10+l),i]
        a[$(6+l),i] = a[$(2+l),i] * a[$(9+l),i] * a[$(9+l),i] * a[$(9+l),i] * a[$(9+l),i]

    end)
    
end

function speed_cur(myprop::HH,l::Int64)

    a[$(14+l),i] = a[$(5+l),i] * (v[i] - vars["ena"][i])
    a[$(15+l),i] = a[$(6+l),i] * (v[i] - vars["ek"][i])
    a[$(7+l),i] = a[$(3+l),i] * (v[i] - a[$(4+l),i])

    nothing
    
end
