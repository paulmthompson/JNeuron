

function gen_current{T<:Channel}(props::Array{T,1})

    confields=Array(Array{Expr,1},length(props))
    curfields=Array(Expr,length(props))
    mylengths=0

    for i=1:length(props)
        confields[i]=speed_con(props[i],mylengths)
        curfields[i]=speed_cur(props[i],mylengths)
        mylengths=length(fieldnames(props[i]))-1
    end
    
    @eval begin

        function con(a::Array{Float64,2},v::Array{Float64,1},temp::Array{Array{Float64,1},1},dt=.025)

            @fastmath @inbounds @simd for i=1:size(a,2)
                $([confields[i][1] for i=1:length(confields)]...)
            end

            $([confields[i][2] for i=1:length(confields)]...)

            @fastmath @simd for i=1:size(a,2)
                $([confields[i][3] for i=1:length(confields)]...)
            end

            $([confields[i][4] for i=1:length(confields)]...)

            @fastmath @simd for i=1:size(a,2)
                $([confields[i][5] for i=1:length(confields)]...)
            end
            
            nothing
        end

        function cur(a::Array{Float64,2},vars::Dict{ASCIIString,Float64},v::Array{Float64,1})

            
            
        end

    end
    
end

function speed_cur(myprop::HH,l::Int64)

    :(begin
    a[$(14+l),i] = a[$(5+l),i] * (v[i] - vars["ena"][i])
    a[$(15+l),i] = a[$(6+l),i] * (v[i] - vars["ek"][i])
    a[$(7+l),i] = a[$(3+l),i] * (v[i] - a[$(4+l),i])

      end)
    
end

function speed_con(myprop::HH,l::Int64)

    ex=Array(Expr,5)
    
    #calculate exp with yeppp

    ex[1] = :(begin

              temp[1][i]=-(v[i]+65.0)/18.0
              temp[2][i]=-(v[i]+65.0)/20.0
              temp[3][i]=-(v[i]+35.0)/10.0
              temp[4][i]=-(v[i]+65.0)/80.0
              temp[5][i]=-(v[i]+40.0)/10.0
              temp[6][i]=-(v[i]+55.0)/10.0

              end)

    ex[2] = :(begin

              Yeppp.exp!(temp[1])
              Yeppp.exp!(temp[2])
              Yeppp.exp!(temp[3])
              Yeppp.exp!(temp[4])
              Yeppp.exp!(temp[5])
              Yeppp.exp!(temp[6])

            end)

    ex[3] = :(begin

              
        #"m" sodium activation system
        @inbounds alpha = .1 * vtrap_s(-(v[i]+40.0),10.0,temp[5][i])
        @inbounds beta =  4.0 * temp[1][i]
        mysum = alpha + beta

        @inbounds a[$(16+l),i] = 1.0/(q10*mysum)
        @inbounds a[$(11+l),i] = alpha/mysum
        @inbounds temp[5][i]=dt*-1.0/a[$(16+l),i]
 
        #"h" sodium inactivation system
        @inbounds alpha = .07 * temp[2][i]
        @inbounds beta = 1.0 / (temp[3][i] + 1.0)
        mysum = alpha + beta
        @inbounds a[$(18+l),i] = 1.0/(q10*mysum)
        @inbounds a[$(13+l),i] = alpha/mysum
        @inbounds temp[6][i]=dt*-1.0/a[$(18+l),i]
    
        #"n" potassium activation system
        @inbounds alpha = .01*vtrap_s(-(v[i]+55.0),10.0,temp[6][i]) 
        @inbounds beta = .125*temp[4][i]
        mysum = alpha + beta
        @inbounds a[$(17+l),i] = 1/(q10*mysum)
        @inbounds a[$(12+l),i] = alpha/mysum
        @inbounds temp[7][i]=dt*-1.0/a[$(17+l),i]
        
        end)

    ex[4] = :(begin
    
    #calculate exptau

              Yeppp.exp!(temp[7])
              Yeppp.exp!(temp[8])
              Yeppp.exp!(temp[9])

              end)

    ex[5] = :(begin
              
        @inbounds a[$(8+l),i]=implicit_euler_s(a[$(8+l),i],temp[7][i],a[$(16+l),i],a[$(11+l),i])
        @inbounds a[$(9+l),i]=implicit_euler_s(a[$(9+l),i],temp[9][i],a[$(17+l),i],a[$(12+l),i])
        @inbounds a[$(10+l),i]=implicit_euler_s(a[$(10+l),i],temp[8][i],a[$(18+l),i],a[$(13+l),i])

        @inbounds a[$(5+l),i] = a[$(1+l),i] * a[$(8+l),i] * a[$(8+l),i] * a[$(8+l),i] * a[$(10+l),i]
        @inbounds a[$(6+l),i] = a[$(2+l),i] * a[$(9+l),i] * a[$(9+l),i] * a[$(9+l),i] * a[$(9+l),i]

        end)

    ex
    
end

function implicit_euler_s(x::Float64,exptau::Float64,xtau::Float64,xinf::Float64)
    x + (1.0 - exptau)*(- ( ( ( xinf ) ) / xtau ) / ( ( ( ( - 1.0) ) ) / xtau ) - x)
end

function vtrap_s(x::Float64,y::Float64,expxy::Float64)

    if abs(x/y) < 1e-6
        trap = y*(1 - x/y/2)
    else
        trap = x/(expxy - 1)
    end
    trap

end

function speedy(a::Array{Float64,2},v::Array{Float64,1})

        #"m" sodium activation system
        @inbounds alpha = .1 * vtrap(-(v[i]+40.0),10.0)
        @inbounds beta =  4.0 * exp(-(v[i]+65.0)/18.0)
        mysum = alpha + beta

        @inbounds a[$(16+l),i] = 1.0/(q10*mysum)
        @inbounds a[$(11+l),i] = alpha/mysum
 
        #"h" sodium inactivation system
        @inbounds alpha = .07 * exp(-(v[i]+65.0)/20.0)
        @inbounds beta = 1.0 / (exp(-(v[i]+35.0)/10.0)+1)
        mysum = alpha + beta
        @inbounds a[$(18+l),i] = 1.0/(q10*mysum)
        @inbounds a[$(13+l),i] = alpha/mysum
    
        #"n" potassium activation system
        @inbounds alpha = .01*vtrap(-(v[i]+55.0),10.0) 
        @inbounds beta = .125*exp(-(v[i]+65.0)/80.0)
        mysum = alpha + beta
        @inbounds a[$(17+l),i] = 1/(q10*mysum)
        @inbounds a[$(12+l),i] = alpha/mysum

        @inbounds a[$(8+l),i]=implicit_euler_s(a[$(8+l),i],dt,a[$(16+l),i],a[$(11+l),i])
        @inbounds a[$(9+l),i]=implicit_euler_s(a[$(9+l),i],dt,a[$(17+l),i],a[$(12+l),i])
        @inbounds a[$(10+l),i]=implicit_euler_s(a[$(10+l),i],dt,a[$(18+l),i],a[$(13+l),i])

        @inbounds a[$(5+l),i] = a[$(1+l),i] * a[$(8+l),i] * a[$(8+l),i] * a[$(8+l),i] * a[$(10+l),i]
        @inbounds a[$(6+l),i] = a[$(2+l),i] * a[$(9+l),i] * a[$(9+l),i] * a[$(9+l),i] * a[$(9+l),i]

end
