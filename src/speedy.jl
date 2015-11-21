

function gen_current()

    @eval begin

        function myspeed(a::Array{Float64,2},v::Array{Float64,1})

            @inbounds @fastmath @simd for i=1:size(a,2)

                $(speedy())
                
            end
            
            nothing
        end

    end
    
end

function speedy()

    :(begin
    
        #"m" sodium activation system
        alpha = .1 * vtrap(-(v[i]+40),10.0)
        beta =  4 * exp(-(v[i]+65)/18)
        mysum = alpha + beta
        a[16,i] = 1/(q10*mysum)
        a[11,i] = alpha/mysum
    
        #"h" sodium inactivation system
        alpha = .07 * exp(-(v[i]+65)/20)
        beta = 1 / (exp(-(v[i]+35)/10) + 1)
        mysum = alpha + beta
        a[18,i] = 1/(q10*mysum)
        a[13,i] = alpha/mysum
    
        #"n" potassium activation system
        alpha = .01*vtrap(-(v[i]+55),10.0) 
        beta = .125*exp(-(v[i]+65)/80)
        mysum = alpha + beta
        a[17,i] = 1/(q10*mysum)
        a[12,i] = alpha/mysum
        
        a[8,i]=implicit_euler(a[8,i],.025,a[16,i],a[11,i])
        a[9,i]=implicit_euler(a[9,i],.025,a[17,i],a[12,i])
        a[10,i]=implicit_euler(a[10,i],.025,a[18,i],a[13,i])

    end)
    
end
