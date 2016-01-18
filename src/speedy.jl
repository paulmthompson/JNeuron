

function speed_const(myprop::Prop,props...)

    ex=Array(Symbol,0)
    
    for i=1:length(props)
        mystart=1
        count=1
        for j=1:props[i].nodevar
            push!(ex,(symbol("p.p[i][$(count)]")))
            mystart+=1
            count+=1
        end 

        for j=mystart:length(fieldnames(props[i]))-1
            push!(ex,(symbol("p_$(count)")))
            count+=1
        end     
    end
    ex    
end

function gen_current(props,myprop::Prop)

    constfields=speed_const(myprop,props)
    
    if issubtype(typeof(props),Tuple)
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
    else
        confields=[speed_con(props,0)]
        curfields=[speed_cur(props,0)]
        initfields=[speed_init(props,0)]
    end
      
    @eval begin

        function con!(p::$(typeof(myprop)),v::Array{Float64,1},internal::Array{Int64,1},dt=.025)

            for i=1:length(internal)
                k=internal[i]
                $([confields[j] for j=1:length(confields)]...)
                p.p[i]=$(typeof(myprop.p))([constfields[j]  for j=1:length(constfields)]...)
            end
            
            nothing
        end

        function cur!(p::$(typeof(myprop)),v::Array{Float64,1},im::Array{Float64,1},internal::Array{Int64,1})
            for i=1:length(internal)
                k=internal[i]
                $([curfields[j] for j=1:length(curfields)]...)
            end

            nothing
        end

        function init!(p::$(typeof(myprop)),v::Array{Float64,1},internal::Array{Int64,1})
            for i=1:length(internal)
                k=internal[i]
                $([initfields[j] for j=1:length(initfields)]...)
                p.p[i]=$(typeof(myprop.p))($(constfields...))
            end
        end        
    end    
end

const ena=50.0
const ek=-77.0
