
import Base: put!, wait, isready, take!



isready(D::MyChannel) = all(D.d)

function wait(D::MyChannel)
    while !isready(D)
        wait(D.cond_take)   
    end
end

function put!(D::MyChannel, k::Int64)
    @inbounds D.d[k] = true
    if isready(D)
        notify(D.cond_take)  
    else
        wait(D)
    end
    D
end

function take!(D::MyChannel)
    if isready(D)
        for i=1:length(D.d)
            @inbounds D.d[i]=false
        end
    end
end
