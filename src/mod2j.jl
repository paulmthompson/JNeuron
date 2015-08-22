
#this will construct types from mod files related to channels or whatever

#=

Everything in the mod files that is defined as a PARAMETER, STATE, or DERIVATIVE ends up being a part of param
Additionally, the last two entries in param seem to always be v, and _g (conductance)?
Things in param can be shared with other Props (like voltage is equal at every Prop for each node)

ppvar contains other values, looks like usually related to ions

So as far as Julia goes, at each node there are state variables unique to one ion channel(e.g. sodium gate m), state variables shared among ion channels (e.g. voltage), constant paramters for each channel
=#

type ParsedHoc
    file::Array{ByteString,1}
    

end

function ParsedHoc(filename::ASCIIString)
    f = open(filename)
    a=readlines(f)
    close(f)
    ParsedHoc()
end

function parsehoc(ph::ParsedHoc)
    linenum=1
    cursection=""
    sections=["COMMENT","UNITS","NEURON","PARAMETER","STATE","ASSIGNED","BREAKPOINT","INITIAL",
              "DERIVATIVE","PROCEDURE","FUNCTION"]

    while linenum<length(ph.file)

        if cursection==""
            for i in sections
                if contains(ph.file[linenum],i)
                    cursection=i
                    break
                end                
            end
            
        else
            if contains(ph.file[linenum],"}")
                cursection=""
            elseif contains(ph.file[linenum],"ENDCOMMENT")
                cursection=""
                #skip blank lines or lines with comments
            elseif cursection=="COMMENT"
            elseif cursection=="UNITS"
                #don't care for now
            elseif cursection=="NEURON"
                #word after READ should be collected to add to myvars dictionary
            elseif cursection=="PARAMETER"
                first=split(ph.file[linenum],"=")
                var=replace(first[1]," ", "")
                num=split(first[2], "(")
                num=replace(num[1], " ", "")
                parsehoc.param[var]=float(num)
            elseif cursection=="STATE"
                #get characters not equal to whitespace
            elseif cursection=="ASSIGNED"
                #get characters not in parenthesis
            elseif cursection=="BREAKPOINT"
                #second half of prop_calc method
            elseif cursection=="INITIAL"
                #prop_init method
            elseif cursection=="DERIVATIVE"
                #first half of prop_calc method
            elseif cursection=="PROCEDURE"
                #for this and function, collect whole thing, figure out what variables from where, and get rid of comments and delcarations
            elseif cursection=="FUNCTION"
            end
       
        end
        
    end
    
end
