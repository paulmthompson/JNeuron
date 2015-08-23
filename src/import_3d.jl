
#=
These functions and types are designed to import 3D files of various formats. Currently only neurolucida3 is supported.

This will mostly address the functionality that was available in Neuron's import3d hoc files. The Import3d had lots of GUI stuff intertwined in it, though and I'm not planning on doing anything GUI related anytime soon, but I'd like to keep it separate if I do
=#

type Section3D
    parent::Array{Section3D,1}
    children::Array{Section3D,1}
    childind::Array{Int64,1}
    mytype::Int64 #Cellbody=1,Axon=2,Dendrite=3,Apical=4
    raw::Array{Float64,2}
    d::Array{Float64,1}
end

function Section3D(ID::Int64)
    Section3D(Array(Section3D,0),Array(Section3D,0),Array(Int64,0),ID,Array(Float64,0,3),Array(Float64,0))         
end

type curxyz
    x::Array{Float64,1}
    y::Array{Float64,1}
    z::Array{Float64,1}
end

function curxyz()
    curxyz(Array(Float64,0),Array(Float64,0),Array(Float64,0))
end

abstract Import3D

function input(morphology::ASCIIString)

    if contains(morphology,".asc")
        import3d=nlcda3(morphology)
    end
        
    parse_file(import3d)

    return import3d
    #connect2soma()

    #should then "instantiate" to create Neuron object and return it
end

function connect2soma(import3d::Import3D)

    #combine somas with overlapping bounding boxes

    #find sections that arent' somas and don't have parents and label them as roots

    #loop through each soma, and find what roots connect to it
    #if inside
    #parentsec=soma section
    #parentx=.5 #where along section of parent it connects
    #somehow incorporate center of soma into root , i think as starting point
    #first=1
    #fid=1

    #If roots are not within any soma bounding box, connect to the closest one   
end

function instantiate(import3d::Import3D)

    neuron=Neuron()

    for i=1:length(import3d.sections)
        sec=Section(import3d.sections[i])
        add_sec(neuron,sec)
    end

    #grow section children
    for i=1:length(import3d.sections)
        push!(neuron.seclist[import3d.sections[i].parent],neuron.seclist[i])
    end
    
    
    #connect them, making adjustments to 3d points as necessary
end

#=
Neurolucida file types
=#

type nlcda3 <: Import3D
    sections::Array{Section3D,1}
    mytypes::Array{Int64,1} # 4x1 array to tally total num of each section type
    file::Array{ByteString,1}
    curxyz::curxyz
    opensecs::Array{Int64,1} #inds of last section at that depth
    depth::Int64
end

function nlcda3(filename::ASCIIString)
    f = open(filename)
    a=readlines(f)
    close(f)
    nlcda3(Array(Section3D,0),zeros(Int64,4),a,curxyz(),zeros(Int64,0),0) 
end

const markers=["Dot","OpenStar","FilledQuadStar","CircleArrow","OpenCircle","DoubleCircle",
         "OpenQuadStar","CircleCross","Cross","Circle1","Flower3","Plus","Circle2",
         "Pinwheel","OpenUpTriangle","Circle3","TexacoStar","OpenDownTriangle","Circle4",
         "ShadedStar","OpenSquare","Circle5","SkiBasket","Asterisk","Circle6","Clock",
         "OpenDiamond","Circle7","ThinArrow","FilledStar","Circle8","ThickArrow","FilledCircle",
         "Circle9","SquareGunSight","FilledUpTriangle","Flower2","GunSight","FilledDownTriangle",
         "SnowFlake","TriStar","FilledSquare","OpenFinial","NinjaStar","FilledDiamond",
         "FilledFinial","KnightsCross","Flower","MalteseCross","Splat"]

const nonsense=["Name", "ImageCoords","Thumbnail","Color","Sections","SSM","dZI","Normal",
                "Low","High","Generated","Incomplete","SSM2","Resolution","\"CellBody\""]

function parse_file{T<:nlcda3}(nlcda::T)
    linenum=1
    leftpar=0
    rightpart=0
    skip=0
    state=0
    while linenum < length(nlcda.file)
        leftpar=length(matchall(r"\(",nlcda.file[linenum]))
        rightpar=length(matchall(r"\)",nlcda.file[linenum]))       
        if leftpar>0
            if contains(nlcda.file[linenum],"(CellBody)")
                state=1               
                newsec(nlcda,state)
                newparent(nlcda)
            elseif contains(nlcda.file[linenum],"(Axon)")
                state=2
                newsec(nlcda,state)
                newparent(nlcda)
            elseif contains(nlcda.file[linenum],"(Dendrite)")
                state=3
                newsec(nlcda,state)
                newparent(nlcda)
            elseif contains(nlcda.file[linenum],"(Apical)")
                state=4
                newsec(nlcda,state)
                newparent(nlcda)
            elseif markerdetect(nlcda, linenum)
                skip+=1
            elseif nonsensedetect(nlcda,linenum)
            else
                if ((leftpar-rightpar)>0) & (state>0) #check for new section
                    nlcda.depth+=1
                    newsec(nlcda,state)
                    newchild(nlcda)
                else        
                    mynums=matchall(r"[-+]?[0-9]*\.?[0-9]+", nlcda.file[linenum])
                    if (state>0) & (length(mynums)>3) & (skip<1)
                        dimadd(nlcda,mynums)
                    end
                end
                
            end
            
        elseif rightpar>0
            if skip>0
                skip-=1
            else
                closesec(nlcda)
            end
        elseif contains(nlcda.file[linenum], " |")
            closesec(nlcda)
            nlcda.depth+=1
            newsec(nlcda,state)
            newchild(nlcda)
        else #no parenthesis so just ignore
            
        end
            linenum+=1
    end
    nothing
end

function newsec(nlcda::nlcda3,state::Int64)
    push!(nlcda.sections,Section3D(state)) #add new section to overall list of 3D
    nlcda.mytypes[state]+=1
    nothing
end

function newparent(nlcda::nlcda3)
    nlcda.curxyz=curxyz()
    nlcda.depth+=1
    push!(nlcda.opensecs,length(nlcda.sections))
    nothing
end

function newchild(nlcda::nlcda3) #new branch off of main section
    nlcda.sections[end-1].raw=vcat(nlcda.sections[end-1].raw,[nlcda.curxyz.x nlcda.curxyz.y nlcda.curxyz.z]) #add points that have accumulated to most recent section
    nlcda.sections[end].raw=nlcda.sections[nlcda.opensecs[nlcda.depth-1]].raw[end,:] #make first point equal to most recent data point (where branch was)
    nlcda.sections[end].d=[nlcda.sections[nlcda.opensecs[nlcda.depth-1]].d[end]]

    push!(nlcda.sections[nlcda.opensecs[nlcda.depth-1]].children,nlcda.sections[end]) #add as child to parent
    push!(nlcda.sections[end].parent, nlcda.sections[nlcda.opensecs[nlcda.depth-1]]) #add as parent to child
    
    nlcda.curxyz=curxyz() #reset xyz container
    push!(nlcda.opensecs,length(nlcda.sections)) #add new section to list of open

    nothing
end

function closesec(nlcda::nlcda3)
    nlcda.sections[end].raw=vcat(nlcda.sections[end].raw,[nlcda.curxyz.x nlcda.curxyz.y nlcda.curxyz.z])
    nlcda.depth-=1
    nlcda.curxyz=curxyz()
    pop!(nlcda.opensecs)

    nothing
end

function dimadd(nlcda::nlcda3,mynums::Array{SubString{UTF8String},1})
    push!(nlcda.curxyz.x,float(mynums[1]))
    push!(nlcda.curxyz.y,float(mynums[2]))
    push!(nlcda.curxyz.z,float(mynums[3]))
    push!(nlcda.sections[end].d,float(mynums[4]))
    nothing
end

function markerdetect(nlcda::nlcda3, linenum::Int64)

    for mark in markers
        if contains(nlcda.file[linenum],mark)
            return true
        end
    end

    return false
end

function nonsensedetect(nlcda::nlcda3,linenum::Int64)
    for non in nonsense
        if contains(nlcda.file[linenum],non)
            return true
        end
    end

    return false
end


