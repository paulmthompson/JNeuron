
#=
These functions and types are designed to import 3D files of various formats. Currently only neurolucida3 is supported.

This will mostly address the functionality that was available in Neuron's import3d hoc files. The Import3d had lots of GUI stuff intertwined in it, though and I'm not planning on doing anything GUI related anytime soon, but I'd like to keep it separate if I do
=#

type Section3D
    is_subsidiary::Bool
    ztrans::Int64
    first::Bool
    fid::Int64
    nameindex::Float64
    parentx::Float64
    volatile::Bool
    volatile2::Bool
    pid::Int64
    iscontour::Bool
    mytype::Int64 #Cellbody=1,Axon=2,Dendrite=3,Apical=4
    centroid_color::Int64
    id::Int64
    raw::Array{Int64,2}
    xyz::Array{Int64,2}
    d::Array{Int64,1}
end

function Section3D(ID::Int64)
    Section(0,0,0,0,0,1,0,0,-1,0,0,2,ID,Array(Int64,0,3),Array(Int64,0,3),Array(Int64,0))         
end

type curxyz #may need to keep this whole list. I think neuron might do that
    x::Array{Float64,1}
    y::Array{Float64,1}
    z::Array{Float64,1}
end

function curxyz()
    curxyz(Array(Float64,0),Array(Float64,0),Array(Float64,0))
end

abstract Import3D

function input(morphology::ASCIISTRING)

    if contains(morphology,".asc")
        import3d=nlcda3(morphology)
    end
        
    parse_file(import3d)
    
    #firstpoints = new Vector(sections.count)
    #set_firstpoints() #don't know how this is different than mytype
    connect2soma()

    #should then "instantiate" to create Neuron object and return it
end


function connect2soma{T<:Import3D}(Import3D::T)
    #move somas to beginning of list

    #combine somas with overlapping bounding boxes

    #find sections that arent' somas and don't have parents and label them as roots

    #loop through each soma, and find what roots connect to it
    #if inside
    #parentsec=soma section
    #parentx=.5
    #somehow incorporate center of soma into root , i think as starting point
    #first=1
    #fid=1

    #If roots are not within any soma bounding box, connect to the closest one   
end

function instantiate(import3d::Import3D)

    neuron=Neuron()

    for i=1:length(import3d.sections)
        newsec=Section(import3d.sections[i])
        add_sec(neuron,newsec)
    end
    

    #connect them, making adjustments to 3d points as necessary
end

#=
Neurolucida file types
=#

type nlcda3 <: Import3D
    cursec::Section3D
    parentsec::Section3D
    sections::Array{Section3D,1}
    mytypes::Array{Int64,1} # 4x1 array to tally total num of each section type
    file::Array{ByteString,1}
    opensec::Array{Section3D,1}
    curxyz::curxyz
end

function nlcda3(filename::ASCIIString)
    f = open(filename)
    a=readlines(f)
    close(f)
    nlcda3() #initialize 
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
    depth=0
    skip=0
    state=0
    while linenum < length(nlcda.file)
        leftpar=length(matchall(r"\(",nlcda.file[linenum]))
        rightpar=length(matchall(r"\)",nlcda.file[linenum]))       
        if leftpar>0
            if contains(nlcda.file[linenum],"CellBody")
                state=1
                newsec(nlcda,state)
                parentsec(nlcda)
            elseif contains(nlcda.file[linenum],"Axon")
                state=2
                newsec(nlcda,state)
                parentsec(nlcda)
            elseif contains(nlcda.file[linenum],"Dendrite")
                state=3
                newsec(nlcda,state)
                parentsec(nlcda)
            elseif contains(nlcda.file[linenum],"Apical")
                state=4
                newsec(nlcda,state)
                parentsec(nlcda)
            elseif markerdetect(nlcda, linenum)
                skip+=1
            elseif nonsensedetect(nlcda,linenum)
            else
                if (leftpar-rightpar)>0 & state>0 #check for new section
                    depth+=(leftpar-rightpar)
                    newsec(nlcda,state)
                    newchild(nlcda)
                else        
                    mynums=matchall(r"[-+]?[0-9]*\.?[0-9]+", nlcda.file[linenum])
                    if state>0 & length(mynums)>3 & skip<1
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
        else #no parenthesis, so can just ignore
        end
            linenum+=1
    end
    nothing
end

function newsec(nlcda::nlcda3,state::Int64)
    append!(nlcda.sections,Section3D(state))
    nlcda.mytype[state]+=1
    nlcda.cursec=nlcda.sections[end]
    nothing
end

function newparent(nlcda::nlcda3)
    nlcda.opensec=[nlcda.sections[end]] #need to have no parent
    nlcda.curxyz=curxyz()
    nothing
end

function newchild(nlcda::nlcda3)
    nlcda.opensec[end-1].raw=vcat(nlcda.opensec[end-1].raw,
[nlcda.curxyz.x;nlcda.curxyz.y;nlcda.curxyz.z])
    nlcda.curxyz=curxyz()
    append!(nlcda.opensec,nlcda.sections[end])
    nlcda.sections[end].parent=length(nlcda.sections)-1
    nothing
end

function closesec(nlcda::nlcda3)
    nlcda.cursec.raw=vcat(nlcda.cursec.raw,[nlcda.curxyz.x;nlcda.curxyz.y;nlcda.curxyz.z])
    nlcda.curxyx()
    pop!(nlcda.opensec)
    nlcda.cursec=nlcda.opensec[end]
    nothing
end

function dimadd(nlcda::nlcda3,mynums::Array{SubString{UTF8String},1})
    append!(nlcda.curxyz.x,float(mynums[1]))
    append!(nlcda.curxyz.y,float(mynums[2]))
    append!(nlcda.curxyz.z,float(mynums[3]))
    append!(nlcda.cursec.d,float(mynums[4]))
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


