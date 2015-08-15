
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

#=
sec type
Cellbody=1
Axon =2
Dendrite =3
Apical = 4
=#

#=
Find the section type

For trees
read in the x, y, z values and diameter for point in that tree
cursec = Section(size(x) before treepoints are read, size(x)_after - size(x)_before)
=#

type Import3d_LexToken
    s
    token
    x
    itok
    iline
end

#=
Need data container that is keeping track of parent sec, current sec, x,y,z,d values etc
=#

type Section
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
    mytype::Int64
    centroid_color::Int64
    id::Int64
    raw::Array{Int64,2}
    xyx::Array{Int64,2}
    d::Array{Int64,1}
end

type nlcda3
    cursec::Section
    parentsec::Section
    sections::Array{Section,1}
    file::Array{ByteString,1}
end

function nlcda3(filename::ASCIIString)
    f = open(filename)
    a=readlines(f)
    close(f)
    nlcda3()
end

function Section(ID::Int64,N::Int64)
    Section(0,0,0,0,0,1,0,0,-1,0,0,2,ID,Array(Int64,N,3),Array(Int64,N,3),Array(Int64,N))         
end

function newsec(start::Int64,finish::Int64)
    first = 0
    n = finish - start
    cursec(start,n)
    cursec.mytype=sectype
    append!(mytype,sectype)
    append!(sections,cursec)
    init_sec(cursec,first,start,n,x,y,z,d)
end

function input(morphology::ASCIISTRING)
    b2serr = new List()
    b2sinfo = new List()
    nspine = 0
    err = 0
    mytype = new Vector() #originally called type
    sections = Array{Section,0)
    alloc(25000, x, y, z, d, iline)
    lines = new List(25000)
    itoken = 0
    depth = 0
    rdfile(morphology)
    firstpoints = new Vector(sections.count)
    set_firstpoints()
    connect2soma()
end

function parse_file(nlcda::nlcda3)
    linenum=1
    leftpar=0
    rightpart=0
    depth=0
    skip=0
    state=0
    opensec=Array(Section,1)
    cursec=Section(0)
    while linenum < length(nlcda.file)
        leftpar=length(matchall(r"\(",nlcda.file[linenum]))
        rightpar=length(matchall(r"\)",nlcda.file[linenum]))       
        if leftpar>0
            if contains(nlcda.file[linenum],"CellBody")
                state=1
                append!(nlcda.sections,Section(1))
                #must be on first level. make top of open sec, set as current sec, and no parent
            elseif contains(nlcda.file[linenum],"Axon")
                state=2
                append!(nlcda.sections,Section(2))
            elseif contains(nlcda.file[linenum],"Dendrite")
                state=3
                append!(nlcda.sections,Section(3))
            elseif contains(nlcda.file[linenum],"Apical")
                state=4
                append!(nlcda.sections,Section(4))
            elseif markerdetect(nlcda, linenum)
                skip+=1
            elseif nonsensedetect(nlcda,linenum)
            else
                if (leftpar-rightpar)>0 & state>0 #check for new section
                    depth+=(leftpar-rightpar)
                    append!(nlcda.sections,Section(state))
                    #add to current sec, set previous depth as parent, add to open sec
                else
                    
                    mynums=matchall(r"[-+]?[0-9]*\.?[0-9]+", nlcda.file[linenum])
                    if state>0 & length(mynums)>3 & skip<1
                        #add numbers to current section
                    end
                end
                
            end
            
        elseif rightpar>0
            if skip>0
                skip-=1
            else
            end
            #close current section and go back to last one at previous depth
        else #no parenthesis, so can just ignore
        end
            linenum+=1
    end      
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

    


