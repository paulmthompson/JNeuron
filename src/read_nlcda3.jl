
proc init() {
	quiet = 0
	debug_on = 0
	gm = new GUIMath()
	filetype = "Neurolucida V3"
	current = new Import3d_LexToken()
	look_ahead = new Import3d_LexToken()
	look_ahead2 = new Import3d_LexToken()
	eof=0
	number=1  leftpar=2  rightpar=3  comma=4  bar=5
	set=6  rgb=7  string=8  label_=9  err_=10
	leftsp=11  rightsp=12
	tokens = new List()
	tokensappend("eof", "number", "leftpar", "rightpar", "comma", "bar")
	tokensappend("set", "rgb", "string", "label", "err")
	tokensappend("leftsp", "rightsp")
	plist = new List()
}

const markers=["Dot","OpenStar","FilledQuadStar","CircleArrow","OpenCircle","DoubleCircle",
         "OpenQuadStar","CircleCross","Cross","Circle1","Flower3","Plus","Circle2",
         "Pinwheel","OpenUpTriangle","Circle3","TexacoStar","OpenDownTriangle","Circle4",
         "ShadedStar","OpenSquare","Circle5","SkiBasket","Asterisk","Circle6","Clock",
         "OpenDiamond","Circle7","ThinArrow","FilledStar","Circle8","ThickArrow","FilledCircle",
         "Circle9","SquareGunSight","FilledUpTriangle","Flower2","GunSight","FilledDownTriangle",
         "SnowFlake","TriStar","FilledSquare","OpenFinial","NinjaStar","FilledDiamond",
         "FilledFinial","KnightsCross","Flower","MalteseCross","Splat"]

const nonsense=["Name", "ImageCoords","Thumbnail","Color","Sections","SSM","dZI","Normal",
                "Low","High","Generated","Incomplete","SSM2","Resolution"]

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

#append()
function init_sec(sec::Section,first,start,n,x,y,z,d)
    sec.raw[:,1]=deepcopy(x)
    sec.raw[:,2]=deepcopy(y)
    sec.raw[:,3]=deepcopy(z)
    sec.d[:]=deepcopy(d)
    nothing 
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
    state=zeros(Int64,1)
    while linenum < length(nlcda.file)
        leftpar=length(matchall(r"\(",nlcda.file[linenum]))
        rightpar=length(matchall(r"\)",nlcda.file[linenum]))
        depth+=(leftpar-rightpar)
        if depth>=length(state)
            append!(state,0)
        end
        
        if leftpar>0
            if contains(nlcda.file[linenum],"CellBody")
                state[depth+1]=1
            elseif contains(nlcda.file[linenum],"Axon")
                state[depth+1]=2
            elseif contains(nlcda.file[linenum],"Dendrite")
                state[depth+1]=3
            elseif contains(nlcda.file[linenum],"Apical")
                state[depth+1]=4
            elseif markerdetect(nlcda, linenum)
                state[depth+1]=0
            elseif nonsensedetect(nlcda,linenum)
                state[depth+1]=0
            else
                mynums=matchall(r"[-+]?[0-9]*\.?[0-9]+", nlcda.file[linenum])
                if state[depth+1]>0 & length(mynums)>3
                    
                elseif
                end
                
            end
            
        elseif rightpar>0
            #close current section and go back to last one at previous depth
        else
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

    


