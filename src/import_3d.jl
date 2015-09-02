
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

type nlcda3 <: Import3D
    sections::Array{Section3D,1}
    mytypes::Array{Int64,1} # 4x1 array to tally total num of each section type
    file::Array{ByteString,1}
    curxyz::curxyz
    opensecs::Array{Int64,1} #inds of last section at that depth
    depth::Int64
end

function input(morphology::ASCIIString)

    if contains(morphology,".asc")
        import3d=nlcda3(morphology)
    end
        
    parse_file(import3d)

    return import3d

    #should then "instantiate" to create Neuron object and return it
end

function instantiate(import3d::Import3D)

    neuron=Neuron()

    somaind=Array(Int64,0)
    rootind=Array(Int64,0)
    childinds=Array(Array{Int64,1},0)
    mapping=Dict{Int64,Int64}(0=>0)
    #skip somas, we'll do them later
    for i=1:length(import3d.sections)
        if import3d.sections[i].mytype>1
            if length(import3d.sections[i].parent)==0
                push!(rootind,i)
            end
            sec=Section(import3d.sections[i])
            add_sec(neuron,sec)
            push!(childinds,import3d.sections[i].childind)
            mapping[i]=length(neuron.secstack)
        elseif import3d.sections[i].mytype==1
            push!(somaind,i)
        end
    end
    
    connect2soma(import3d,neuron,somaind,rootind)

    for i=1:length(childinds)
        for j in childinds[i]
            push!(neuron.secstack[i].child,neuron.secstack[mapping[j]])
        end
    end
    
    for i=1:length(neuron.secstack)  
        neuron.secstack[i].refcount=i
    end
    
    neuron
    #connect them, making adjustments to 3d points as necessary
end

function connect2soma(import3d::Import3D,neuron::Neuron,somaind::Array{Int64,1},rootind::Array{Int64,1})

    #We need to
    #1) determine how many somas are described in file
    #2) assign a sensical 3D shape to it
    #3) find the centroid of that 3D shape

    centroids=Array(Float64,3,0)
    somas=Array(Section,0)
    #Determine how soma is described in this file

    if reduce(&, [mostly_constantz(import3d,somaind[i]) for i=1:length(somaind)])
        #if soma contours have a constant z

        if length(somaind)==1
            #approximate as sphere with diameter taken from section

            centroids=hcat(centroids,sphere_approx(import3d, somaind[1],somas))
            
        else
            overlaps=stack_overlap(import3d, somaind) #find which soma sections overlap
            for i=1:length(overlaps)
                if length(overlaps[i])==1 #if a soma doesn't overlap with any other sections
                    
                    #approximate as sphere
                    centroids=hcat(centroids,sphere_approx(import3d, overlaps[i],somas))
                else
                    #multiple outlines of cell at different z positions - "stack of pancakes"
                    #find principle axis through stack
                    #approximate as "carrot" looking thing along principle axis
                    #need to make new soma Section
                end
            end
        end

    else
        #xyzd measurements along estimated centroid - "stack of frusta"
    end

    #calculate root distance to each centroid
    #connect to closest centroid (first point in root)
    #add root as child to soma
    ds=zeros(Float64,length(somas))
    for i in rootind
        for j=1:length(somas)
            ds[j]=sqrt((centroids[1,j]-neuron.secstack[i].pt3d[1].x)^2 + (centroids[2,j]-neuron.secstack[i].pt3d[1].y)^2 + (centroids[3,j]-neuron.secstack[i].pt3d[1].z)^2)
        end
        push!(neuron.secstack[i].child,somas[indmin(ds)]) #add root as child of soma
    end

    append!(neuron.secstack,somas)

    nothing
    
end

function mostly_constantz(import3d::Import3D, ind::Int64)
    if std(import3d.sections[ind].raw[:,3])/mean(import3d.sections[ind].raw[:,3])<0.1
        return true
    else
        return false
    end
end

function stack_overlap(import3d::Import3D, somaind::Array{Int64,1})
    #this is the ugliest thing i have ever created. shame on me
 
    n = length(somaind)
    combos=[[list[c[i]-i+1] for i=1:length(c)] for c in combinations(collect(1:(n+k-1)),2)]

    overlaps=Array(Array{Int64,1},1)
    overlaps[1]=[1]
    for i=1:length(combos)
        bb1=boundbox(import3d, combos[i][1])
        bb2=boundbox(import3d, combos[i][2])

        if max(bb1[1],bb2[1])<min(bb1[2],bb2[2])
            if max(bb1[3],bb2[3])<min(bb1[4],bb2[4]) #if bounding boxes of these slices overlap
                flag=false
                for j=1:length(overlaps) #check groups of overlapping slices from previous iterations
                    if length(intersect(overlaps[j],combos[i]))>0
                        append!(overlaps[j],setdiff(combos[i],a[j]))
                        flag=true
                    end
                end
                if flag==false #if no group has yet been determined, assign to new group
                    append!(overlaps,combos[i])
                end
            end
        end       
    end

    return overlaps
   
end

function boundbox(import3d::Import3D, ind::Int64)

    xmin=min(import3d.sections[ind].raw[:,1])
    xmax=max(import3d.sections[ind].raw[:,1])
    ymin=min(import3d.sections[ind].raw[:,2])
    ymax=max(import3d.sections[ind].raw[:,2])

    bb=collect(xmin, xmax, ymin, ymax)
    
end

function sphere_approx(import3d::Import3D, ind::Int64,somas::Array{Section,1})
    #going to make 3 3D points to describe 1 constant z section through the soma as a cylinder with L=D
    #trace perimeter
    perimeter=0.0
    for i=2:size(import3d.sections[ind].raw,1)
        perimeter+=sqrt((import3d.sections[ind].raw[i,1]-import3d.sections[ind].raw[i-1,1])^2 +
                        (import3d.sections[ind].raw[i,2]-import3d.sections[ind].raw[i-1,1])^2)
    end

    perimeter+=sqrt((import3d.sections[ind].raw[1,1]-import3d.sections[ind].raw[end,1])^2 +
                    (import3d.sections[ind].raw[1,2]-import3d.sections[ind].raw[end,1])^2)

    radi=perimeter/(2*pi)

    centroid=[mean(import3d.sections[ind].raw[:,1]); mean(import3d.sections[ind].raw[:,2]); mean(import3d.sections[ind].raw[:,3])]

    mypoints=[Pt3d(centroid[1]-radi,centroid[2]-radi,centroid[3],2*radi,0.0);
    Pt3d(centroid...,2*radi,.5); Pt3d(centroid[1]+radi,centroid[2]+radi,centroid[3],2*radi,1.0)]


    sec=Section(1,1,Array(Node,0),Array(Section,0),mypoints,0.0,0.0,2*radi)
    push!(somas,sec)
    
    centroid
    
end

                   
#=
Neurolucida file types
=#

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
    push!(nlcda.sections[nlcda.opensecs[nlcda.depth-1]].childind,length(nlcda.sections))
    
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


function Section(section3d::Section3D) #like new_section
    sec=Section(1,section3d.mytype,Array(Node,0),Array(Section,0),Array(Pt3d,size(section3d.raw,1)),0.0,0.0,0.0)

    #add 3d points from 3d
    for i=1:length(section3d.d)
        if i>1
            mylength=sqrt((section3d.raw[i,1]-section3d.raw[i-1,1])^2+(section3d.raw[i,2]-section3d.raw[i-1,2])^2+(section3d.raw[i,3]-section3d.raw[i-1,3])^2)
            sec.pt3d[i]=Pt3d(section3d.raw[i,:]...,section3d.d[i],mylength)
            sec.length+=sec.pt3d[i].arc
        else
            sec.pt3d[i]=Pt3d(section3d.raw[i,:]...,section3d.d[i],0.0)
        end
    end

    #normalize arc length
    for i=2:size(sec.pt3d,1)
        sec.pt3d[i].arc=sec.pt3d[i].arc/sec.length
    end

    sec.parentx=1.0
    
    sec
end
