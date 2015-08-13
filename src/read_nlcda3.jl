
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


type Import3d_LexToken
    s
    token
    x
    itok
    iline
end


function input(morphology::ASCIISTRING)
    b2serr = new List()
    b2sinfo = new List()
    nspine = 0
    err = 0
    mtype = new Vector() #originally called type
    sections = Array{something,1000) #don't know why they said 1000. doesn't matter
    alloc(25000, x, y, z, d, iline)
    lines = new List(25000)
    itoken = 0
    depth = 0
    rdfile(morphology) #seems important
    firstpoints = new Vector(sections.count)
    set_firstpoints()
    connect2soma()
end

#sets firstpoints equal to object
function set_firstpoints(sections::Array{something,1000))
    
end

function rdfile(morphology::ASCIISTRING)
    iline_ = 0
    file = new File(morphology)

    #finds length of file (i=number of lines-1)
    for (i=0; !file.eof(); i += 1) {
	file.gets(line)
	}
    end
    
    alloc(i, x, y, z, d, iline) #make vectors the length of lines in the file
    file.close
    
    lines = new List(25000)
    line=""
    if (!quiet) {
	printf("\n")
	}
    end
    
    file.ropen()
    p_file()
    file.close
            
end
            
function p_file() {
	look_ahead2.token = eof
	look_ahead.token = eof
	if (lex(current) != eof) {
		if (lex(look_ahead) != eof) {
			lex(look_ahead2)
		}
	            }
        end
end

	enter("p_file")
	objects()
	leave("p_file")
            }
end

proc objects() {
	enter("objects")
	object()
	while(1) {
		optionalcomma()
		if (current.token != leftpar) {
			break
		}
		object()
	}
	leave("objects")
            }


proc object() {local i
	i = current.itok
	enter("object")
	if (current.token == leftpar) {
		plook()		
		if (look_ahead.token == string) {
			contour()
		}else if (look_ahead.token == label_) {
			marker_or_property()
		}else if (look_ahead.token == leftpar) {
			tree_or_text()
		}else if (look_ahead.token == set) {
			p_set()
		}else{
			p_err()
		}
	}else{
		p_err()
	}
	leave("object")
	if (i == current.itok) {
		print "internal error: ", "object consumed no tokens"
		stop
	}
}
