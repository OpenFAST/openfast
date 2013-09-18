def read_input_file(fileName,init):
    f = open(fileName, 'r')

    line_offset = []
    offset = 0
    for line in f:
        line_offset.append(offset)
        offset += len(line)

    f.seek(0)

    i = 0
    for line in f:
        #print line
        words = line.split()
        if words[0] == "LineType":
            next(f)
            LineType_ref = i
        elif words[0] == "Node":
            next(f)
            Node_ref = i
        elif words[0] == "Element":
            next(f)
            Element_ref = i 
        elif words[0] == "Option":
            next(f)
            Option_ref = i 

        i+=1
    
    f.seek(line_offset[LineType_ref+2])
    
    for line in f:
        if line[0] == "-":
            break
        else:
            #print line
            init.SetCableLibraryData(line)

    f.seek(line_offset[Node_ref+3])
        
    for line in f:
        if line[0] == "-":
            break
        else:
            #print line
            init.SetNodeData(line)
            
    f.seek(line_offset[Element_ref+4])

    for line in f:
        if line[0] == "-":
            break
        else:
            #print line
            init.SetElementData(line)

    f.seek(line_offset[Option_ref+5])

    for line in f:
        if line[0] != "!":
            #print line
            init.SetSolverOptions(line)
    
    

if __name__ == '__main__':  
    import sys
    import MAP
    import numpy as np

    #map_file = "../Test_Cases/input12.map"
    map_file = sys.argv[1];
    
    Init = MAP.MAP_InitInputType()
    InitOut = MAP.MAP_InitOutputType()

    read_input_file( map_file , Init )
    Init.SetDepth     ( -200 );
    Init.SetGravity   ( 9.81 );
    Init.SetSeaDensity( 1025 );
    
    d   = MAP.MAP_OtherStateType()
    u   = MAP.MAP_InputType()
    p   = MAP.MAP_ParameterType()
    z   = MAP.MAP_ConstraintStateType()
    y   = MAP.MAP_OutputType()
    msg = MAP.MAP_Message()
    err = MAP.MAP_ErrStat()

    MAP.MSQS_Init(Init,u,p,None,None,z,d,y,None,InitOut,err,msg)
    if err.error_status() != 0 :
        print msg.status()
    
	print "here is the converge reason : " #msg.converge_reason()
    	
    # run every time step
    MAP.MSQS_UpdateStates(0,0,u,p,None,None,z,d,err,msg)
    if int(err.error_status()) !=0 :
        print msg.status()
	 
    MAP.MSQS_CalcOutput( 0 , u,p,None,None,z,d,y,err,msg)
    if err.error_status() != 0 :
        print msg.status()

    #print d.summary()
    #d.plot( err,msg )

    # Destroy objects
    MAP.MSQS_End(u,p,None,None,z,d,y,err,msg)
    if err.error_status() != 0 :
        print msg.status()		
    
    
	
    #u_new = u.get() 
    #u.set( u_new )
    #print u.details()
