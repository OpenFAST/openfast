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
            init.setCableLibraryData(line)

    f.seek(line_offset[Node_ref+3])
        
    for line in f:
        if line[0] == "-":
            break
        else:
            #print line
            init.setNodeData(line)
            
    f.seek(line_offset[Element_ref+4])

    for line in f:
        if line[0] == "-":
            break
        else:
            #print line
            init.setElementData(line)

    f.seek(line_offset[Option_ref+5])

    for line in f:
        if line[0] != "!":
            #print line
            init.setSolverOptions(line)
    
    

if __name__ == '__main__':  
    import sys
    import MAP
    import numpy as np
    
    map_file = "../MAP_Fortran_Binding/input6.map"
    #map_file = sys.argv[1];
    
    InitIn = MAP.MAP_InitInputType()
    InitOut  = MAP.MAP_InitOutputType()

    read_input_file( map_file , InitIn )
    InitIn.setDepth     ( "-350" );
    InitIn.setGravity   ( "9.81" );
    InitIn.setSeaDensity( "1020" );
    
    d   = MAP.MAP_OtherStateType()
    u   = MAP.MAP_InputType()
    p   = MAP.MAP_ParameterType()
    z   = MAP.MAP_ConstraintStateType()
    y   = MAP.MAP_OutputType()
    msg = MAP.MAP_Message()
    err = MAP.MAP_ErrStat()

#    print 'Running MAP v2'

    MAP.MSQS_Init(InitIn,u,p,None,None,z,d,y,None,InitOut,err,msg)
    if err.error_status() != 0 :
        print msg.status()
    
    InitIn  = None;
    InitOut = None;
        
    MAP.MSQS_CalcOutput( 0 , u,p,None,None,z,d,y,err,msg)
    if err.error_status() != 0 :
        print msg.status()
    
    # run every time step
    MAP.MSQS_UpdateStates(0,0,u,p,None,None,z,d,err,msg)
    if int(err.error_status()) !=0 :
        print msg.status()

#    print d.summary()
    d.plot( err,msg )

#    # Destroy objects
#    MAP.MSQS_End(u,p,None,None,z,d,y,err,msg)
#    if err.error_status() != 0 :
#        print msg.status()
