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
    
