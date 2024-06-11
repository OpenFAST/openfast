import numpy as np

def yaml_read(filename=None, dictIn=None, lines=None, text=None):
    """
    read yaml files only supports:
       - Key value pairs: 
             key: value
       - Key with lists of lists:
             key:  
               - [0,1]
               - [0,1]
       - Comments are stripped based on first # found (in string or not)
       - Keys are found based on first : found (in string or not)
    """
    # --- swtich depending on what the user provided
    if filename is not None:
        # Read all lines at once
        with open(filename, 'r', errors="surrogateescape") as f:
            lines=f.read().splitlines()
    elif text is not None:
        lines = text.split('\n')
    elif lines is not None:
        # OK
        pass

    if dictIn is None:
        d=dict()
    else:
        d=dictIn

    # --- Loop on lines
    i=0
    while i<len(lines):
        l=_cleanComment(lines[i])
        i+=1;
        if len(l)==0:
            continue
        sp=l.split(':')
        if len(sp)==2 and len(sp[1].strip())==0:
            key=sp[0]
            array,i=_readDashList(i, lines)
            d[key]=array
        elif len(sp)==2:
            key=sp[0]
            val=sp[1].strip()
            #
            if val[0]=='[':
                d[key], _ = _readInlineList(val)
            else:
                try:
                    d[key]=int(val)
                except:
                    try:
                        d[key]=float(val)
                    except:
                        d[key]=val.strip() # string
        else:
            raise Exception('Line {:d} has colon, number of splits is {}, which is not supported'.format(len(sp)))
    return d

def _cleanComment(l):
    """ remove comments from a line"""
    return l.split('#')[0].strip()

def _readDashList(iStart, lines):
    """ """
    i=iStart
    # Count number of lines that starts with dash
    while i<len(lines):
        l = lines[i].strip()
        if len(l)==0:
            iEnd=i-1
            break
        if l[0]=='-':
            iEnd=i
            i+=1
        else:
            iEnd=i-1
            break
    nLines=iEnd-iStart+1

    # determine type based on first line
    firstLine = lines[iStart].lstrip()[1:]
    FirstElems, mytype= _readInlineList(firstLine)
    M = np.zeros((nLines,len(FirstElems)), mytype)
    if len(FirstElems)>0:
        for i in np.arange(iStart,iEnd+1):
            L,_ =  _readInlineList(lines[i].lstrip()[1:], mytype=mytype)
            M[i-iStart,:] = L
    return M, iEnd+1

def _readInlineList(line, mytype=None):
    """ 
    Parse a simple list of int, float or string
     [a, b, c]
     [a, b, c,]
    """
    L = _cleanComment(line.replace(']','').replace('[','')).strip().rstrip(',').strip()
    if len(L)==0:
        return np.array([]), float 
    L=L.split(',')
    L = np.asarray(L)
    if mytype is None:
        ## try to detect type
        try: 
            L=L.astype(int)
            mytype=int
        except:
            try: 
                L=L.astype(float)
                mytype=float
            except:
                try: 
                    L=L.astype(str)
                    mytype=str
                    L=np.array([c.strip() for c in L])
                except:
                    raise Exception('Cannot parse list from string: >{}<'.format(line)) 
    else:
        try:
            L = L.astype(mytype)
            if mytype==str:
                L=np.array([c.strip() for c in L])
        except:
            raise Exception('Cannot parse list of type {} from string: >{}<'.format(mytype, line)) 

    #
    return L, mytype


if __name__=='__main__':
#     #d=yaml_read('test.SD.sum.yaml')
    text = """
# Comment
IS1:     40567 # int scalar
FS1:     40567.32 # float scalar
FA1: # Array1
  - [  3.97887E+07,  0.00000E+00,  0.00000E+00]
  - [  0.00000E+00,  3.97887E+07,  0.00000E+00]
  - [  0.00000E+00,  0.00000E+00,  0.00000E+00]
FA2: # Array2
  - [  1.E+00,  0.E+00,  0.E+00,]
  - [  0.E+00,  1.E+00,  0.E+00,]
  - [  0.E+00,  0.E+00,  1.E+00,]
FL1: [       0.0,  0.0, 1.0 ] # FloatList1
FL2: [       0.0,  0.0, 1.0,] # FloatList2
SL2: [  aa, bb , cc, dd, ] #string list
"""
    text="""
EL1: [ ] #empty list
EL2: #empty list2
 - [ ] # empty list
"""
    d=yaml_read(text=text)
    print(d)
    for k,v, in d.items():
        if hasattr(v,'__len__'):
            if len(v)>0:
                print('{:12s} {:20s} {}'.format(k, str(type(v[0]))[6:], v[0]))
            else:
                print('{:12s} {:20s} {}'.format(k, str(type(v))[6:], v))
        else:
            print('{:12s} {:20s} {}'.format(k, str(type(v))[6:], v))
    


