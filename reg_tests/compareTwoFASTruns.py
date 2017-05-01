import fast_io
import numpy as np
import sys, os

if __name__=="__main__":
    fin1 = sys.argv[1]
    fin2 = sys.argv[2]
    maxTol = sys.argv[3]

    try:
        d1, i1 = fast_io.load_binary_output(fin1)
    except:
        print "Problem reading file ", fin1

    try:
        d2, i2 = fast_io.load_binary_output(fin2)
    except:
        print "Problem reading file ", fin2

    nCols = np.size(d1,1)
    var = np.ones(nCols)
    for j in range(nCols):
        var[j] = (d1[:,j]-d2[:,j]).var()

    if (max(var) < maxTol):
        sys.exit(0)
    else:
        print "var = ", var
        print i1['attribute_names']
        sys.exit(1)
