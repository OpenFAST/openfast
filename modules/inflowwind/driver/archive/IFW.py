'''
This is the main python script that calls inflowWind (IFW) for testing

'''

# Import modules
import os
import numpy as np
from ctypes import (
	CDLL,
	POINTER,
	pointer,
	c_bool,
	c_char,
	c_char_p,
    c_wchar_p, 
	c_int,
	c_float,
	c_double,
	c_longdouble,
	create_string_buffer,
	byref,
	get_errno,
	sizeof
)

# Run inflow wind - via the stand-alone driver
os.system("inflowwind_driver -tsteps[10] IFW.dvr")
# calls IEA-15-240-RWT-UMaineSemi_InflowFile.dat input file
# outputs data in file named "IFW.WindGrid.out" 
