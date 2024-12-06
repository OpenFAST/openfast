from __future__ import print_function
from rosco.toolbox.ofTools.fast_io.output_processing import output_processing
import rosco.toolbox

def FAST_IO_timeseries(fname):
    # interface to FAST_IO data load
    try:
        test = rosco.toolbox.__file__
    except:
        print('WARNING: rosco.toolbox required for wisdem.aeroelasticse.FAST_post.FAST_IO_timeseries')
    
    fast_out = output_processing()
    fast_data = fast_out.load_fast_out(fname, verbose=True)[0]
    return fast_data
