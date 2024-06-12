""""
Test 01
=======
Testing message integrity of the communication between the real-time interactor and the FAST model.


"""
import os
import pickle 
import gzip 
import numpy as np
import matplotlib.pyplot as plt

from openfast_toolbox.io import FASTOutputFile

def test01_integrity(path_zmq_out: str, path_fast_out: str):
    
    # Load zmq output
    with gzip.open(path_zmq_out, 'rb') as f: 
        zmq_data = pickle.load(f)
        
    # Load FAST output
    fast_data = FASTOutputFile(path_fast_out).toDataFrame()
    
    # # Compare
    # assert len(zmq_data) == len(fast_data)
    
    outputs = fast_data.columns
    
    for output in outputs:
        output_zmq = ' ' + output.split('_')[0]
        assert np.allclose(zmq_data[output_zmq], fast_data[output], atol=1e-6)
        
    print('Test 01 passed.')
    
    return True

if __name__ == '__main__':
        cwd = os.getcwd()
        path_zmq_out = os.path.join(cwd, '../zmq_logs/OFZMQ_test01.pkl.gz')
        path_fast_out = os.path.join(cwd, '../../templateDir/Main01.outb')
        test01_integrity(path_zmq_out, path_fast_out)
    