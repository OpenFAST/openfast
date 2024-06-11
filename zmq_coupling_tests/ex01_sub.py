'''
Example 1: Simple subscriber to FAST

This example shows how to create a simple subscriber that connects to a FAST to
get results at runtime.
'''
from zmq_python_toolbox.real_fast_interactor import RFInteractor
import os 
import multiprocessing as mp
import concurrent.futures

dir = os.path.dirname(__file__)
fstdir = os.path.join(dir, 'templateDir')


# Create a subscriber
sub = RFInteractor(ZmqInAddress=None, ZmqOutAddress="tcp://127.0.0.1:5557", verbose=True)


def run_openfast():
    return os.system('openfast ./templateDir/Main02.fst')

def run_sub():
    return sub.fast_sub(25)

def run_in_parallel():
    with concurrent.futures.ProcessPoolExecutor() as executor:
        p1 = executor.submit(run_openfast)
        p2 = executor.submit(run_sub)
        print(p1.result())
        print(p2.result())

if __name__ == '__main__':
    # create two processes, one for zmq sub and other for openfast
    run_in_parallel()
    
    
    