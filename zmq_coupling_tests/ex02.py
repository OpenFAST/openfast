import zmq
import time
from openfast_toolbox.io import FASTInputFile, FASTOutputFile
import numpy as np 

def main(rep_port: str = 'tcp://127.0.0.1:5555'): 
    context = zmq.Context()

    responder = context.socket(zmq.REP)
    responder.bind(rep_port)
    
    poller = zmq.Poller()
    poller.register(responder, zmq.POLLIN)
    
    limit = 0
    n_communications = 0 
    
    while True:
        events = dict(poller.poll(1000))

        if responder in events and events[responder] == zmq.POLLIN:
            message = responder.recv_string()
            print(f"Received request: {message}")
            
            msg_dict = {'VelH': 23, 'VelV': 0, 'BlPitchCom1': np.deg2rad(35)}
            response = ';'.join(map(str, msg_dict.values())) + ';'
            responder.send_string(response)
            
            n_communications += 1
            
        else: 
            limit += 1
            if limit == 5:
                print('Closing REP socket')
                responder.close()
                break
    
    fastin = FASTInputFile('./templatesDir/OFZMQ_test03/OFZMQ_test03.fst') 

    n_expected_communications = int(fastin['TMax'] / fastin['ZmqInDT']) + 1

    assert n_communications == n_expected_communications
    print('Test 03A ZMQ - Communication Integrity (nbr of communications): PASSED')
    
    fastout = FASTOutputFile('./templatesDir/OFZMQ_test03/OFZMQ_test03.outb').toDataFrame()
    assert np.allclose(fastout['Wind1VelX_[m/s]'][50:], [23] * len(fastout['Wind1VelX_[m/s]'][50:]))
    assert np.allclose(fastout['Wind1VelY_[m/s]'][50:], [0] * len(fastout['Wind1VelY_[m/s]'][50:]))
    assert np.allclose(fastout['BlPitchC1_[deg]'][50:], [35] * len(fastout['BlPitchC1_[deg]'][50:]))

    print('Test 03B ZMQ - Message Integrity: PASSED')
        
            
if __name__ == "__main__":
    main()