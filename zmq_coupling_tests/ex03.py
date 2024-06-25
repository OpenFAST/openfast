import zmq
import threading
import time
from openfast_toolbox.io import FASTInputFile, FASTOutputFile
import numpy as np
import pickle 

def update_dict(update_, dict_):
    with threading.Lock():
        for key, value in update_.items():
            if key not in dict_:
                dict_[key] = []
                dict_[key].append(value)
            else:
                dict_[key].append(value)
                
    return dict_


def pub_sub_handler(context, sub_port, subscribed_messages):
    subscriber = context.socket(zmq.SUB)
    subscriber.bind(sub_port)
    subscriber.setsockopt_string(zmq.SUBSCRIBE, "")  # Subscribe to all messages
    
    while True:
        message = subscriber.recv_json()
        # print(f"Received message from subscriber: {message}")
        update_dict(message, subscribed_messages)
        
        # Assuming you want to close the subscriber socket based on specific conditions in the message
        if list(message.values())[2:] == [0.0] * (len(message) - 2):
            print('Closing SUB socket')
            subscriber.close()
            
            with open('subscribed_messages.pkl', 'wb') as f:
                pickle.dump(subscribed_messages, f)
            
            
            break

def req_rep_handler(context, rep_port, subscribed_messages):
    responder = context.socket(zmq.REP)
    responder.bind(rep_port)
    n_communications = 0
    
    poller = zmq.Poller()
    poller.register(responder, zmq.POLLIN)
    
    limit = 0
    
    while True:
        
        socks = dict(poller.poll(1000))
        
        if responder in socks and socks[responder] == zmq.POLLIN:
            message = responder.recv_string()
            # print(f"Received request: {message}")
            
            msg_dict = {'VelH': 8, 'BlPitchCom1': np.deg2rad(15)}
            response = ';'.join(map(str, msg_dict.values())) + ';'
            responder.send_string(response)
            
            n_communications += 1
        else:
            limit += 1
            if limit == 5:
                print('Closing REP socket')
                responder.close()
                break        

def main(sub_port: str = 'tcp://127.0.0.1:5557', rep_port: str = 'tcp://127.0.0.1:5555'):
    context = zmq.Context()  # Create a single context
    subscribed_messages = {}

    # Create and start threads
    pub_sub_thread = threading.Thread(target=pub_sub_handler, args=(context, sub_port, subscribed_messages))
    req_rep_thread = threading.Thread(target=req_rep_handler, args=(context, rep_port, subscribed_messages))
    
    pub_sub_thread.start()
    req_rep_thread.start()
    
    # Join threads to the main thread to keep the main program running
    pub_sub_thread.join()
    req_rep_thread.join()

    time.sleep(1)  # Wait for the threads to close
    
    fastin = FASTInputFile('./templatesDir/OFZMQ_test04/OFZMQ_test04.fst')
    fastout = FASTOutputFile('./templatesDir/OFZMQ_test04/OFZMQ_test04.outb').toDataFrame()

    assert np.allclose(fastout['Wind1VelX_[m/s]'][50:], [8] * len(fastout['Wind1VelX_[m/s]'][50:]))
    assert np.allclose(fastout['BlPitchC1_[deg]'][50:], [15] * len(fastout['BlPitchC1_[deg]'][50:]))

    print('Test 04A ZMQ - Message Integrity (REP): PASSED')
    
    assert np.allclose(fastout['Wind1VelX_[m/s]'][1:-2], subscribed_messages[' Wind1VelX_(m/s)'][:-1], atol=1e-5)
    assert np.allclose(fastout['Azimuth_[deg]'][1:-2], subscribed_messages[' Azimuth_(deg)'][:-1], atol=1e-5)
    assert np.allclose(fastout['GenTq_[kN-m]'][1:-2], subscribed_messages[' GenTq_(kN-m)'][:-1], atol=1e-5)
    assert np.allclose(fastout['GenPwr_[kW]'][1:-2], subscribed_messages[' GenPwr_(kW)'][:-1], atol=1e-5)
    
    print('Test 04B ZMQ - Message Integrity (PUB): PASSED')
    
    return True
    
if __name__ == "__main__":
    main()