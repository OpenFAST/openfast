import zmq 
from openfast_toolbox.io import FASTOutputFile
import numpy as np 

def update_dict(update_, dict_):
    for key, value in update_.items():
        if key not in dict_:
            dict_[key] = []  
            dict_[key].append(value)
        else: 
            dict_[key].append(value)
            
    return dict_


def main(sub_port: str = 'tcp://127.0.0.1:5557'): 
    context = zmq.Context()
    
    subscriber = context.socket(zmq.SUB)
    subscriber.bind(sub_port)
    subscriber.setsockopt_string(zmq.SUBSCRIBE, "") # for now subscribe to all

    poller = zmq.Poller()
    poller.register(subscriber, zmq.POLLIN)
    
    subscribed_messages = {}
        
    while True:
        events = dict(poller.poll(1000))
        
        if subscriber in events and events[subscriber] == zmq.POLLIN:
            message = subscriber.recv_json()
            print(f"Received message: {message}")
            update_dict(message, subscribed_messages)
            
            
            # Assuming you want to close the subscriber socket based on specific conditions in the message
            if list(message.values())[2:] == [0.0] * (len(message) - 2):
                print('Closing SUB socket')
                subscriber.close()
                break
    
    
    fastout = FASTOutputFile('./templatesDir/OFZMQ_test01/OFZMQ_test01_2.outb').toDataFrame()
    
    assert np.allclose(fastout['Wind1VelX_[m/s]'][:-2], subscribed_messages[' Wind1VelX_(m/s)'][:-1], atol=1e-5)
    assert np.allclose(fastout['Azimuth_[deg]'][:-2], subscribed_messages[' Azimuth_(deg)'][:-1], atol=1e-5)
    assert np.allclose(fastout['GenTq_[kN-m]'][:-2], subscribed_messages[' GenTq_(kN-m)'][:-1], atol=1e-5)
    assert np.allclose(fastout['GenPwr_[kW]'][:-2], subscribed_messages[' GenPwr_(kW)'][:-1], atol=1e-5)
    
    print('Test 01 ZMQ - Message Integrity: PASSED')
    
    return True   


if __name__ == "__main__":
    main()
    