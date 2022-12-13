import argparse
import time
import zmq

class Comms:
    def __init__(self):
        self.ctx = zmq.Context()
        self.socket = self.ctx.socket(zmq.PUB)
        self.socket.connect("tcp://127.0.0.1:53270")

    def send_msg(self, msg):
        self.socket.send(bytes(msg, 'utf-8'))
        print("Sent msg")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-n", '--name', type=str, default='aoeu')
    args = parser.parse_args()
    
    comms = Comms()
    while True:
        time.sleep(1)
        comms.send_msg(args.name)
        
    
