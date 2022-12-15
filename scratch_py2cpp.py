import argparse
import time
import zmq

class Comms:
    def __init__(self):
        self.ctx = zmq.Context()
        self.pub_sock = self.ctx.socket(zmq.PUB)
        self.pub_sock.connect(f"tcp://{args.server}:53270")

    def send_msg(self, msg):
        self.pub_sock.send(bytes(msg, 'utf-8'))
        print("Sent msg")

    def send_file(self, path):
        with open(path, 'r') as f:
            msg = f.readlines()
            msg = ''.join(msg)
        self.send_msg(msg)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-n", '--name', type=str, default='aoeu')
    parser.add_argument('-s', '--server', type=str, default='127.0.0.1')
    args = parser.parse_args()
    
    comms = Comms()
    while True:
        time.sleep(1)
        # comms.send_msg(args.name)
        comms.send_file('dna.yaml')
        
    
