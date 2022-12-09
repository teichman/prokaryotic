import numpy as np
import struct
import ipdb
import time
import zmq

context = zmq.Context()
socket = context.socket(zmq.SUB)
#socket.bind("tcp://*:5555")
socket.connect("tcp://127.0.0.1:53269")
socket.subscribe("")
print("Connected.")


class MessageInterpreter:

    ArrayXd = 10
    ArrayXXd = 11
    
    def __init__(self):
        pass

    def interpret(self, msg):
        typecode = struct.unpack('B', msg[:1])[0]
        if typecode == self.ArrayXXd:
            return self.interpretArrayXXd(msg[1:])
        assert False

    def interpretArrayXXd(self, msg):
        rows, cols = struct.unpack('ii', msg[:8])
        msg = msg[8:]
        print(f"{rows=} {cols=}")
        mat = np.zeros((rows, cols))
        for c in range(cols):
            mat[:, c] = struct.unpack('d'*rows, msg[:8*rows])
            msg = msg[8*rows:]
        return mat
    


mi = MessageInterpreter()

while True:
    #  Wait for next request from client
    time.sleep(1)
    print("Waiting for a msg.")
    msg = socket.recv()
    print(f"Received msg: {type(msg)} {msg}")
    print(mi.interpret(msg))

    
    # # https://docs.python.org/3/library/struct.html#format-characters
    # typecode = struct.unpack('B', msg[:1])[0]
    # print(f"{typecode}")
    # msg = msg[1:]

    # rows = struct.unpack('i', msg[:4])[0]
    # print(f"{rows=}")
    # msg = msg[4:]

    # values = struct.unpack('d'*rows, msg[:8*rows])
    # print(f"{values=}")

    
    
    # breakpoint()

    #  Do some 'work'


    # #  Send reply back to client
    # socket.send_string("World")
