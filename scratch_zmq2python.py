from rich import print
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
    String = 9
    ArrayXd = 10
    ArrayXXd = 11
    
    def __init__(self):
        pass

    def interpret(self, msg):
        magic = struct.unpack('B', msg[:1])[0]
        assert magic == 13

        result = {}
        idx = 1
        while True:
            name, val, idx = self.interpretField(msg, idx)
            result[name] = val
            print(f"{idx=} {len(msg)=}")
            if idx == len(msg):
                break

        return result
            
    def interpretField(self, msg, idx):
        field_name_length = struct.unpack_from('i', msg, offset=idx)[0]
        idx += 4
        field_name = struct.unpack_from(f'{field_name_length}s', msg, offset=idx)[0].decode("UTF-8")
        print(f"{field_name=}")
        idx += field_name_length
        typecode = struct.unpack_from('B', msg, offset=idx)[0]
        print(f"{typecode=}")
        idx += 1
        if typecode == self.ArrayXXd:
            val, idx = self.interpretArrayXXd(msg, idx)
            return field_name, val, idx
        if typecode == self.ArrayXd:
            val, idx = self.interpretArrayXd(msg, idx)
            return field_name, val, idx
        assert False

    def interpretArrayXd(self, msg, idx):
        rows = struct.unpack_from('i', msg, offset=idx)[0]
        idx += 4
        print(f"{rows=}")
        vec = np.zeros((rows,))
        vec[:] = struct.unpack_from('d'*rows, msg, offset=idx)
        idx += rows*8
        return vec, idx

    def interpretArrayXXd(self, msg, idx):
        rows, cols = struct.unpack_from('ii', msg, offset=idx)
        idx += 8
        print(f"{rows=} {cols=}")
        mat = np.zeros((rows, cols))
        for c in range(cols):
            mat[:, c] = struct.unpack_from('d'*rows, msg, offset=idx)
            idx += rows*8
        return mat, idx

mi = MessageInterpreter()

while True:
    msg = socket.recv()
    print(f"Received msg: {type(msg)} {msg}")
    result = mi.interpret(msg)
    print(f"{result=}")
    print(result['matrix_something'])
    time.sleep(0.1)
