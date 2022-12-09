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


while True:
    #  Wait for next request from client
    print("Waiting for a msg.")
    msg = socket.recv()
    print(f"Received msg: {type(msg)} {msg}")
    # https://docs.python.org/3/library/struct.html#format-characters
    typecode = struct.unpack('B', msg[:1])[0]
    print(f"{typecode}")
    msg = msg[1:]

    rows = struct.unpack('i', msg[:4])[0]
    print(f"{rows=}")
    msg = msg[4:]

    values = struct.unpack('d'*rows, msg[:8*rows])
    print(f"{values=}")
    
    # breakpoint()

    #  Do some 'work'
    time.sleep(1)

    # #  Send reply back to client
    # socket.send_string("World")
