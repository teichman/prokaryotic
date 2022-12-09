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
    print("Waiting for a message.")
    message = socket.recv()
    print(f"Received request: {message}")

    #  Do some 'work'
    time.sleep(1)

    # #  Send reply back to client
    # socket.send_string("World")
