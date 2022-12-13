import time
import threading
import readchar

class KeypressListener:
    
    def __init__(self, callback):
        self.callback = callback
        self.start()

    def __del__(self):
        print("KeypressListener.__del__ starting")
        self.thread.join()
        print("KeypressListener.__del__ done")
        
    def _loopfn(self):
        while self.running:
            key = readchar.readkey()
            retval = self.callback(key)
            if retval == 'stop':
                self.running = False
    
    def start(self):
        self.running = True
        self.thread = threading.Thread(target=self._loopfn)
        self.thread.start()

    def stop(self):
        print(f"KeypressListener.stop")
        self.running = False
        self.thread.join()


class Controller:
    def __init__(self):
        self.stop = False
        self.listener = KeypressListener(self.handle_keypress)        

    def __del__(self):
        print(f"{self.__class__.__name__}.__del__")

    def handle_keypress(self, key):
        print(f"Controller got key: {key}")
        if key == 'q':
            self.stop = True
            return 'stop'
        
    def run(self):
        while not self.stop:
            print(f"{self.stop=}")
            time.sleep(0.1)
        print(f"Controller.run is ending.")

def handle_keypress_standalone(key):
    print(f"Standalone got key: {key}")
        
        
if __name__ == "__main__":
    con = Controller()
    con.run()
    
