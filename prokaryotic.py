from pynput.keyboard import Key, Listener
import threading
import math
import argparse
import pickle
import sys
from rich.markdown import Markdown
from rich import print
from rich.console import Console
from rich.table import Table
from rich.live import Live
import numpy as np
import struct
import ipdb
import time
import zmq

class MessageInterpreter:

    String = 10
    Strings = 11
    ArrayXd = 12
    ArrayXXd = 13
    
    def __init__(self):
        pass

    def interpret(self, msg):
        magic = struct.unpack('B', msg[:1])[0]
        assert magic == 13

        result = {}
        idx = 1
        while True:
            name, val, idx = self._interpretField(msg, idx)
            result[name] = val
            logger.log(f"{idx=} {len(msg)=}")
            if idx == len(msg):
                break

        return result

    def _interpretField(self, msg, idx):
        # field_name, idx = self._interpretString(self, msg, idx)
        field_name_length = struct.unpack_from('i', msg, offset=idx)[0]
        idx += 4
        field_name = struct.unpack_from(f'{field_name_length}s', msg, offset=idx)[0].decode("UTF-8")
        logger.log(f"{field_name=}")
        idx += field_name_length
        typecode = struct.unpack_from('B', msg, offset=idx)[0]
        logger.log(f"{typecode=}")
        idx += 1
        if typecode == self.ArrayXXd:
            val, idx = self._interpretArrayXXd(msg, idx)
        elif typecode == self.ArrayXd:
            val, idx = self._interpretArrayXd(msg, idx)
        elif typecode == self.Strings:
            val, idx = self._interpretStrings(msg, idx)
        else:
            assert False
            
        return field_name, val, idx

    def _interpretString(self, msg, idx):
        length = struct.unpack_from('i', msg, offset=idx)[0]
        idx += 4
        string = struct.unpack_from(f'{length}s', msg, offset=idx)[0].decode("UTF-8")
        logger.log(f"{string=}")
        idx += length
        return string, idx
    
    def _interpretStrings(self, msg, idx):
        num_strings = struct.unpack_from('i', msg, offset=idx)[0]
        idx += 4
        logger.log(f"{num_strings=}")
        strings = []
        for _ in range(num_strings):
            string, idx = self._interpretString(msg, idx)
            strings.append(string)
        return strings, idx
    
    def _interpretArrayXd(self, msg, idx):
        rows = struct.unpack_from('i', msg, offset=idx)[0]
        idx += 4
        logger.log(f"{rows=}")
        vec = np.zeros((rows,))
        vec[:] = struct.unpack_from('d'*rows, msg, offset=idx)
        idx += rows*8
        return vec, idx

    def _interpretArrayXXd(self, msg, idx):
        rows, cols = struct.unpack_from('ii', msg, offset=idx)
        idx += 8
        logger.log(f"{rows=} {cols=}")
        mat = np.zeros((rows, cols))
        for c in range(cols):
            mat[:, c] = struct.unpack_from('d'*rows, msg, offset=idx)
            idx += rows*8
        return mat, idx


class Comms:
    def __init__(self, callback=None):
        self.callback = callback
        self.ctx = zmq.Context()
        self.socket = self.ctx.socket(zmq.SUB)
        #socket.bind("tcp://*:5555")
        self.socket.connect("tcp://127.0.0.1:53269")
        self.socket.subscribe("")
        self.mi = MessageInterpreter()
        logger.log("Connected.")

        
    # blocking
    def receive(self):
        msg = self.socket.recv()
        # print(f"Received msg: {type(msg)} {msg}")
        result = self.mi.interpret(msg)
        # print(f"{result=}")
        return result

    def start(self):
        self.thread = threading.Thread(target=self._loopfn)
        self.thread.start()
        logger.log("started thread for Comms")
    
    def _loopfn(self):
        while True:
            msgdict = self.receive()
            logger.log("Comms got a msg")
            if self.callback:
                logger.log("Comms is calling the callback")
                self.callback(msgdict)
    
class View:
    def __init__(self, start=False):
        self.console = Console()
        self.grid = None
        self.thread = None
        if start:
            self.start()

    def _loopfn(self):        
        with Live("loading...", refresh_per_second=4, screen=True) as live:
            while True:
                if self.grid:
                    live.update(self.grid, refresh=True)
                    self.grid = None
                time.sleep(0.4)

    def start(self):
        self.thread = threading.Thread(target=self._loopfn)
        self.thread.start()
        logger.log("started thread")
        
    def display(self, msgdict):
        self.grid = self.generate_grid(msgdict)
        logger.log(f"display sees self.grid: {self.grid} {id(self.grid)}")
        time.sleep(0.4)
                
    def generate_grid(self, msgdict):
        # Cytosol contents table
        cct = Table(title="Molecule stats", title_style="black on rgb(255,255,255)", show_lines=False, highlight=True, row_styles=["", ""])
        cct.add_column("Molecule", justify="right", style="cyan", min_width=10)
        cct.add_column("Cytosol Contents (#)".replace(' ', '\n'), justify="center", style="magenta")
        cct.add_column("Cytosol Concentrations (mM)".replace(' ', '\n'), justify="center", style="magenta")
        for idx, mname in enumerate(msgdict["molecule_names"]):
            cct.add_row(mname,
                        val2str(msgdict["cytosol_contents_hist_avg"][idx]),
                        val2str(msgdict["cytosol_concentration_hist_avg"][idx]))

        # Protein IO table
        piot = Table(title="Protein IO", title_style="black on rgb(255,255,255)", show_lines=False, highlight=True, row_styles=["", ""])
        piot.add_column("", justify="right", style="cyan", no_wrap=True)  # min_width=10, width=10
        for mname in msgdict["molecule_names"]:
            # col_display_name = mname.replace(' ', '\n')
            max_num_chars = 100
            col_display_name = '\n'.join([s[:max_num_chars] for s in mname.split(' ')])
            piot.add_column(col_display_name, justify="center", style="cyan", width=8)

        for idx, mname in enumerate(msgdict["molecule_names"]):
            mat = msgdict["protein_io_flux"]
            entries = [mname]
            entries += [val2str_colored(val) for val in mat[idx, :]]
            # for val in mat[idx, :]:
            #     if val > 0:
            #         entries.append(f'[green]{val:4.2}[/]')
            #     elif val < 0:
            #         entries.append(f'[red]{val:4.2}[/]')
            #     else:
            #         entries.append('-')
            # entries += [f'{val:4.2}' if abs(val) > 0 else '' for val in mat[idx, :]] 
            piot.add_row(*entries)

        num_grid_cols = 2
        # grid = Table.grid(*([""] * num_grid_cols))
        # title="Prokaryotic v0.1", title_style="white on blue"
        grid = Table(style="on black", show_header=False, show_edge=False)
        # grid.add_column("")
        # grid.add_column("")
        grid.add_row(cct)
        grid.add_row(piot)

        grid.add_row("Prokaryotic v0.1", style="white on blue")
        
        return grid

    
class Controller:
    def __init__(self, view):
        self.view = view
        self.comms = Comms(self.handle_msgdict)
        self.comms.start()
    
    def run(self):
        while True:
            time.sleep(0.1)

    def handle_msgdict(self, msgdict):
        logger.log("handle_msgdict")
        self.view.display(msgdict)
        
def val2str(val):
    if val == 0:
        return '-'
    elif abs(val) > 1000 or (abs(val) < 1e-2):
        return f'{val:02.1e}'
    else:
        return f'{val:4.2f}'

def val2str_colored(val):
    color = 'white'
    if val > 0:
        color = 'green'
    elif val < 0:
        color = 'red'
    return f"[{color}]{val2str(val)}[/]"

class Logger:
    def __init__(self, path):
        self.path = path
        self.f = open(path, 'w')

    def log(self, msg):
        self.f.write(msg + '\n')

    def __del__(self):
        self.f.close()

if __name__ == "__main__":
    logger = Logger(".prokaryotic.log")
    
    parser = argparse.ArgumentParser()
    parser.add_argument("--savemsg", type=str)
    parser.add_argument("-m", "--msg", type=str)
    args = parser.parse_args()

    # Getting keypresses requires perms.
    # def on_press(key):
    #     logger.log('{0} pressed'.format(
    #         key))
    # def on_release(key):
    #     logger.log('{0} release'.format(
    #         key))
    #     if key == Key.esc:
    #         # Stop listener
    #         return False    
    # with Listener(on_press=on_press, on_release=on_release) as listener:
    #     listener.join()

    if args.savemsg:
        comms = Comms()
        msgdict = comms.receive()
        with open(args.savemsg, 'wb') as f:
            pickle.dump(msgdict, f)
            logger.log(f"Saved message to {args.savemsg}")
        sys.exit(0)

    elif args.msg:
        with open(args.msg, 'rb') as f:
            msgdict = pickle.load(f)
        view = View(start=True)
        while True:
            logger.log("Calling display()")
            view.display(msgdict)
            time.sleep(0.4)

    else:
         view = View(start=True)
         controller = Controller(view)
         controller.run()

    
