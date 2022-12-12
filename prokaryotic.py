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
            #print(f"{idx=} {len(msg)=}")
            if idx == len(msg):
                break

        return result

    def _interpretField(self, msg, idx):
        # field_name, idx = self._interpretString(self, msg, idx)
        field_name_length = struct.unpack_from('i', msg, offset=idx)[0]
        idx += 4
        field_name = struct.unpack_from(f'{field_name_length}s', msg, offset=idx)[0].decode("UTF-8")
        #print(f"{field_name=}")
        idx += field_name_length
        typecode = struct.unpack_from('B', msg, offset=idx)[0]
        #print(f"{typecode=}")
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
        #print(f"{string=}")
        idx += length
        return string, idx
    
    def _interpretStrings(self, msg, idx):
        num_strings = struct.unpack_from('i', msg, offset=idx)[0]
        idx += 4
        #print(f"{num_strings=}")
        strings = []
        for _ in range(num_strings):
            string, idx = self._interpretString(msg, idx)
            strings.append(string)
        return strings, idx
    
    def _interpretArrayXd(self, msg, idx):
        rows = struct.unpack_from('i', msg, offset=idx)[0]
        idx += 4
        #print(f"{rows=}")
        vec = np.zeros((rows,))
        vec[:] = struct.unpack_from('d'*rows, msg, offset=idx)
        idx += rows*8
        return vec, idx

    def _interpretArrayXXd(self, msg, idx):
        rows, cols = struct.unpack_from('ii', msg, offset=idx)
        idx += 8
        #print(f"{rows=} {cols=}")
        mat = np.zeros((rows, cols))
        for c in range(cols):
            mat[:, c] = struct.unpack_from('d'*rows, msg, offset=idx)
            idx += rows*8
        return mat, idx


class Comms:
    def __init__(self):
        self.ctx = zmq.Context()
        self.socket = self.ctx.socket(zmq.SUB)
        #socket.bind("tcp://*:5555")
        self.socket.connect("tcp://127.0.0.1:53269")
        self.socket.subscribe("")
        self.mi = MessageInterpreter()
        print("Connected.")

    # blocking
    def receive(self):
        msg = self.socket.recv()
        # print(f"Received msg: {type(msg)} {msg}")
        result = self.mi.interpret(msg)
        # print(f"{result=}")
        return result
    
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
                print("in _loopfn")
                print(f"_loopfn sees self.grid: {self.grid} {id(self.grid)}")
                if self.grid:
                    print("got new grid")
                    live.update(self.grid, refresh=True)
                    self.grid = None
                time.sleep(0.4)

    def start(self):
        self.thread = threading.Thread(target=self._loopfn)
        self.thread.start()
        print("started thread")
        
    def display(self, msgdict):
        self.grid = self.generate_grid(msgdict)
        print(f"display sees self.grid: {self.grid} {id(self.grid)}")
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
        piot.add_column("."*10, justify="right", style="cyan on black", min_width=10, width=10)
        for mname in msgdict["molecule_names"]:
            # col_display_name = mname.replace(' ', '\n')
            col_display_name = '\n'.join([s[:7] for s in mname.split(' ')])
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

        menu = ""
        menu += "`p` :: Protein IO table"
        menu += "     "
        menu += "`a` :: Actions"
        menu = Markdown(menu)

        grid.add_row(menu)
        grid.add_row("Prokaryotic v0.1", style="white on blue")
        
        return grid

    
class Controller:
    def __init__(self, view):
        self.comms = Comms()
        self.view = view

    def run(self):
        while True:
            msgdict = self.comms.receive()
            print("got msgdict")
            self.view.display(msgdict)
        
def val2str(val):
    if val == 0:
        return '-'
    if abs(val) > 1000 or (abs(val) < 1e-3 and abs(val) > 0):
        return f'{val:02.1e}'
    else:
        return f'{val:4.2f}'

def val2str_colored(val):
    if val > 0:
        # return f'[green]{val:4.2}[/]'
        return f'[green]{val:02.1e}[/]'
    elif val < 0:
        #return f'[red]{val:4.2}[/]'
        return f'[red]{val:02.1e}[/]'
    else:
        return '-'
        
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--savemsg", type=str)
    parser.add_argument("-m", "--msg", type=str)
    args = parser.parse_args()

    # Getting keypresses requires perms.
    # def on_press(key):
    #     print('{0} pressed'.format(
    #         key))
    # def on_release(key):
    #     print('{0} release'.format(
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
            print(f"Saved message to {args.savemsg}")
        sys.exit(0)

    elif args.msg:
        with open(args.msg, 'rb') as f:
            msgdict = pickle.load(f)
        view = View(start=True)
        while True:
            print("Calling display()")
            view.display(msgdict)
            time.sleep(0.4)

    else:
         view = View(start=True)
         controller = Controller(view)
         controller.run()

    
