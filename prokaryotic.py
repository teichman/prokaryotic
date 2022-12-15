import copy
from rich.prompt import Prompt
import readchar 
from pynput.keyboard import Key, Listener
import threading
import math
import argparse
import pickle
import sys
from rich.markdown import Markdown
from rich import print
from rich.progress import Progress
from rich.tree import Tree
from rich.console import Console
from rich.table import Table
from rich.live import Live
from rich.text import Text
from rich.layout import Layout
from rich.syntax import Syntax
from rich.panel import Panel
import numpy as np
import struct
import ipdb
import time
import zmq

g_stop = False

class MessageInterpreter:

    Double = 8
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
        elif typecode == self.Double:
            val, idx = self._interpretDouble(msg, idx)
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
    
    def _interpretDouble(self, msg, idx):
        val = struct.unpack_from('d', msg, offset=idx)[0]
        idx += 8
        return val, idx
    
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



class KeypressListener:
    
    def __init__(self, callback):
        self.callback = callback
        self.start()

    def __del__(self):
        self.thread.join()
        
    def _loopfn(self):
        while self.running:
            key = readchar.readkey()
            retval = self.callback(key)
            # We can't just call stop() because readchar.readkey() is blocking, ugh really?
            # Anyway this gets it done.
            if retval == 'stop':
                self.running = False
    
    def start(self):
        self.running = True
        self.thread = threading.Thread(target=self._loopfn)
        self.thread.start()

    def stop(self):
        self.running = False
        self.thread.join()
    
    
class Comms:
    def __init__(self, callback=None, msgfile=None):
        self.callback = callback
        self.msgdict = None
        self.msgfile = msgfile
        if self.msgfile:
            with open(self.msgfile, 'rb') as f:
                self.msgdict = pickle.load(f)
                print(self.msgdict)
        
        self.ctx = zmq.Context()
        self.sub_sock = self.ctx.socket(zmq.SUB)
        self.sub_sock.connect(f"tcp://{args.server}:53269")
        self.sub_sock.subscribe("")
        self.pub_sock = self.ctx.socket(zmq.PUB)
        self.pub_sock.connect(f"tcp://{args.server}:53270")
        self.mi = MessageInterpreter()
        logger.log("Connected.")
        self.start()

    def send_msg(self, msg):
        self.pub_sock.send(bytes(msg, 'utf-8'))

    def send_file(self, path):
        with open(path, 'r') as f:
            msg = f.readlines()
            msg = ''.join(msg)
        self.send_msg(msg)
        
    def advance(self):
        # Tell the server that we're done.
        # Include a copy of DNA programming every time, since it may have changed.
        # For now, the only action is to send dna.yaml up to the server,
        # and it interprets that as "i'm ready to advance"
        self.send_file('dna.yaml')
        
    def __del__(self):
        self.thread.join()
        
    # blocking
    def receive(self):
        msg = self.sub_sock.recv()
        # print(f"Received msg: {type(msg)} {msg}")
        result = self.mi.interpret(msg)
        # print(f"{result=}")
        return result

    def start(self):
        self.running = True
        self.thread = threading.Thread(target=self._loopfn)
        self.thread.start()
        logger.log("started thread for Comms")
        
    def stop(self):
        self.running = False
        self.thread.join()
        
    def _loopfn(self):
        while self.running:
            if self.msgdict:
                self.callback(self.msgdict)
                logger.log("Processing cached msg and exiting Comms._loopfn")
                logger.log(self.msgdict)
                return
                
            msgdict = self.receive()
            logger.log("Comms got a msg")
            if self.callback:
                logger.log("Comms is calling the callback")
                self.callback(msgdict)
    
class View:
    def __init__(self):
        self.console = Console()
        self.layout = Layout(ratio=1.0, name="layout")
        self.layout.split_column(Layout(name="header", ratio=1),
                                 Layout(name="mainrow", ratio=10),
                                 Layout(name="footer", ratio=1))
        self.layout["mainrow"].split_row(Layout(name="left_sidebar", ratio=1),
                                         Layout(name="main", ratio=4),
                                         Layout(name="log", ratio=1))
        self.layout["mainrow"]["left_sidebar"].split_column(Layout(name="cmds"),
                                                            Layout(name="metadata"))

        self.mname_map = dict()
        self.msgdict = None

        self.draw_cmds()
        self.draw_header()
        self.draw_placeholder_in_main()
        self.draw_planetary_metadata()

        #self.start()

    # def __del__(self):
    #     print("View.__del__ starting")
    #     self.thread.join()
    #     print("View.__del__ done")
        
    # def _loopfn(self):
    #     while self.running:
    #         if self.grid:
    #             self.console.clear()
    #             self.console.print(self.grid)
    #             self.grid = None
    #         time.sleep(0.1)
                
    # def start(self):
    #     self.running = True
    #     self.thread = threading.Thread(target=self._loopfn)
    #     self.thread.start()
    #     logger.log("started thread")

    # def stop(self):
    #     print(f"View.stop")
    #     self.running = False
    #     self.thread.join()

    def draw_placeholder_in_main(self):
        self.layout["mainrow"]["main"].update(Panel("Nothing to see here yet...",
                                                    border_style="dim white"))

    def draw_planetary_metadata(self, msgdict=None):
        if msgdict:
            dtstr = f"{msgdict['division_hours']:2.1f}"
        else:
            dtstr = "TBD"
            
        pm = self.layout["mainrow"]["left_sidebar"]["metadata"]
        st = Table()
        st.add_column("Species")
        st.add_column("Doubling\ntime\n(hours)")
        st.add_row("E. coli", dtstr)
        st.add_row("M. magneticum", "1.42")

        pm.update(Panel(st,
                        title="Planetary metadata",
                        border_style="dim white"))
        self.draw()

    def mname_to_idx(self, mname):
        return self.msgdict["molecule_names"].index(mname)
        
    def draw_transformation_chain(self, key):
        mname = self.mname_map[key]
        midx = self.mname_to_idx(mname)
        logger.log(f"Drawing transformation chain for {mname}, molecule idx {midx}")
        
    def draw_main_dashboard(self):
        if self.msgdict is None:
            return
        self.grid = self.generate_grid(self.msgdict)
        self.layout["mainrow"]["main"].update(Panel(self.grid,
                                                    border_style="dim white",
                                                    title='[bold]Cytosol dashboard[/]'))
        self.draw()

    def draw_dna_programming(self):
        layout = Layout(name="dna programming")
        layout.split_column(Layout(Panel(Syntax.from_path("dna.yaml"))),
                            Layout(Panel("DNA programming file: aoeu/snth.yaml"), ratio=0.1))
        self.layout["mainrow"]["main"].update(Panel(Syntax.from_path("dna.yaml"),
                                                    title="[bold]DNA Programming (dna.yaml)[/]",
                                                    border_style="dim white"))
        self.draw()
        
    def draw(self):
        self.console.clear()
        self.console.print(self.layout)

    def draw_cmds(self):
        cmds = Tree("Commands")
        cmds.add("d   :: Cytosol dashboard")
        cmds.add("p   :: View DNA Programming")
        cmds.add("RET :: Continue simulation")        
        cmds.add("[dim]i   :: Inspect pathways[/]")
        cmds.add("[dim]h   :: Molecule history[/]")
        cmds.add("[dim]b   :: Select biome[/]")
        cmds.add("q   :: Quit")
        self.layout["mainrow"]["left_sidebar"]["cmds"].update(Panel(cmds, border_style="dim white"))
        self.draw()

    def draw_header(self):
        # self.layout["header"].update(Panel(title="Prokaryotic v0.0001", style="bold white on blue"))
        # self.layout["header"].update(Text(text="\n\nProkaryotic v0.0001", style="bold white on blue", justify='center'))
        self.layout["header"].update(Text(text="", style=None, justify='center'))
        self.draw()

    def update_progress(self, msgdict):
        self.layout["footer"].update(Panel(f"Simulation progress: {msgdict['step_progress'] * 100:4.1f}%", border_style="dim white"))
        logger.log(f"{msgdict['step_progress']=}")
        self.draw()
        
    def new_step_data(self, msgdict):
        self.layout["footer"].update(Panel("", border_style="dim white"))
        self.draw_planetary_metadata(msgdict)
        self.msgdict = msgdict
        self.draw_main_dashboard()
        
    # def display(self, msgdict):
    #     self.msgdict = msgdict
    #     self.grid = self.generate_grid(self.msgdict)
    #     logger.log(f"display sees self.grid: {self.grid} {id(self.grid)}")
    #     self.draw()
                
    def generate_grid(self, msgdict):
        logger.log(f"In generate_grid")
        logger.log(f"{msgdict['cytosol_contents_denatured_hist_avg']=}")
        # Cytosol contents table
        cct = Table(title="Molecule Stats (average)",
                    title_style="black on rgb(255,255,255)", show_lines=False, highlight=True, row_styles=["", ""])
        cct.add_column("Molecule", justify="right", style="cyan", min_width=10)
        cct.add_column("Key".replace(' ', '\n'), style='cyan')
        cct.add_column("Cytosol Contents (#)".replace(' ', '\n'), justify="center")
        cct.add_column("Cytosol Concentrations (mM)".replace(' ', '\n'), justify="center")
        cct.add_column("Proteasome\nRecycling of\nDenatured\n(#/s)", justify="center")
        cct.add_column("Denatured (#)".replace(' ', '\n'), justify="center")
        self.mname_map.clear()
        for idx, mname in enumerate(msgdict["molecule_names"]):
            key = chr(ord('A') + idx)
            self.mname_map[key] = mname
            cct.add_row(mname,
                        f"{key}",
                        val2str(msgdict["cytosol_contents_hist_avg"][idx]),
                        val2str(msgdict["cytosol_concentration_hist_avg"][idx]),
                        val2str(msgdict["proteasome_action"][idx]),
                        val2str(msgdict["cytosol_contents_denatured_hist_avg"][idx]))
            

        # Protein IO table
        piot = Table(title="Protein IO (#/s)",
                     title_style="black on rgb(255,255,255)", show_lines=False, highlight=True, row_styles=["", ""])
        piot.add_column("", justify="right", style="cyan", no_wrap=True)  # min_width=10, width=10
        col_indices_to_keep = []
        for idx, mname in enumerate(msgdict["molecule_names"]):
            if msgdict["protein_io_flux"][:, idx].sum() == 0:
                continue
            col_indices_to_keep.append(idx)
            max_num_chars = 100
            col_display_name = '\n'.join([s[:max_num_chars] for s in mname.split(' ')])
            piot.add_column(col_display_name, justify="center", style="cyan", width=8)
        cols_stripped = msgdict["protein_io_flux"][:, col_indices_to_keep]

        for idx, mname in enumerate(msgdict["molecule_names"]):
            if msgdict["protein_io_flux"][idx, :].sum() == 0:
                continue
            entries = [mname]
            entries += [val2str_colored(val) for val in cols_stripped[idx, :]]
            piot.add_row(*entries)

        num_grid_cols = 2
        grid = Table(style="on black", show_header=False, show_edge=False)
        grid.add_row(cct)
        grid.add_row(piot)
        
        return grid

    
class Controller:
    def __init__(self, msgfile=None):
        self.view = View()
        self.comms = Comms(self.handle_msgdict, msgfile)
        self.keypress_listener = KeypressListener(self.handle_keypress)
        

    def __del__(self):
        pass

    def handle_keypress(self, key):
        # C^a comes through with ord(key) == 1.
        logger.log(f"Controller.handle_keypress got {key}")
        try:
            logger.log(f"Controller.handle_keypress got {ord(key)=}")
        except:
            pass
        
        if key == 'q':
            self.running = False
            return 'stop'
        elif key == 'p':
            self.view.draw_dna_programming()
        elif key == 'd':
            self.view.draw_main_dashboard()
        elif len(key) == 1 and ord(key) >= ord('A') and ord(key) <= ord('Z'):
            self.view.draw_transformation_chain(key)
        elif ord(key) == 10:  # RET
            self.comms.advance()
        # elif key == "":
        #     self.view.draw_main_dashboard()
        # elif key == 'm':
        #     usermsg = Prompt.ask("Message to send")
        #     print(f"Sending {usermsg}")
        #     self.view.draw()
        # elif key == 's':
        #     self.view.draw_molecule_stats()
        
    def run(self):
        self.running = True
        
        while self.running:
            time.sleep(0.1)

        #self.view.stop()
        self.comms.stop()

    def handle_msgdict(self, msgdict):
        logger.log(f"In Controller.handle_msgdict with {msgdict}")
        if 'protein_io_flux' in msgdict:
            self.view.new_step_data(msgdict)
        elif 'step_progress' in msgdict:
            self.view.update_progress(msgdict)
        
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
        self.f.write(f"{msg}\n")

    def __del__(self):
        self.f.close()

if __name__ == "__main__":
    logger = Logger(".prokaryotic.log")
    
    parser = argparse.ArgumentParser()
    parser.add_argument("--savemsg", type=str)
    parser.add_argument("-m", "--msg", type=str)
    parser.add_argument('-s', '--server', type=str, default='127.0.0.1')
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

    # elif args.msg:
    #     with open(args.msg, 'rb') as f:
    #         msgdict = pickle.load(f)
    #     view = View(start=True)
    #     while True:
    #         logger.log("Calling display()")
    #         view.display(msgdict)
    #         time.sleep(0.4)

    else:
        controller = Controller(args.msg)
        controller.run()
    
