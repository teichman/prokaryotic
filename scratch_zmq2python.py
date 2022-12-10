from rich import print
from rich.console import Console
from rich.table import Table
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
            print(f"{idx=} {len(msg)=}")
            if idx == len(msg):
                break

        return result

    def _interpretField(self, msg, idx):
        # field_name, idx = self._interpretString(self, msg, idx)
        field_name_length = struct.unpack_from('i', msg, offset=idx)[0]
        idx += 4
        field_name = struct.unpack_from(f'{field_name_length}s', msg, offset=idx)[0].decode("UTF-8")
        print(f"{field_name=}")
        idx += field_name_length
        typecode = struct.unpack_from('B', msg, offset=idx)[0]
        print(f"{typecode=}")
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
        print(f"{string=}")
        idx += length
        return string, idx
    
    def _interpretStrings(self, msg, idx):
        num_strings = struct.unpack_from('i', msg, offset=idx)[0]
        idx += 4
        print(f"{num_strings=}")
        strings = []
        for _ in range(num_strings):
            string, idx = self._interpretString(msg, idx)
            strings.append(string)
        return strings, idx
    
    def _interpretArrayXd(self, msg, idx):
        rows = struct.unpack_from('i', msg, offset=idx)[0]
        idx += 4
        print(f"{rows=}")
        vec = np.zeros((rows,))
        vec[:] = struct.unpack_from('d'*rows, msg, offset=idx)
        idx += rows*8
        return vec, idx

    def _interpretArrayXXd(self, msg, idx):
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

    piot = Table(title="Protein IO", show_lines=False, highlight=True, row_styles=["", ""])
    piot.add_column("", justify="center", style="cyan", no_wrap=False)
    for mname in result["molecule_names"]:
        # col_display_name = mname.replace(' ', '\n')
        col_display_name = '\n'.join([s[:7] for s in mname.split(' ')])
        piot.add_column(col_display_name, justify="center", style="cyan", no_wrap=False, width=8)

    for idx, mname in enumerate(result["molecule_names"]):
        mat = result["protein_io_flux"]
        entries = [mname]
        for val in mat[idx, :]:
            if val > 0:
                entries.append(f'[green]{val:4.2}[/]')
            elif val < 0:
                entries.append(f'[red]{val:4.2}[/]')
            else:
                entries.append('-')
                
        # entries += [f'{val:4.2}' if abs(val) > 0 else '' for val in mat[idx, :]] 
        piot.add_row(*entries)

    console = Console()
    console.print(piot)
    
    time.sleep(0.1)
