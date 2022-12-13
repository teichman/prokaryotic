import pickle
from time import monotonic

from rich.table import Table
from textual.app import App, ComposeResult
from textual.containers import Container, Grid
from textual.widgets import Button, Header, Footer, Static, Input
from textual.widget import Widget
from textual.reactive import reactive

from prokaryotic import View, Logger, val2str_colored

class TimeDisplay(Static):
    """A widget to display elapsed time."""
    start_time = reactive(monotonic)
    time = reactive(0.0)

    def on_mount(self) -> None:
        """Event handler called when widget is added to the app."""
        self.set_interval(1 / 60, self.update_time)

    def update_time(self) -> None:
        """Method to update the time to the current time."""
        self.time = monotonic() - self.start_time

    def watch_time(self, time: float) -> None:
        """Called when the time attribute changes."""
        minutes, seconds = divmod(time, 60)
        hours, minutes = divmod(minutes, 60)
        self.update(f"{hours:02,.0f}:{minutes:02.0f}:{seconds:05.2f}")

class Stopwatch(Static):
    """A stopwatch widget."""

    def on_button_pressed(self, event: Button.Pressed) -> None:
        """Event handler called when a button is pressed."""
        if event.button.id == "start":
            self.add_class("started")
        elif event.button.id == "stop":
            self.remove_class("started")

    def compose(self) -> ComposeResult:
        """Create child widgets of a stopwatch."""
        yield Button("Start", id="start", variant="success")
        yield Button("Stop", id="stop", variant="error")
        yield Button("Reset", id="reset")
        yield TimeDisplay("00:00:00.00")


class StopwatchApp(App):
    """A Textual app to manage stopwatches."""

    CSS_PATH = "stopwatch03.css"
    BINDINGS = [("d", "toggle_dark", "Toggle dark mode"),
                ("p", "switch_to_protein_io", "Protein IO")]

    def compose(self) -> ComposeResult:
        """Create child widgets for the app."""
        # yield Header()
        yield Footer()
        yield Container(Stopwatch(), Stopwatch(), Stopwatch())

    def action_toggle_dark(self) -> None:
        """An action to toggle dark mode."""
        self.dark = not self.dark

        

class Hello(Widget):
    
    def render(self):
        with open("example_msg.pkl", 'rb') as f:
            msgdict = pickle.load(f)
        
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
        
        return piot
    

        

class ProkaryoticApp(App):
    # CSS_PATH = "prokaryotic.css"
    BINDINGS = [('q', 'quit', 'Quit'),
                ("p", "switch_to_protein_io", "Protein IO")]

    def compose(self):
        with open("example_msg.pkl", 'rb') as f:
            msgdict = pickle.load(f)
        view = View()
        self.grid = view.generate_grid(msgdict)
            
        yield Header()
        yield Static(self.grid, markup=False)
        yield Hello()
        # yield Input(placeholder="Enter something")
        yield Footer()

    def on_key(self, event):
        logger.log(f"Got key: {event}")

    def action_switch_to_protein_io(self):
        # self.screen.styles.background = "white"
        self.query_one(Static).update(renderable=self.grid)

    # def on_input_submitted(self, event):
    #     logger.log(f"Got input: {event}")
    #     self.query_one(Input).value = ""


class KeypressWrapper(App):
    def on_key(self, event):
        print(f"Got key: {event}")
        logger.log(f"Got key: {event}")
    
    def render(self):
        return None

        
if __name__ == "__main__":
    #app = StopwatchApp()

    logger = Logger(".prokaryotic.log")    
    app = KeypressWrapper()
    app.run()

    # logger = Logger(".prokaryotic.log")    
    # app = ProkaryoticApp()
    # app.run()

