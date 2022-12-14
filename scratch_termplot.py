import termplotlib as tpl
import numpy as np
from rich.layout import Layout
from rich.console import Console
from rich.panel import Panel

x = np.linspace(0, 2 * np.pi, 24 * 60 * 60)
y = np.sin(x) * 1e8

fig = tpl.figure()
fig.plot(x / (2*np.pi), y, width=100, height=30)

console = Console()
layout = Layout()
layout.split_column(Layout(name='header', ratio=1), Layout(name='main', ratio=10))
layout['main'].split_row(Layout(name='left'), Layout(name='right'))
layout['main']['left'].split_column(Layout(name='top'), Layout(name='bottom'))
layout['main']['right'].split_column(Layout(name='top'), Layout(name='bottom'))
layout['main']['right']['top'].update(Panel(fig.get_string(), title="Molecule history during cell lifecycle"))
console.print(layout)
# fig.show()
