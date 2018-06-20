import time
import random

from bokeh.client import push_session, ClientSession
from bokeh.io import curdoc
from bokeh.models import ColumnDataSource
from bokeh.plotting import figure, show

# prepare output to server

source = ColumnDataSource(data=dict(x=[], y=[]))
doc = curdoc()
p1 = figure(plot_width=400, plot_height=400)
p1.line(x='x', y='y', source=source)
doc.add_root(p1)

c = ClientSession()
c.push(doc)
c.show()

for i in range(0, 100000):
	xx = [random.randrange(1, 10) for _ in range(0, 4)]
	yy = [random.randrange(1, 10) for _ in range(0, 4)]
	source.stream({'x': xx, 'y': yy})