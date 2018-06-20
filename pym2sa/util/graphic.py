import logging
import random
from typing import List, Tuple, TypeVar

from bokeh.application import Application
from bokeh.application.handlers import FunctionHandler
from bokeh.client import ClientSession
from bokeh.io import curdoc
from bokeh.layouts import column
from bokeh.models import ColumnDataSource, HoverTool
from bokeh.server.server import Server
from bokeh.plotting import Figure, show
from tornado import gen

logger = logging.getLogger(__name__)
import threading
S = TypeVar('S')

# https://groups.google.com/a/continuum.io/forum/#!topic/bokeh/zVzLyvJFWtE


class ScatterPlot(threading.Thread):

	def __init__(self, title: str, circle_size: int=3, tools: []=None):
		super().__init__()

		if tools is None:
			tools = [HoverTool(tooltips=[("index", "$index"), ("(x,y)", "($x, $y)")]), 'save', 'lasso_select']

		self.figure_title = title
		self.plot_tools = tools
		self.client = ClientSession()
		self.plot_circle_size = circle_size
		self.source = ColumnDataSource(data=dict(x=[], y=[]))
		self.doc = curdoc()

	def plot(self, solution_list: List[S], shw=True) -> None:
		print("Opening Bokeh application on http://localhost:{0}/".format(5001))
		x_values, y_values, z_values = self.__get_objectives(solution_list)
		self.source.stream({'x': x_values, 'y': y_values})

		# Set-up figure
		self.figure = Figure(output_backend="webgl", sizing_mode='scale_width',
		                     title=self.figure_title, tools=self.plot_tools)
		self.figure.circle(x='x', y='y', name=self.figure_title, size=self.plot_circle_size, source=self.source)

		# Add to currdoc
		self.doc.add_root(column(self.figure))
		self.client.push(self.doc)

		if shw:
			self.client.show()

	def update(self, solution_list: List[S]):
		x_values, y_values, z_values = self.__get_objectives(solution_list)
		self.source.stream({'x': x_values, 'y': y_values}, rollover=len(solution_list))

	def __get_objectives(self, solution_list: List[S]) -> Tuple[list, list, list]:
		""" Get coords (x,y,z) from a solution_list.

		:return: A tuple with (x,y,z) values. The third might be empty if working with a problem with 2 objectives."""
		if solution_list is None:
			raise Exception("Solution list is none!")

		points = list(solution.objectives for solution in solution_list)

		x_values, y_values = [point[0] for point in points], [point[1] for point in points]
		z_values = []

		return x_values, y_values, z_values
