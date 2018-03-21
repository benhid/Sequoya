import logging

from bokeh.application import Application
from bokeh.application.handlers import FunctionHandler
from bokeh.layouts import column
from bokeh.models import ColumnDataSource, HoverTool
from bokeh.server.server import Server
from bokeh.plotting import Figure, show

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class ScatterPlot:

    def __init__(self, title: str, circle_size: int = 10, period_milliseconds: int = 50, tools: [] = None):
        if tools is None:
            tools = [HoverTool(tooltips=[("index", "$index"), ("(x,y)", "($x, $y)")]), 'save']

        self.figure_title = title
        self.app_title = 'pym2sa'
        self.plot_tools = tools
        self.plot_circle_size = circle_size
        self.period_milliseconds = period_milliseconds

        self.source = ColumnDataSource({'x': [], 'y': [], 'seq_names': []})
        self.points = []

    def simple_plot(self, x: list, y: list, name: str) -> None:
        fig = Figure(output_backend="webgl", sizing_mode='scale_width',
                     title=self.figure_title, tools=self.plot_tools)
        fig.circle(x=x, y=y, name=name, size=self.plot_circle_size)
        show(fig)

    def interactive_plot(self, port: int = 5001) -> None:
        logger.info("Opening Bokeh application on http://localhost:{0}/".format(port))

        app = {'/' + self.app_title: Application(FunctionHandler(self.__make_interactive_document))}

        server = Server(app, port=port)
        server.show('/' + self.app_title)
        server.run_until_shutdown()

    def replace_points(self, points: list) -> None:
        print("updating points")
        self.points = points

    def __update_data_with_callback(self) -> None:
        if self.points is not []:
            new = {'x': [msa.objectives[0] for msa in self.points],
                   'y': [msa.objectives[1] for msa in self.points],
                   'seq_names': ['a' for _ in self.points]}

            # keep only new solutions, remove old ones
            self.source.stream(new, rollover=len(self.points))

    def __make_interactive_document(self, doc) -> None:
        fig = Figure(output_backend="webgl", sizing_mode='scale_width',
                     title=self.figure_title, tools=self.plot_tools)
        fig.circle(source=self.source, x='x', y='y', name='seq_names', size=self.plot_circle_size)

        # update plot
        doc.add_periodic_callback(self.__update_data_with_callback, self.period_milliseconds)

        # put the plot in a layout and add to the document
        doc.add_root(column(fig))
