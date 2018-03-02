import logging
from functools import partial
from typing import List

from bokeh.application import Application
from bokeh.application.handlers import FunctionHandler
from bokeh.layouts import column
from bokeh.models import ColumnDataSource, HoverTool
from bokeh.server.server import Server
from bokeh.plotting import Figure

from pym2sa.core.solution import MSASolution

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class ScatterPlotMSA:

    def __init__(self, title: str, circle_size: int = 10, period_milliseconds: int = 50, tools: [] = None):
        if tools is None:
            tools = [HoverTool(tooltips=[("index", "$index"), ("(x,y)", "($x, $y)")]), 'save']

        self.plot_title = title
        self.plot_tools = tools
        self.plot_circle_size = circle_size
        self.period_milliseconds = period_milliseconds

        self.source = ColumnDataSource({'x': [], 'y': [], 'seq_names': []})
        self.population = []

    def interactive_plot(self, port: int = 5001) -> None:
        logger.info("Opening Bokeh application on http://localhost:{0}/".format(port))

        app = {'/pym2sa': Application(FunctionHandler(self.__make_document))}

        server = Server(app, port=port)
        server.show('/pym2sa')
        server.run_until_shutdown()

    '''
    def update_data(self, solutions: List[MSASolution]) -> None:
        new = {'x': [msa.objectives[0] for msa in solutions],
               'y': [msa.objectives[1] for msa in solutions],
               'seq_names': ['a' for msa in solutions]}

        # keep only new solutions, remove old ones
        self.source.stream(new, rollover=len(solutions))
    '''

    def replace_population(self, population: List[MSASolution]) -> None:
        self.population = population

    def __update_data_with_callback(self, population: List[MSASolution]) -> None:
        if population is not []:
            new = {'x': [msa.objectives[0] for msa in population],
                   'y': [msa.objectives[1] for msa in population],
                   'seq_names': ['a' for _ in population]}

            # keep only new solutions, remove old ones
            self.source.stream(new, rollover=len(population))

    def __make_document(self, doc) -> None:
        fig = Figure(output_backend="webgl", sizing_mode='scale_width',
                     title=self.plot_title, tools=self.plot_tools)
        fig.circle(source=self.source, x='x', y='y', name='seq_names', size=self.plot_circle_size)

        doc.add_periodic_callback(partial(self.__update_data_with_callback, population=self.population),
                                  self.period_milliseconds)

        # put the plot in a layout and add to the document
        doc.add_root(column(fig))
