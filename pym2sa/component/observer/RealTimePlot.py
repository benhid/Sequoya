import logging

from jmetal.util.observable import Observer

from pym2sa.util.graphic import ScatterPlot

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class RealTimePlot(Observer):
    def __init__(self, title: str, ) -> None:
        self.figure = ScatterPlot(title)
        self.figure.interactive_plot()

    def update(self, *args, **kwargs):
        new_population = kwargs["population"]
        self.figure.replace_points(new_population)
