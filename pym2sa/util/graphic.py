import logging
import json
from typing import TypeVar, List
from os.path import dirname, join

from bokeh.embed import file_html
from bokeh.resources import CDN
from bokeh.client import ClientSession
from bokeh.io import curdoc, reset_output
from bokeh.layouts import column, row
from bokeh.models import HoverTool, ColumnDataSource, TapTool, CustomJS, WheelZoomTool, Markup
from bokeh.plotting import Figure
from jinja2 import Environment, FileSystemLoader

from jmetal.util.graphic import Plot


logger = logging.getLogger(__name__)
S = TypeVar('S')

BASE_PATH = dirname(join(dirname(__file__)))


class ScatterMSA(Plot):

    def __init__(self, plot_title: str, number_of_objectives: int,
                 xaxis_label: str='', yaxis_label: str='', ws_url: str='localhost:5006'):
        super().__init__(plot_title, number_of_objectives, xaxis_label=xaxis_label, yaxis_label=yaxis_label)

        if self.number_of_objectives == 2:
            self.source = ColumnDataSource(data=dict(x=[], y=[], str=[]))
        else:
            raise Exception('Wrong number of objectives: {0}'.format(number_of_objectives))

        self.client = ClientSession(websocket_url="ws://{0}/ws".format(ws_url))
        self.doc = curdoc()
        self.doc.title = plot_title
        self.figure_xy = None

        self.__initialize()

    def __initialize(self) -> None:
        """ Set-up tools for plot. """
        code = '''
            selected = source.selected['1d']['indices'][0]
            var str = source.data.str[selected]
            console.log(str)
            
            // read msa            
            var opts = {
                el: document.getElementById("rootDiv"),
                seqs: msa.io.fasta.parse(str),
                vis: {
                    conserv: false,
                    overviewbox: false,
                    seqlogo: true
                },
                conf: {
                    dropImport: true,
                    debug: false,
                },
                zoomer: {
                    menuFontsize: "12px",
                    autoResize: true
                }
            };

            // init msa
            var m = new msa.msa(opts);
            renderMSA();
            
            function renderMSA() {
                // the menu is independent to the MSA container
                var menuOpts = {
                    el: document.getElementById('menuDiv'),
                    msa: m
                };
                var defMenu = new msa.menu.defaultmenu(menuOpts);
                m.addView("menu", defMenu);
    
                // call render at the end to display the whole MSA
                m.render();
            }
        '''

        callback = CustomJS(args=dict(source=self.source), code=code)
        self.plot_tools = [TapTool(callback=callback), WheelZoomTool(), 'save', 'pan',
                           HoverTool(tooltips=[("index", "$index"), ("(x,y)", "($x, $y)")])]

    def plot(self, solution_list: List[S], output: str='') -> None:
        # This is important to purge data (if any) between calls
        reset_output()

        # Set up figure
        self.figure_xy = Figure(output_backend='webgl', sizing_mode='scale_width', title=self.plot_title,
                                tools=self.plot_tools)
        self.figure_xy.scatter(x='x', y='y', legend='solution', fill_alpha=0.7, source=self.source)
        self.figure_xy.xaxis.axis_label = self.xaxis_label
        self.figure_xy.yaxis.axis_label = self.yaxis_label

        x_values, y_values, _ = self.get_objectives(solution_list)

        # Push data to server
        self.source.stream({'x': [-v for v in x_values], 'y': [-v for v in y_values],
                            'str': [s.__str__() for s in solution_list]})
        self.doc.add_root(column(self.figure_xy))

        self.client.push(self.doc)

        if output:
            self.__save(output)

    def __save(self, file_name: str):
        env = Environment(loader=FileSystemLoader(BASE_PATH + '/util/'))
        env.filters['json'] = lambda obj: Markup(json.dumps(obj))

        html = file_html(models=self.doc, resources=CDN, template=env.get_template('file.html'))
        with open(file_name + '.html', 'w') as of:
            of.write(html)

    def disconnect(self):
        if self.is_connected():
            self.client.close()

    def is_connected(self) -> bool:
        return self.client.connected
