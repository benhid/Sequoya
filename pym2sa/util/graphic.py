import os
import webbrowser
import logging
from typing import List, TypeVar

from jmetal.util.graphic import ScatterPlot

from pym2sa.core.solution import MSASolution

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

S = TypeVar('S')
R = TypeVar(List[S])


class ScatterPlotMSA(ScatterPlot):

    def retrieve_info(self, solution: MSASolution):
        """ Let us specify the information that we want to retrieve from the points by clicking. """
        logger.info('Saving into file...')

        file_name = "solution-" + str(solution.objectives[0]) + "-" + str(solution.objectives[1])

        # Save solution to file
        with open(file_name + ".txt", 'w') as output:
            for i in range(0, solution.number_of_variables):
                output.write('>{0}\n{1}\n'.format(solution.header[i], solution.variables[i]))
            output.write('Scores: Sum of pairs = {0}, Percentage of non-gaps = {1}%'
                         .format(solution.objectives[0], solution.objectives[1]))

        # Save solution to .html
        with open(file_name + ".html", 'w') as output:
            fasta = ""

            for i in range(0, solution.number_of_variables):
                fasta = fasta + ('>{0}\n{1}\n'.format(solution.header[i], solution.variables[i]))

            output.write("""
                            <!DOCTYPE html>
                            <html lang="en">
                                <head>
                                    <meta charset="UTF-8">
                                    <meta name="description" content="MSA Viewer" />

                                    <!-- include MSA js + css -->
                                    <script src="http://cdn.bio.sh.s3.eu-central-1.amazonaws.com/msa/latest/msa.min.gz.js"></script>
                                    <!-- Bootstrap Core CSS -->
                                    <link rel="stylesheet" href="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.7/css/bootstrap.min.css" integrity="sha384-BVYiiSIFeK1dGmJRAkycuHAHRg32OmUcww7on3RYdg4Va+PmSTsz/K68vbdEjh4u" crossorigin="anonymous">
                                    <style>
                                        body {
                                            display: block;
                                            margin: 8px;
                                        }
                                    </style>
                                </head>

                                <body>
                                    <div class="container-fluid">
                                        <div class="page-header">
                                            <h1>MSA <small>Multiple Sequence Alignment</small></h1>
                                        </div>
                                        <div class="alert alert-info" role="alert">
                                          <span class="glyphicon glyphicon-flash" aria-hidden="true"></span> This document was automatically generated by <a href="https://github.com/jMetal/jMetalPy" class="alert-link">jMetalPy</a>. It represents one solution of a MSA problem.
                                        </div>

                                        <div class="panel panel-info"> 
                                            <div class="panel-heading"> 
                                                <h3 class="panel-title">Scores</h3> 
                                            </div> 
                                            <div class="panel-body"> 
                                                <p>Sum of pairs: <span class="badge">""" + str(solution.objectives[0]) + """</span></p>
                                                <p>Percentage of non-gaps: <span class="badge">""" + str(solution.objectives[1]) + """</span></p>
                                            </div>
                                        </div>

                                        <div class="panel panel-success">
                                            <div class="panel-heading">
                                                <div class="row">
                                                  <div class="col col-lg-6 col-md-6 col-xs-6">
                                                    <h3 class="panel-title">Viewer</h3>
                                                  </div>
                                                  <div class="col col-lg-6 col-md-6 col-xs-6">
                                                    <a id="download" download="image.png">
                                                      <button type="button" class="btn btn-default btn-xs pull-right" onclick="toPNG()">Export MSA image (as PNG)</button>
                                                    </a>
                                                  </div>
                                                </div>
                                            </div>
                                            <div class="panel-body">
                                                <div id="msa" style="overflow-y: auto; overflow-x: hidden;">Loading solution... </div>
                                            </div>
                                        </div>
                                    </div>

                                    <pre style="display: none" id="fasta-file">""" + fasta + """</pre>

                                    <script>
                                        var fasta = document.getElementById("fasta-file").innerText;

                                        var opts = {
                                            el: document.getElementById("yourDiv"),
                                            vis: {
                                                conserv: true,
                                                seqlogo: true
                                            },
                                            conf: {
                                                dropImport: true
                                            },
                                            el: document.getElementById("msa"),
                                            seqs: msa.io.fasta.parse(fasta)
                                        };

                                        // init msa
                                        var m = new msa.msa(opts);
                                        m.render();
                                    </script>
                                    <script>
                                        function toPNG(){
                                          var download = document.getElementById("download");
                                          var sequences = document.getElementsByClassName("biojs_msa_seqblock")[0].toDataURL("image/png")
                                                      .replace("image/png", "image/octet-stream");
                                          download.setAttribute("href", sequences);
                                          download.click;
                                        }
                                    </script>
                                </body>
                            </html>
                        """)

            logger.info("Opening .HTML ...")

            url = "file:///" + os.path.abspath(os.path.join(os.path.dirname(__file__), '..')) \
                  + "\\" + file_name + ".html"

            webbrowser.open(url, new=2)