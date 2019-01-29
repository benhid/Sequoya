import logging
from typing import TypeVar

from jmetal.util.visualization import InteractivePlot

logger = logging.getLogger('pyM2SA')
S = TypeVar('S')


class MSAPlot(InteractivePlot):

    def export_to_html(self, filename: str) -> None:
        html_string = '''
        <!DOCTYPE html>
        <html lang="en">
            <head>
                <meta charset="utf-8">
                <title>pyM2SA</title>
                <script src="https://cdn.plot.ly/plotly-latest.min.js"></script>
                <script src="https://unpkg.com/sweetalert2@7.7.0/dist/sweetalert2.all.js"></script>
                <!-- include MSA js + css -->
                <script src="http://cdn.bio.sh.s3.eu-central-1.amazonaws.com/msa/latest/msa.min.gz.js"></script>
                <!-- Style -->
                <link href="https://fonts.googleapis.com/css?family=Open+Sans" rel="stylesheet">
                <link href="https://raw.githubusercontent.com/benhid/pyM2SA/gh-pages/lib/css/site.css" media="screen" rel="stylesheet" type="text/css" />
                <style>
                    * {
                      box-sizing: border-box
                    }

                    html {
                      -webkit-font-smoothing: antialiased;
                      -moz-osx-font-smoothing: grayscale;
                      font-family: 'Open Sans', sans-serif;
                      line-height: 1.6;
                      color: #666;
                      background: #F6F6F6;
                    }

                    img {
                      display: block;
                      max-width: 100%;
                    }

                    .logo{
                      margin: 0.5rem auto;
                      padding: 0.5rem 2.5rem;
                    }

                    h1 {
                      text-align: center;
                      padding: 0.5rem 2.5rem;
                      background-image: linear-gradient(120deg, #d64744 0%, #ff7c00 100%);
                      margin: 0 0 2rem 0;
                      font-size: 1.1rem;
                      color: white;
                    }


                    h2 {
                      text-align: center;
                      padding: 0.5rem 2.5rem;
                      background-image: linear-gradient(120deg, #990099 0%, #b72075 100%);
                      margin: 0 0 2rem 0;
                      font-size: 1.1rem;
                      color: white;
                    }


                    p {
                      padding: 0 2.5rem 1.5rem;
                      margin: 0;
                    }

                    .card {
                      margin: 1rem;
                      background: white;
                      box-shadow: 2px 8px 45px rgba(0, 0, 0, .15);
                      border-radius: 12px;
                      overflow: hidden;
                      transition: all .2s linear;
                    }

                    .card:hover {
                      box-shadow: 2px 8px 45px rgba(0, 0, 0, .20);
                    }

                    .btn {
                      border-radius: 5px;
                      padding: 2px 20px;
                      text-decoration: none;
                      color: #fff;
                      position: relative;
                      display: inline-block;
                    }

                    .btn:active {
                      transform: translate(0px, 3px);
                      -webkit-transform: translate(0px, 3px);
                      box-shadow: 0px 1px 0px 0px;
                    }

                    .orange {
                      background-image: linear-gradient(120deg, #e15631 0%, #ff7c00 100%);
                      box-shadow: 0px 3px 0px 0px #e15631;
                    }

                    .purple {
                      background-image: linear-gradient(120deg, #990099 0%, #b72075 100%);
                      box-shadow: 0px 3px 0px 0px #990099;
                    }

                    .grid-container {
                      margin: 2rem auto;
                      display: grid;
                      grid-template-columns: 100%; /*25% 10% auto auto;*/
                      grid-template-rows: auto;
                    }

                    .grid-container > div {
                    }

                    .float{
                      position:fixed;
                      right:40px;
                      bottom:40px;
                    }

                    .bk-root{
                        padding: 0 2.5rem 1.5rem;
                        margin: 0;
                    }

                    .msaviewer{
                        padding: 0 2.5rem 1.5rem;
                        margin: 0;
                    }

                    .smenubar .smenubar_alink {
                        background: #b72075;
                        border-radius: 5px;
                        padding: 2px 20px;
                        text-decoration: none;
                        background-image: none;
                        color: #fff;
                        padding: 3px 10px;
                        margin-left: 10px;
                        text-decoration: none;
                    }
                </style>
            </head>
            <body>
                <div class="grid-container">
                    <div class="item" style="align-self: center">
                        <a href="https://benhid.github.io/pyM2SA/"><img src="https://raw.githubusercontent.com/benhid/pyM2SA/gh-pages/lib/img/sequoya.png" class="logo"/></a>
                    </div>
                </div>

                <div class="container">
                    <div class="card">
                        <h1>Plot</h1>
                        ''' + self.export_to_div(include_plotlyjs=False) + '''
                    </div>

                    <div class="card">
                        <h2>MSA Viewer</h2>
                        <p>
                            <div id="menuDiv" class="msaviewer"></div>
                            <div id="rootDiv" class="msaviewer">Select solution.</div>
                        </p>
                    </div>
                </div>

                <script>                
                    var myPlot = document.getElementsByClassName('plotly-graph-div js-plotly-plot')[0];
                    myPlot.on('plotly_click', function(data){
                        var pts = '';

                        for(var i=0; i < data.points.length; i++){
                            pts = '(x, y) = ('+data.points[i].x +', '+ data.points[i].y.toPrecision(4)+')';
                            multiple_seq = data.points[i].customdata
                            if(multiple_seq == undefined) multiple_seq = "";
                        }

                        // read msa            
                        var opts = {
                            el: document.getElementById("rootDiv"),
                            seqs: msa.io.fasta.parse(multiple_seq),
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
                    });
                </script>
            </body>
        </html>'''

        with open(filename + '.html', 'w') as outf:
            outf.write(html_string)
