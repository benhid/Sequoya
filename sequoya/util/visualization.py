import logging
from typing import TypeVar

from jmetal.util.visualization import InteractivePlot

LOGGER = logging.getLogger('Sequoya')
S = TypeVar('S')


class MSAPlot(InteractivePlot):

    def export_to_html(self, filename: str) -> None:
        html_string = '''
        <!DOCTYPE html>
        <html lang="en">
        
        <head>
          <meta charset="utf-8">
          <title>Sequoya dashboard</title>
          <!-- Plotly -->
          <script src="https://cdn.plot.ly/plotly-latest.min.js"></script>
          <!-- Bulma stylesheet -->
          <link href="https://cdnjs.cloudflare.com/ajax/libs/bulma/0.7.4/css/bulma.min.css" rel="stylesheet" />
          <style>
            .smenubar a.smenubar_alink {
              background: #00d1b2;
              border-radius: 5px;
              padding: 2px 20px;
              text-decoration: none;
              background-image: none;
              color: #fff;
              padding: 3px 10px;
              margin-left: 10px;
              text-decoration: none;
            }
        
            .box {
              overflow: hidden;
            }
        
            .biojs_msa_div {
              margin: 5px;
              overflow: hidden;
            }
          </style>
        </head>
        
        <body>
          <section class="hero is-primary">
            <div class="hero-body">
              <div class="container">
                <figure class="image is-128x128">
                  <img src="https://raw.githubusercontent.com/benhid/Sequoya/develop/docs/sequoya-white.png" />
                </figure>
              </div>
            </div>
          </section>
        
          <section class="section">
            <div class="container">
              ''' + self.export_to_div(include_plotlyjs=False) + '''
            </div>
          </section>
        
          <section class="section">
            <div class="container">
              <h1 class="subtitle">MSA viewer</h1>
              <div class="box">
                <div id="menuDiv"></div>
                <div id="rootDiv">(select solution from figure)</div>
              </div>
            </div>
          </section>
        
          <footer class="footer">
            <div class="content has-text-centered">
              <p>
                <strong>Sequoya</strong> by <a href="https://benhid.com">Antonio Benítez-Hidalgo</a>.
              </p>
            </div>
          </footer>
        
          <script src="http://cdn.bio.sh.s3.eu-central-1.amazonaws.com/msa/latest/msa.min.gz.js"></script>
          <script>
            var myPlot = document.getElementsByClassName('plotly-graph-div js-plotly-plot')[0];
            myPlot.on('plotly_click', function (data) {
              var pts = '';
        
              for (var i = 0; i < data.points.length; i++) {
                pts = '(x, y) = (' + data.points[i].x + ', ' + data.points[i].y.toPrecision(4) + ')';
                multiple_seq = data.points[i].customdata
                if (multiple_seq == undefined) multiple_seq = "";
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
