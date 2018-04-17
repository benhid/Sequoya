<p align="center">
  <br/>
  <img src=resources/pym2sa.png alt="pyM2SA">
  <br/>
</p>

# Solving Multiple Sequence Alignments with Python
[![Build Status](https://travis-ci.com/benhid/pyM2SA.svg?token=6p1jsqj6w1XK5sV6iy3j&branch=master)](https://travis-ci.com/benhid/pyM2SA)
[![PyPI](https://img.shields.io/pypi/l/pyM2SA.svg)]()
[![PyPI](https://img.shields.io/pypi/v/pyM2SA.svg)]()


pyM<sup>2</sup>SA is an open source software tool aimed at for solving
*M*ultiple *S*equence *A*lignment problems with multi-objective metaheuristics.

## Features
* The **scores** that are currently available are those from [pyMSA](https://github.com/benhid/pyMSA) (v0.2):
    * Sum of pairs,
    * Star,
    * Minimum entropy,
    * Percentage of non-gaps,
    * Percentage of totally conserved columns,
    * STRIKE.
* The **algorithms** that are currently available are those from [jMetalPy](https://github.com/Metal/MetalPy) (v0.2):
    * NSGA-II,
    * SMPSO.
* **Crossover operator**:
    * Single-point crossover.
* **Mutation operator**:
    * Random gap insertion,
    * Closed gap shifting.

## Downloading
To download pyM<sup>2</sup>SA just clone the Git repository hosted in GitHub:

```bash
$ git clone https://github.com/benhid/pyM2SA.git
$ python setup.py install
```

## Requirements
pyM<sup>2</sup>SA has been developed with Python 3.6.0 :: [Anaconda](https://www.continuum.io) 4.3.0 (x86_64).

To install all dependencies use:

```bash
$ pip install -r requirements.txt
```

## Usage
Examples of running pyM<sup>2</sup>SA are located in the [`runner`](pym2sa/runner/) folder.

### Plotting solutions

It's possible to visualize in real time the _Pareto Frontier_ using [`ScatterPlot()`](pym2sa/util/graphic.py).
This class is able to deploy a [Bokeh](https://bokeh.pydata.org/en/latest/) Server to allow visualization in real-time using your web browser _or_ just plot the
final pareto frontier.

## Authors
### Active development team
* Antonio Benítez <antonio.benitez@lcc.uma.es>
* Antonio J. Nebro <antonio@lcc.uma.es>

## License
This project is licensed under the terms of the MIT - see the [LICENSE](LICENSE) file for details.
