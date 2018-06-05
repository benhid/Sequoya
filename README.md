<p align="center">
  <br/>
  <img src=resources/pym2sa.png alt="pyM2SA">
  <br/>
</p>

# Solving Multiple Sequence Alignments with Python
[![Build Status](https://travis-ci.com/benhid/pyM2SA.svg?token=6p1jsqj6w1XK5sV6iy3j&branch=master?style=flat-square)](https://travis-ci.com/benhid/pyM2SA)
[![PyPI License](https://img.shields.io/pypi/l/pyM2SA.svg?style=flat-square)]()
[![PyPI Python version](https://img.shields.io/pypi/pyversions/pyM2SA.svg?style=flat-square)]()

pyM<sup>2</sup>SA is an open source software tool aimed at for solving
*M*ultiple *S*equence *A*lignment problems with multi-objective metaheuristics.

This tool implements the [M2Align](https://github.com/KhaosResearch/M2Align) algorithm as shown in:

> "M2Align: parallel multiple sequence alignment with a multi-objective metaheuristic". Cristian Zambrano-Vega, Antonio J. Nebro José García-Nieto, José F. Aldana-Montes. Bioinformatics, Volume 33, Issue 19, 1 October 2017, Pages 3011–3017 ([DOI](https://doi.org/10.1093/bioinformatics/btx338)).

## Features
* The **scores** that are currently available are those from [pyMSA](https://github.com/benhid/pyMSA) (v1.0.0):
    * Sum of pairs,
    * Star,
    * Minimum entropy,
    * Percentage of non-gaps,
    * Percentage of totally conserved columns,
    * STRIKE.
* The **algorithms** that are currently available are those from [jMetalPy](https://github.com/Metal/MetalPy) (v1.0.0):
    * NSGA-II,
    * SMPSO.
* **Crossover operator**:
    * Single-point crossover (`GapSequenceSolutionSinglePoint`).
* **Mutation operator**:
    * Shift gap group (`ShiftGapGroup`),
    * Random gap insertion (`OneRandomGapInsertion`),
    * Merge two random adjacent gaps group (`TwoRandomAdjacentGapGroup`).

## Install
To download and install pyM<sup>2</sup>SA just clone the Git repository hosted in GitHub:

```bash
$ git clone https://github.com/benhid/pyM2SA.git
$ python setup.py install
```

Or via pip:

```bash
$ pip install pym2sa
```

## Requirements
pyM<sup>2</sup>SA has been developed with Python 3.6.0 :: [Anaconda](https://www.continuum.io) 4.3.0 (x86_64).

To install all dependencies use:

```bash
$ pip install -r requirements.txt
```

## Usage
Examples of running pyM<sup>2</sup>SA are located in the [`runner`](pym2sa/runner/) folder.

## Authors
### Active development team
* Antonio Benítez-Hidalgo <antonio.b@uma.es>
* Antonio J. Nebro <antonio@lcc.uma.es>

## License
This project is licensed under the terms of the MIT - see the [LICENSE](LICENSE) file for details.
