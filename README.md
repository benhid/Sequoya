# Solving Multiple Sequence Alignments with Python
[![Build Status](https://img.shields.io/travis/benhid/Sequoya.svg?style=flat-square)](https://travis-ci.org/benhid/pyM2SA)
[![PyPI License](https://img.shields.io/pypi/l/Sequoya.svg?style=flat-square)]()
[![PyPI Python version](https://img.shields.io/pypi/pyversions/Sequoya.svg?style=flat-square)]()

Sequoya is an open source software tool aimed at for solving
*M*ultiple *S*equence *A*lignment problems with multi-objective metaheuristics.

This tool implements the [M2Align](https://github.com/KhaosResearch/M2Align) algorithm as shown in:

> "M2Align: parallel multiple sequence alignment with a multi-objective metaheuristic". Cristian Zambrano-Vega, Antonio J. Nebro José García-Nieto, José F. Aldana-Montes. Bioinformatics, Volume 33, Issue 19, 1 October 2017, Pages 3011–3017 ([DOI](https://doi.org/10.1093/bioinformatics/btx338)).

## Features
* The **scores** that are currently available are those from [pyMSA](https://github.com/benhid/pyMSA) (v0.5.1):
    * Sum of pairs,
    * Star,
    * Minimum entropy,
    * Percentage of non-gaps,
    * Percentage of totally conserved columns,
    * STRIKE.
* The **algorithm** that is currently available is:
    * NSGA-II
* **Crossover operator**:
    * Single-point crossover (`GapSequenceSolutionSinglePoint`).
* **Mutation operators**:
    * Shift closest gap group (`ShiftClosedGapGroups`),
    * Shift gap group (`ShiftGapGroup`),
    * Random gap insertion (`OneRandomGapInsertion`),
    * Merge two random adjacent gaps group (`TwoRandomAdjacentGapGroup`),
    * Multiple mutation (`MultipleMSAMutation`).

## Install
To download and install Sequoya just clone the Git repository hosted in GitHub:

```bash
$ git clone https://github.com/benhid/Sequoya.git
$ cd Sequoya
$ python setup.py install
```

Or via *pip*:

```bash
$ pip install Sequoya
```

## Usage
Examples of running Sequoya are located in the [`examples`](examples/) folder:

## Authors
### Active development team
* [Antonio Benítez-Hidalgo](https://benhid.github.io/about/) <antonio.b@uma.es>
* [Antonio J. Nebro](http://www.lcc.uma.es/%7Eantonio/) <antonio@lcc.uma.es>

## License
This project is licensed under the terms of the MIT - see the [LICENSE](LICENSE) file for details.
