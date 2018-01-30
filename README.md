<p align="center">
  <br/>
  <img src=resources/pym2sa.png alt="PyM2SA">
  <br/>
</p>

# Scoring Multiple Sequence Alignments with Python

PyM<sup>2</sup>SA is an open source software tool aimed at for solving
*M*ultiple *S*equence *A*lignment problems with multi-objective metaheuristics.

## Features
* The **scores** that are currently available are those from [pyMSA](https://github.com/benhid/pyMSA):
    * Sum of pairs,
    * Star,
    * Minimum entropy,
    * Percentage of non-gaps,
    * Percentage of totally conserved columns,
    * STRIKE
* The **algorithms** that are currently available are:
    * NSGA-II,
    * SMPSO
* **Crossover operator**:
    * Single-point crossover
* **Mutation operator**:
    * Random gap insertion

## Downloading
To download PyM<sup>2</sup>SA just clone the Git repository hosted in GitHub:

```bash
$ git clone https://github.com/benhid/pyM2SA.git
$ python setup.py install
```

## Requirements
PyM<sup>2</sup>SA has been developed with Python 3.6.0 :: [Anaconda](https://www.continuum.io) 4.3.0 (x86_64).

To install all dependencies use:

```bash
$ pip install -r requirements.txt
```

## Usage
Examples of running PyM<sup>2</sup>SA are located in the [`runner`](pym2sa/runner/) folder.

## Authors
### Active development team
* Antonio Benítez <antonio.benitez@lcc.uma.es>

## License
This project is licensed under the terms of the MIT - see the [LICENSE](LICENSE) file for details.
