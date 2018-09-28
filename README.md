<p align="center">
  <br/>
  <img src=resources/pym2sa-new.png alt="pyM2SA">
  <br/>
</p>

# Solving Multiple Sequence Alignments with Python
[![Build Status](https://img.shields.io/travis/benhid/pyM2SA.svg?style=flat-square)](https://travis-ci.org/benhid/pyM2SA)
[![Read the Docs](https://img.shields.io/readthedocs/pym2sa.svg?style=flat-square)](https://readthedocs.org/projects/pym2sa/)
[![PyPI License](https://img.shields.io/pypi/l/pyM2SA.svg?style=flat-square)]()
[![PyPI Python version](https://img.shields.io/pypi/pyversions/pyM2SA.svg?style=flat-square)]()

pyM<sup>2</sup>SA is an open source software tool aimed at for solving
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
* The **algorithms** that are currently available are:
    * NSGA-II
    * dNSGA-II
* **Crossover operator**:
    * Single-point crossover (`GapSequenceSolutionSinglePoint`).
* **Mutation operators**:
    * Shift closest gap group (`ShiftClosedGapGroups`),
    * Shift gap group (`ShiftGapGroup`),
    * Random gap insertion (`OneRandomGapInsertion`),
    * Merge two random adjacent gaps group (`TwoRandomAdjacentGapGroup`),
    * Multiple mutation (`MultipleMSAMutation`).

## Install
To download and install pyM<sup>2</sup>SA just clone the Git repository hosted in GitHub:

```bash
$ git clone https://github.com/benhid/pyM2SA.git
$ cd pyM2SA
$ python setup.py install
```

Or via *pip*:

```bash
$ pip install pym2sa
```

## Usage
Examples of running pyM<sup>2</sup>SA are located in the [`examples`](examples/) folder:

```python
# Defines the problem (as an instance from BAliBase 3.0)
score_list = [SumOfPairs(PAM250()), PercentageOfTotallyConservedColumns()]
problem = BAliBASE(instance='BB12010', balibase_path='resources', score_list=score_list)
problem.obj_labels = ['TC', 'SOP']

# Defines the algorithm
algorithm = NSGA2MSA[MSASolution, List[MSASolution]](
    problem=problem,
    population_size=100,
    max_evaluations=25000,
    mutation=MultipleMSAMutation([ShiftGapGroup(1.0), TwoRandomAdjacentGapGroup(1.0)], global_probability=0.2),
    crossover=SPXMSA(probability=0.8),
    selection=BinaryTournamentSelection(comparator=RankingAndCrowdingDistanceComparator())
)

# Register several observables
algorithm.observable.register(observer=WriteFrontToFileObserver('FUN_' + problem.instance))
algorithm.observable.register(observer=WriteSequencesToFileObserver('VAR_' + problem.instance))

# Run the algorithm
algorithm.run()
result = algorithm.get_result()

# Plot the solution
pareto_front = MSAPlot(plot_title='NSGAII for ' + problem.instance, axis_labels=problem.obj_labels)
pareto_front.plot(front)
pareto_front.to_html(filename='NSGAII-' + problem.instance)
```

<p align="center">
  <br/>
  <img src=resources/msaviewer.png alt="Interactive MSA viewer">
  <br/>
</p>

## Authors
### Active development team
* [Antonio Benítez-Hidalgo](https://benhid.github.io/about/) <antonio.b@uma.es>
* [Antonio J. Nebro](http://www.lcc.uma.es/%7Eantonio/) <antonio@lcc.uma.es>

## License
This project is licensed under the terms of the MIT - see the [LICENSE](LICENSE) file for details.
