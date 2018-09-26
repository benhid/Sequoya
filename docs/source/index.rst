pyM2SA: Solving Multiple Sequence Alignments with Python
========================================================

.. warning:: Documentation is WIP!! Some information may be missing.

.. toctree::
   :maxdepth: 1
   :caption: Contents:

   about
   api/pym2sa

Installation steps
------------------

Via pip:

.. code-block:: console

    $ pip install pym2sa

Via Github:

.. code-block:: console

    $ git clone https://github.com/benhid/pyM2SA.git
    $ cd pyM2SA
    $ python setup.py install

Features
--------

* The **scores** that are currently available are those from [pyMSA](https://github.com/benhid/pyMSA) (v0.5.0):
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