import logging
from typing import List

import os
from jmetal.component.observer import BasicAlgorithmConsumer
from jmetal.operator.selection import BinaryTournament
from jmetal.util.comparator import RankingAndCrowdingDistanceComparator
from pymsa.core.score import PercentageOfNonGaps, PercentageOfTotallyConservedColumns

from pym2sa.algorithm.multiobjective.nsgaii import NSGA2MSA
from pym2sa.problem.dynamic_msa import MSA
from pym2sa.core.solution import MSASolution
from pym2sa.operators.crossover import SinglePointMSA
from pym2sa.operators.mutation import RandomGapInsertion
from pym2sa.util.graphic import ScatterPlotMSA

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def main() -> None:

    problem = MSA()
    problem.score_list = [PercentageOfNonGaps(), PercentageOfTotallyConservedColumns()]

    algorithm = NSGA2MSA[MSASolution, List[MSASolution]](
        problem=problem,
        population_size=100,
        initial_population_path=os.path.dirname(__file__)+'/files',
        max_evaluations=3000,
        mutation=RandomGapInsertion(probability=0.1),
        crossover=SinglePointMSA(probability=0.8),
        selection=BinaryTournament(comparator=RankingAndCrowdingDistanceComparator()))

    basic_consumer = BasicAlgorithmConsumer()
    algorithm.observable.register(observer=basic_consumer)

    algorithm.run()
    result = algorithm.get_result()

    # Plotting results
    plt = ScatterPlotMSA(plot_title=problem.get_name())
    plt.interactive_plot(solution_list=result)

    logger.info("Algorithm (MSA problem): " + algorithm.get_name())
    logger.info("Problem: " + problem.get_name())


if __name__ == '__main__':
    main()