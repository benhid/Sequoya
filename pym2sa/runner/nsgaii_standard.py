import os
import logging
from typing import List

from jmetal.component.observer import AlgorithmObserver
from jmetal.operator.selection import BinaryTournament
from jmetal.util.comparator import RankingAndCrowdingDistanceComparator
from pymsa.core.score import PercentageOfNonGaps, PercentageOfTotallyConservedColumns

from pym2sa.algorithm.multiobjective.nsgaii import NSGA2MSA
from pym2sa.component.observer import RealTimePlot
from pym2sa.problem.dynamic_msa import MSA
from pym2sa.core.solution import MSASolution
from pym2sa.operators.crossover import GapSequenceSolutionSinglePoint
from pym2sa.operators.mutation import OneRandomGapInsertion

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def main() -> None:

    problem = MSA()
    problem.score_list = [PercentageOfNonGaps(), PercentageOfTotallyConservedColumns()]

    algorithm = NSGA2MSA[MSASolution, List[MSASolution]](
        problem=problem,
        population_size=100,
        initial_population_path=os.path.dirname(__file__)+'/dummy_files/',
        max_evaluations=30000,
        mutation=OneRandomGapInsertion(probability=0.1),
        crossover=GapSequenceSolutionSinglePoint(probability=0.8),
        selection=BinaryTournament(comparator=RankingAndCrowdingDistanceComparator()))

    #graphic_consumer = RealTimePlot(title="NSGA-II")
    graphic_consumer = AlgorithmObserver(animation_speed=1 * 10e-8)
    algorithm.observable.register(observer=graphic_consumer)

    algorithm.run()

    logger.info("Algorithm (MSA problem): " + algorithm.get_name())
    logger.info("Problem: " + problem.get_name())


if __name__ == '__main__':
    main()