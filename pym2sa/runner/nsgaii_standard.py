import os
import logging
from typing import List

from jmetal.operator.selection import BinaryTournamentSelection
from jmetal.util.comparator import RankingAndCrowdingDistanceComparator
from jmetal.util.solution_list_output import SolutionListOutput

from pym2sa.algorithm.multiobjective.nsgaii import NSGA2MSA
from pym2sa.component.observer import SpyPopulation
from pym2sa.problem.MSA import MSA
from pym2sa.core.solution import MSASolution
from pym2sa.operators.crossover import GapSequenceSolutionSinglePoint
from pym2sa.operators.mutation import OneRandomGapInsertion

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def main() -> None:
    problem = MSA(100)

    algorithm = NSGA2MSA[MSASolution, List[MSASolution]](
        problem=problem,
        population_size=100,
        initial_population_path=os.path.dirname(__file__)+'/dummy_files/',
        max_evaluations=1000,
        mutation=OneRandomGapInsertion(probability=0.1),
        crossover=GapSequenceSolutionSinglePoint(probability=0.8),
        selection=BinaryTournamentSelection(comparator=RankingAndCrowdingDistanceComparator()))
    #algorithm.observable.register(observer=SpyPopulation())
    algorithm.run()

    result = algorithm.get_result()
    SolutionListOutput[MSASolution].plot_scatter_to_screen(result)

    print(result[0].decode_alignment_as_list_of_pairs())

    logger.info("Algorithm: " + algorithm.get_name())
    logger.info("Problem: " + problem.get_name())


if __name__ == '__main__':
    main()