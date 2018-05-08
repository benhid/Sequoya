import logging
from typing import List

from jmetal.operator.selection import BinaryTournamentSelection
from jmetal.util.comparator import RankingAndCrowdingDistanceComparator
from jmetal.util.solution_list_output import SolutionListOutput
from jmetal.component.observer import AlgorithmObserver

from pym2sa.algorithm.multiobjective.nsgaii import NSGA2MSA
from pym2sa.component.observer import SpyPopulation
from pym2sa.problem.BalibaseMSA import BalisebaseMSA
from pym2sa.core.solution import MSASolution
from pym2sa.operators.crossover import GapSequenceSolutionSinglePoint
from pym2sa.operators.mutation import OneRandomGapInsertion, TwoRandomAdjacentGapGroup

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def main() -> None:
    problem = BalisebaseMSA(instance='BB11006', number_of_variables=100)

    algorithm = NSGA2MSA[MSASolution, List[MSASolution]](
        problem=problem,
        population_size=100,
        max_evaluations=25000,
        mutation=TwoRandomAdjacentGapGroup(probability=0.2),
        crossover=GapSequenceSolutionSinglePoint(probability=0.9),
        selection=BinaryTournamentSelection(comparator=RankingAndCrowdingDistanceComparator()))

    algorithm.observable.register(observer=AlgorithmObserver(1*10e-3))
    algorithm.observable.register(observer=SpyPopulation())
    algorithm.run()

    result = algorithm.get_result()
    SolutionListOutput[MSASolution].plot_scatter_to_screen(result)

    print(result[0].decode_alignment_as_list_of_pairs())

    logger.info("Algorithm: " + algorithm.get_name())
    logger.info("Problem: " + problem.get_name())


if __name__ == '__main__':
    main()