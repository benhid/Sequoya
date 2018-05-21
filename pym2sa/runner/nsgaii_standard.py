import logging
from typing import List

from jmetal.operator.selection import BinaryTournamentSelection
from jmetal.util.comparator import RankingAndCrowdingDistanceComparator
from jmetal.component.observer import VisualizerObserver, WriteFrontToFileObserver
from jmetal.util.solution_list_output import SolutionListOutput

from pym2sa.algorithm.multiobjective.nsgaii import NSGA2MSA
from pym2sa.core.solution import MSASolution
from pym2sa.problem.BalibaseMSA import BalisebaseMSA
from pym2sa.operators.crossover import SPXMSA
from pym2sa.operators.mutation import OneRandomGapInsertion, TwoRandomAdjacentGapGroup

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def main() -> None:
    problem = BalisebaseMSA(instance='BB12024', number_of_variables=100)

    algorithm = NSGA2MSA[MSASolution, List[MSASolution]](
        problem=problem,
        population_size=100,
        max_evaluations=25000,
        mutation=TwoRandomAdjacentGapGroup(probability=0.2),
        crossover=SPXMSA(probability=0.9, remove_gap_columns=True),
        selection=BinaryTournamentSelection(comparator=RankingAndCrowdingDistanceComparator()))

    algorithm.observable.register(observer=VisualizerObserver(1*10e-3))
    algorithm.observable.register(observer=WriteFrontToFileObserver("FUN"))
    algorithm.run()

    result = algorithm.get_result()
    SolutionListOutput[MSASolution].plot_frontier_to_screen(result)


    logger.info("Algorithm: " + algorithm.get_name())
    logger.info("Problem: " + problem.get_name())


if __name__ == '__main__':
    main()
