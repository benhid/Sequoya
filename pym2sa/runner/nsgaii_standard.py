import logging
from typing import List

from jmetal.operator.selection import BinaryTournamentSelection
from jmetal.util.comparator import RankingAndCrowdingDistanceComparator
from jmetal.component.observer import VisualizerObserver, WriteFrontToFileObserver
from pymsa.core.score import SumOfPairs, PercentageOfTotallyConservedColumns
from pymsa.core.substitution_matrix import PAM250

from pym2sa.algorithm.multiobjective.nsgaii import NSGA2MSA
from pym2sa.component.observer import WriteSequencesToFileObserver
from pym2sa.core.solution import MSASolution
from pym2sa.problem.BalibaseMSA import BAliBaseMSA
from pym2sa.operators.crossover import SPXMSA
from pym2sa.operators.mutation import TwoRandomAdjacentGapGroup, ShiftGapGroup, MultipleMSAMutation
from pym2sa.util.graphic import ScatterMSA


def main() -> None:
    score_list = [SumOfPairs(PAM250()), PercentageOfTotallyConservedColumns()]
    problem = BAliBaseMSA(instance='BB12010', score_list=score_list)

    algorithm = NSGA2MSA[MSASolution, List[MSASolution]](
        problem=problem,
        population_size=100,
        max_evaluations=25000,
        mutation=MultipleMSAMutation([ShiftGapGroup(1.0), TwoRandomAdjacentGapGroup(1.0)], global_probability=0.2),
        crossover=SPXMSA(probability=0.8),
        selection=BinaryTournamentSelection(comparator=RankingAndCrowdingDistanceComparator())
    )

    algorithm.observable.register(observer=VisualizerObserver(problem))
    algorithm.observable.register(observer=WriteFrontToFileObserver('FUN_' + problem.get_name()))
    algorithm.observable.register(observer=WriteSequencesToFileObserver('VAR_' + problem.get_name()))

    algorithm.run()
    result = algorithm.get_result()

    pareto_front = ScatterMSA(plot_title='NSGAII for ' + problem.get_name(), number_of_objectives=problem.number_of_objectives,
                              xaxis_label='SOP', yaxis_label='TC')
    pareto_front.plot(result, output='plot-msa-' + problem.get_name())


if __name__ == '__main__':
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(threadName)-12.12s] [%(levelname)-5.5s]  %(message)s",
        handlers=[
            logging.FileHandler('jmetalpy.log'),
            logging.StreamHandler()
        ]
    )

    main()