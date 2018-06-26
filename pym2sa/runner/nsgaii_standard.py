import logging
from typing import List

from jmetal.operator.selection import BinaryTournamentSelection
from jmetal.util.comparator import RankingAndCrowdingDistanceComparator
from jmetal.component.observer import VisualizerObserver, WriteFrontToFileObserver

from pym2sa.algorithm.multiobjective.nsgaii import NSGA2MSA
from pym2sa.component.observer import WriteSequencesToFileObserver
from pym2sa.core.solution import MSASolution
from pym2sa.problem.BalibaseMSA import BAliBaseMSA
from pym2sa.operator.crossover import SPXMSA
from pym2sa.operator.mutation import TwoRandomAdjacentGapGroup, ShiftGapGroup, MultipleMSAMutation
from pym2sa.util.graphic import ScatterMSA


def main() -> None:
    problem = BAliBaseMSA(instance='BB12010')

    algorithm = NSGA2MSA[MSASolution, List[MSASolution]](
        problem=problem,
        population_size=100,
        max_evaluations=25000,
        mutation=MultipleMSAMutation([ShiftGapGroup(1.0), TwoRandomAdjacentGapGroup(1.0)], global_probability=0.2),
        crossover=SPXMSA(probability=0.8),
        selection=BinaryTournamentSelection(comparator=RankingAndCrowdingDistanceComparator())
    )

    algorithm.observable.register(observer=VisualizerObserver(problem))
    algorithm.observable.register(observer=WriteFrontToFileObserver("FUN"))
    algorithm.observable.register(observer=WriteSequencesToFileObserver("VAR"))

    algorithm.run()
    result = algorithm.get_result()

    # Plot final population
    pareto_front = ScatterMSA(plot_title='NSGAII for BB12010', number_of_objectives=problem.number_of_objectives,
                              xaxis_label='SOP', yaxis_label='TC')
    pareto_front.plot(result, output='output')


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
