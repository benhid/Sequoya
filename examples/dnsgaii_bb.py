import logging

from dask.distributed import Client
from jmetal.operator.selection import BinaryTournamentSelection
from jmetal.util.comparator import RankingAndCrowdingDistanceComparator
from pymsa.core.score import SumOfPairs, PercentageOfTotallyConservedColumns
from pymsa.core.substitution_matrix import PAM250

from pym2sa.algorithm.multiobjective.dnsgaii import dNSGA2MSA
from pym2sa.problem.BalibaseMSA import BAliBASE
from pym2sa.operator.crossover import SPXMSA
from pym2sa.operator.mutation import TwoRandomAdjacentGapGroup, ShiftGapGroup, MultipleMSAMutation


def main() -> None:
    score_list = [SumOfPairs(PAM250()), PercentageOfTotallyConservedColumns()]
    problem = BAliBASE(instance='BB12010', score_list=score_list)

    algorithm = dNSGA2MSA(
        problem=problem,
        population_size=100,
        max_evaluations=10000,
        mutation=MultipleMSAMutation([ShiftGapGroup(1.0), TwoRandomAdjacentGapGroup(1.0)], probability=0.2),
        crossover=SPXMSA(probability=0.8),
        selection=BinaryTournamentSelection(comparator=RankingAndCrowdingDistanceComparator()),
        number_of_cores=8,
        client=Client('192.168.1.136:8786')
    )

    algorithm.run()
    print("Computing time: " + str(algorithm.total_computing_time))


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