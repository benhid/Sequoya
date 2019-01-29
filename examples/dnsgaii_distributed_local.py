from jmetal.operator import BinaryTournamentSelection
from pymsa.core.score import SumOfPairs, PercentageOfTotallyConservedColumns
from dask.distributed import Client, LocalCluster

from jmetal.util.comparator import RankingAndCrowdingDistanceComparator
from jmetal.util.observer import ProgressBarObserver
from jmetal.util.termination_criterion import StoppingByEvaluations
from jmetal.util.visualization import Plot
from pym2sa.algorithm.multiobjective.nsgaii import DistributedNSGAII
from pym2sa.problem import BAliBASE
from pym2sa.operator import SPXMSA, ShiftClosedGapGroups


if __name__ == '__main__':
    # Creates the problem
    problem = BAliBASE(balibase_instance='BB12001', balibase_path='../resources',
                       score_list=[SumOfPairs(), PercentageOfTotallyConservedColumns()])

    # Setup Dask client (web interface will be initialized at http://127.0.0.1:8787/workers)
    cluster = LocalCluster(n_workers=8, processes=True)
    client = Client(cluster)

    # Creates the algorithm
    max_evaluations = 10000
    algorithm = DistributedNSGAII(
        problem=problem,
        population_size=10,
        mutation=ShiftClosedGapGroups(probability=0.2),
        crossover=SPXMSA(probability=0.0),
        selection=BinaryTournamentSelection(comparator=RankingAndCrowdingDistanceComparator()),
        termination_criterion=StoppingByEvaluations(max=max_evaluations),
        number_of_cores=8,
        client=client
    )

    algorithm.observable.register(observer=ProgressBarObserver(max=max_evaluations))

    algorithm.run()
    front = algorithm.get_result()

    # Plot front
    plot_front = Plot(plot_title='Pareto front approximation', axis_labels=['%TC', '%SOP'])
    plot_front.plot(front, label='NSGAIII-BB12010', filename='NSGAIII-BB12010')
