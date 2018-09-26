from jmetal.component import RankingAndCrowdingDistanceComparator, VisualizerObserver
from jmetal.operator import BinaryTournamentSelection
from pymsa.core.score import SumOfPairs, PercentageOfTotallyConservedColumns, PercentageOfNonGaps
from dask.distributed import Client, LocalCluster

from pym2sa.algorithm import dNSGA2BAliBASE
from pym2sa.problem import BAliBASE
from pym2sa.operator import SPXMSA, TwoRandomAdjacentGapGroup
from pym2sa.util.graphic import MSAPlot


if __name__ == '__main__':
    # Creates the problem
    problem = BAliBASE(instance='BB12010', balibase_path='../resources',
                       score_list=[SumOfPairs(), PercentageOfTotallyConservedColumns(), PercentageOfNonGaps()])
    problem.obj_labels = ['SOP', '%TC', '%NonGaps']

    # Setup Dask client (web interface will be initialized at http://127.0.0.1:8787/workers)
    cluster = LocalCluster(n_workers=8, processes=True)
    client = Client(cluster)

    # Creates the algorithm
    algorithm = dNSGA2BAliBASE(
        problem=problem,
        population_size=100,
        max_evaluations=1000,
        mutation=TwoRandomAdjacentGapGroup(probability=0.2),
        crossover=SPXMSA(probability=0.8),
        selection=BinaryTournamentSelection(comparator=RankingAndCrowdingDistanceComparator()),
        number_of_cores=4,
        client=client
    )

    visualizer = VisualizerObserver()
    algorithm.observable.register(observer=visualizer)

    algorithm.run()
    front = algorithm.get_result()

    # Plot the solution
    pareto_front = MSAPlot(plot_title='NSGAII for ' + problem.instance, axis_labels=problem.obj_labels)
    pareto_front.plot(front)
    pareto_front.to_html(filename='NSGAII-' + problem.instance)

    print('Computing time: ' + str(algorithm.total_computing_time))