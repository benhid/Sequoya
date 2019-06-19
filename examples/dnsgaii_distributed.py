import matplotlib

matplotlib.use('TkAgg')

from dask.distributed import Client
from jmetal.operator import BinaryTournamentSelection
from jmetal.util.comparator import RankingAndCrowdingDistanceComparator
from jmetal.util.observer import ProgressBarObserver, VisualizerObserver
from jmetal.util.termination_criterion import StoppingByEvaluations
from jmetal.util.visualization import Plot
from pymsa.core.score import SumOfPairs, PercentageOfTotallyConservedColumns

from sequoya.algorithm.multiobjective.nsgaii import DistributedNSGAII
from sequoya.operator import SPXMSA, ShiftClosedGapGroups
from sequoya.problem import BAliBASE
from sequoya.util.visualization import MSAPlot

if __name__ == '__main__':
    # setup Dask client
    client = Client('192.168.213.3:8786')
    client.restart()

    ncores = sum(client.ncores().values())
    print(f'{ncores} cores available on cluster')

    # creates the problem
    problem = BAliBASE(instance='BB12010', path='../resources',
                       score_list=[SumOfPairs(), PercentageOfTotallyConservedColumns()])

    # creates the algorithm
    max_evaluations = 25000

    algorithm = DistributedNSGAII(
        problem=problem,
        population_size=100,
        mutation=ShiftClosedGapGroups(probability=0.3),
        crossover=SPXMSA(probability=0.7),
        selection=BinaryTournamentSelection(comparator=RankingAndCrowdingDistanceComparator()),
        termination_criterion=StoppingByEvaluations(max=max_evaluations),
        number_of_cores=ncores,
        client=client
    )

    algorithm.observable.register(observer=ProgressBarObserver(max=max_evaluations))
    #algorithm.observable.register(observer=VisualizerObserver())

    algorithm.run()
    front = algorithm.get_result()

    # plot front
    plot_front = Plot(plot_title='Pareto front approximation', axis_labels=['%SOP', '%TC'])
    plot_front.plot(front, label='NSGAIII-BB50012', filename='NSGAIII-BB50012')

    # plot interactive front
    pareto_front = MSAPlot(plot_title='Pareto front approximation', axis_labels=['%SOP', '%TC'])
    pareto_front.plot(front, label='NSGAIII-BB50011', filename='NSGAIII-BB50011')

    print('Computing time: ' + str(algorithm.total_computing_time))
