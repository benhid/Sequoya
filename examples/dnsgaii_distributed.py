from dask.distributed import Client
from jmetal.operator import BinaryTournamentSelection
from jmetal.util.comparator import RankingAndCrowdingDistanceComparator
from jmetal.util.observer import ProgressBarObserver
from jmetal.util.termination_criterion import StoppingByEvaluations
from jmetal.util.visualization import Plot
from pymsa.core.score import SumOfPairs, PercentageOfTotallyConservedColumns

from sequoya.algorithm.multiobjective.nsgaii import DistributedNSGAII
from sequoya.operator import SPXMSA, ShiftClosedGapGroups
from sequoya.problem import BAliBASE
from sequoya.util.visualization import MSAPlot


def install_dependencies():
    """ Install packages on worker nodes. Note that this may take a while (the first time)! """
    import os
    os.system('pip install pymsa jmetalpy sequoya')


if __name__ == '__main__':
    # creates the problem
    problem = BAliBASE(balibase_instance='BB20019', balibase_path='../resources',
                       score_list=[SumOfPairs(), PercentageOfTotallyConservedColumns()])

    # creates the Dask client
    client = Client('127.0.0.1:8786')
    client.run(install_dependencies)

    # creates the algorithm
    max_evaluations = 1000

    algorithm = DistributedNSGAII(
        problem=problem,
        population_size=10,
        mutation=ShiftClosedGapGroups(probability=0.2),
        crossover=SPXMSA(probability=0.7),
        selection=BinaryTournamentSelection(comparator=RankingAndCrowdingDistanceComparator()),
        termination_criterion=StoppingByEvaluations(max=max_evaluations),
        number_of_cores=4,
        client=client
    )

    algorithm.observable.register(observer=ProgressBarObserver(max=max_evaluations))

    algorithm.run()
    front = algorithm.get_result()

    # plot front
    plot_front = Plot(plot_title='Pareto front approximation', axis_labels=['%SOP', '%TC'])
    plot_front.plot(front, label='NSGAIII-BB20019', filename='NSGAIII-BB20019')

    # plot interactive front
    pareto_front = MSAPlot(plot_title='Pareto front approximation', axis_labels=['%SOP', '%TC'])
    pareto_front.plot(front, label='NSGAIII-BB20019', filename='NSGAIII-BB20019')

    print('Computing time: ' + str(algorithm.total_computing_time))
