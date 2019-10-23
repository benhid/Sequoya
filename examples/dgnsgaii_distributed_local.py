from dask.distributed import Client, LocalCluster
from jmetal.util.observer import ProgressBarObserver, PlotFrontToFileObserver
from jmetal.util.solutions.comparator import GDominanceComparator
from jmetal.util.termination_criterion import StoppingByEvaluations
from pymsa.core.score import SumOfPairs, PercentageOfTotallyConservedColumns

from sequoya.algorithm.multiobjective.nsgaii import DistributedNSGAII
from sequoya.operator import SPXMSA, ShiftClosedGapGroups
from sequoya.problem import BAliBASE

if __name__ == '__main__':
    # setup Dask client (web interface will be initialized at http://127.0.0.1:8787/workers)
    cluster = LocalCluster(n_workers=4, processes=True)
    client = Client(cluster)

    ncores = sum(client.ncores().values())
    print(f'{ncores} cores available')

    # creates the problem
    problem = BAliBASE(instance='BB12005', path='../resources',
                       score_list=[SumOfPairs(), PercentageOfTotallyConservedColumns()])

    # creates the algorithm
    max_evaluations = 20000
    reference_point = [4.15, 6300]

    algorithm = DistributedNSGAII(
        problem=problem,
        population_size=100,
        mutation=ShiftClosedGapGroups(probability=0.3),
        crossover=SPXMSA(probability=0.7),
        termination_criterion=StoppingByEvaluations(max=max_evaluations),
        dominance_comparator=GDominanceComparator(reference_point),
        number_of_cores=ncores,
        client=client
    )

    algorithm.observable.register(observer=ProgressBarObserver(max=max_evaluations))
    algorithm.observable.register(observer=PlotFrontToFileObserver(reference_point=[-6300, -4.15],
                                                                   output_directory='fronts_bb12005'))
    #algorithm.observable.register(observer=VisualizerObserver(reference_point=[-170000, -1.2]))

    algorithm.run()
    front = algorithm.get_result()

    print('Computing time: ' + str(algorithm.total_computing_time))
