from jmetal.component import RankingAndCrowdingDistanceComparator, VisualizerObserver
from jmetal.operator import BinaryTournamentSelection
from pymsa.core.score import SumOfPairs, PercentageOfTotallyConservedColumns
from dask.distributed import Client

from pym2sa.algorithm import dNSGAII
from pym2sa.problem import BAliBASE
from pym2sa.operator import SPXMSA, ShiftClosedGapGroups
from pym2sa.util.graphic import MSAPlot


def setup_distributed_client(address: str):
    new_client = Client(address)

    # This method will send (and import) the module up to all worker nodes in the cluster
    # Note: this file must be created by running `python setup.py install`
    new_client.upload_file('./pym2sa.egg')

    # Also, each worker should install the dependencies
    def install_dependencies_on_workers():
        import os
        os.system('pip install pymsa jmetalpy')
    new_client.run(install_dependencies_on_workers)

    return new_client


if __name__ == '__main__':
    # Creates the problem
    problem = BAliBASE(balibase_instance='BB12010', balibase_path='../resources',
                       score_list=[SumOfPairs(), PercentageOfTotallyConservedColumns()])
    problem.obj_labels = ['SOP', '%TC']

    # Setup Dask client
    client = setup_distributed_client('<dask-scheduler-ip>:8786')

    # Creates the algorithm
    algorithm = dNSGAII(
        problem=problem,
        population_size=100,
        max_evaluations=1000,
        mutation=ShiftClosedGapGroups(probability=0.2),
        crossover=SPXMSA(probability=0.8),
        selection=BinaryTournamentSelection(comparator=RankingAndCrowdingDistanceComparator()),
        number_of_cores=8,
        client=client
    )

    visualizer = VisualizerObserver()
    algorithm.observable.register(observer=visualizer)

    algorithm.run()
    front = algorithm.get_result()

    # Plot the solution
    pareto_front = MSAPlot(plot_title='NSGAII for ' + problem.balibase_instance, axis_labels=problem.obj_labels)
    pareto_front.plot(front)
    pareto_front.to_html(filename='NSGAII-' + problem.balibase_instance)

    print('Computing time: ' + str(algorithm.total_computing_time))