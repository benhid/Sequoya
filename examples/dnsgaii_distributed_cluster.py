from jmetal.component import RankingAndCrowdingDistanceComparator
from jmetal.operator import BinaryTournamentSelection
from pymsa.core.score import SumOfPairs, PercentageOfTotallyConservedColumns
from dask.distributed import Client

from pym2sa.algorithm import dNSGA2MSA
from pym2sa.problem import BAliBASE
from pym2sa.operator import SPXMSA, ShiftClosedGapGroups
from pym2sa.util.graphic import MSAPlot


if __name__ == '__main__':
    # Creates the problem
    problem = BAliBASE(instance='BB12001', balibase_path='../resources',
                       score_list=[SumOfPairs(), PercentageOfTotallyConservedColumns()])
    problem.obj_labels = ['TC', 'SOP']

    # Setup Dask client
    client = Client('<dask-scheduler-ip>:8786')

    # This method will send (and import) the module up to all worker nodes in the cluster
    # Note: this file must be created by running `python setup.py bdist_egg`
    client.upload_file('./pym2sa.egg')

    # Creates the algorithm
    algorithm = dNSGA2MSA(
        problem=problem,
        max_evaluations=25000,
        mutation=ShiftClosedGapGroups(probability=0.2),
        crossover=SPXMSA(probability=0.8),
        selection=BinaryTournamentSelection(comparator=RankingAndCrowdingDistanceComparator()),
        number_of_cores=100,
        client=client
    )

    algorithm.run()
    front = algorithm.get_result()

    # Plot the solution
    pareto_front = MSAPlot(plot_title='NSGAII for ' + problem.instance, axis_labels=problem.obj_labels)
    pareto_front.plot(front)
    pareto_front.to_html(filename='NSGAII-' + problem.instance)

    print('Computing time: ' + str(algorithm.total_computing_time))