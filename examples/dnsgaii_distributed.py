from dask.distributed import Client
from jmetal.lab.visualization import Plot
from jmetal.util.observer import ProgressBarObserver, VisualizerObserver
from jmetal.util.termination_criterion import StoppingByEvaluations
from pymsa.core.score import SumOfPairs, PercentageOfTotallyConservedColumns

from sequoya.algorithm.multiobjective.nsgaii import DistributedNSGAII
from sequoya.operator import SPXMSA, ShiftClosedGapGroups
from sequoya.problem import BAliBASE
from sequoya.util.solution import get_representative_set
from sequoya.util.visualization import MSAPlot

if __name__ == '__main__':
    # setup Dask client
    client = Client('<ip>:<port>')
    client.restart()

    ncores = sum(client.ncores().values())
    print(f'{ncores} cores available on cluster')

    # creates the problem
    problem = BAliBASE(instance='BB50011', path='../resources',
                       score_list=[SumOfPairs(), PercentageOfTotallyConservedColumns()])

    # creates the algorithm
    max_evaluations = 25000

    algorithm = DistributedNSGAII(
        problem=problem,
        population_size=100,
        mutation=ShiftClosedGapGroups(probability=0.3),
        crossover=SPXMSA(probability=0.7),
        termination_criterion=StoppingByEvaluations(max_evaluations=max_evaluations),
        number_of_cores=ncores,
        client=client
    )

    algorithm.observable.register(observer=ProgressBarObserver(max=max_evaluations))
    algorithm.observable.register(observer=VisualizerObserver())

    algorithm.run()
    front = algorithm.get_result()

    # plot front
    plot_front = Plot(title='Pareto front approximation', axis_labels=['%SOP', '%TC'])
    plot_front.plot(front, label='NSGAII-BB50011', filename='NSGAII-BB50011')

    pareto_front = MSAPlot(title='Pareto front approximation', axis_labels=['%SOP', '%TC'])
    pareto_front.plot(front, label='NSGAII-BB50011', filename='NSGAII-BB50011')

    # find extreme solutions
    solutions = get_representative_set(front)
    for solution in solutions:
        print(solution.objectives)

    print('Computing time: ' + str(algorithm.total_computing_time))
