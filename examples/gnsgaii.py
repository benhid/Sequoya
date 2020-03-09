from jmetal.algorithm.multiobjective.nsgaii import NSGAII
from jmetal.lab.visualization import Plot
from jmetal.util.observer import ProgressBarObserver, PlotFrontToFileObserver
from jmetal.util.solutions.comparator import GDominanceComparator
from jmetal.util.termination_criterion import StoppingByEvaluations
from pymsa.core.score import SumOfPairs, PercentageOfTotallyConservedColumns

from sequoya.operator import SPXMSA, ShiftClosedGapGroups
from sequoya.problem import BAliBASE
from sequoya.util.visualization import MSAPlot

if __name__ == '__main__':
    # creates the problem
    problem = BAliBASE(instance='BB20019', path='../resources',
                       score_list=[SumOfPairs(), PercentageOfTotallyConservedColumns()])

    # creates the algorithm
    max_evaluations = 25000
    reference_point = [-0.6, 11000]

    algorithm = NSGAII(
        problem=problem,
        population_size=100,
        offspring_population_size=100,
        mutation=ShiftClosedGapGroups(probability=0.3),
        crossover=SPXMSA(probability=0.7),
        dominance_comparator=GDominanceComparator(reference_point),
        termination_criterion=StoppingByEvaluations(max_evaluations=max_evaluations)
    )

    algorithm.observable.register(observer=ProgressBarObserver(max=max_evaluations))
    algorithm.observable.register(observer=PlotFrontToFileObserver(output_directory='fronts_zdt1'))

    algorithm.run()
    front = algorithm.get_result()

    # plot front
    plot_front = Plot(title='Pareto front approximation', axis_labels=['%SOP', '%TC'])
    plot_front.plot(front, label='NSGAII-BB20019', filename='NSGAII-BB20019')

    # plot interactive front
    pareto_front = MSAPlot(title='Pareto front approximation', axis_labels=['%SOP', '%TC'])
    pareto_front.plot(front, label='NSGAII-BB20019', filename='NSGAII-BB20019')

    print('Computing time: ' + str(algorithm.total_computing_time))
