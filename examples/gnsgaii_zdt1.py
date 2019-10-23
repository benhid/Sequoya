from jmetal.algorithm.multiobjective.nsgaii import NSGAII
from jmetal.lab.visualization import Plot
from jmetal.operator import PolynomialMutation, SBXCrossover
from jmetal.problem import ZDT1
from jmetal.util.observer import ProgressBarObserver, VisualizerObserver, PlotFrontToFileObserver
from jmetal.util.solutions.comparator import GDominanceComparator
from jmetal.util.termination_criterion import StoppingByEvaluations
from pymsa.core.score import SumOfPairs, PercentageOfTotallyConservedColumns

from sequoya.operator import SPXMSA, ShiftClosedGapGroups
from sequoya.problem import BAliBASE
from sequoya.util.visualization import MSAPlot

if __name__ == '__main__':
    # creates the problem
    problem = ZDT1()

    # creates the algorithm
    max_evaluations = 25000
    reference_point = [0.2, 0.6]

    algorithm = NSGAII(
        problem=problem,
        population_size=100,
        offspring_population_size=100,
        mutation=PolynomialMutation(probability=1.0 / problem.number_of_variables, distribution_index=20),
        crossover=SBXCrossover(probability=1.0, distribution_index=20),
        dominance_comparator=GDominanceComparator(reference_point),
        termination_criterion=StoppingByEvaluations(max=max_evaluations)
    )

    algorithm.observable.register(observer=ProgressBarObserver(max=max_evaluations))
    algorithm.observable.register(observer=PlotFrontToFileObserver(reference_point=reference_point,
                                                                   output_directory='fronts_zdt1'))
    #algorithm.observable.register(observer=VisualizerObserver(reference_point=reference_point))

    algorithm.run()
    front = algorithm.get_result()

    # plot front
    plot_front = Plot(plot_title='Pareto front approximation', axis_labels=['%SOP', '%TC'])
    plot_front.plot(front, label='NSGAII-BB20019', filename='NSGAII-BB20019')

    # plot interactive front
    pareto_front = MSAPlot(plot_title='Pareto front approximation', axis_labels=['%SOP', '%TC'])
    pareto_front.plot(front, label='NSGAII-BB20019', filename='NSGAII-BB20019')

    print('Computing time: ' + str(algorithm.total_computing_time))
