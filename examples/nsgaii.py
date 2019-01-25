from jmetal.algorithm.multiobjective.nsgaii import NSGAII
from pymsa.core.score import SumOfPairs, PercentageOfTotallyConservedColumns

from jmetal.operator import BinaryTournamentSelection
from jmetal.util.comparator import RankingAndCrowdingDistanceComparator
from jmetal.util.observer import VisualizerObserver, ProgressBarObserver
from jmetal.util.solution_list.evaluator import MultiprocessEvaluator
from jmetal.util.termination_criterion import StoppingByEvaluations
from jmetal.util.visualization import Plot
from pym2sa.operator import SPXMSA, ShiftClosedGapGroups
from pym2sa.problem import BAliBASE
from pym2sa.util.visualization import MSAPlot

if __name__ == '__main__':
    # Creates the problem
    problem = BAliBASE(balibase_instance='BB11030', balibase_path='../resources',
                       score_list=[SumOfPairs(), PercentageOfTotallyConservedColumns()])

    # Creates the algorithm
    max_evaluations = 10000

    algorithm = NSGAII(
        problem=problem,
        population_size=100,
        offspring_population_size=100,
        mutation=ShiftClosedGapGroups(probability=0.2),
        crossover=SPXMSA(probability=0.0),
        selection=BinaryTournamentSelection(comparator=RankingAndCrowdingDistanceComparator()),
        termination_criterion=StoppingByEvaluations(max=max_evaluations),
        population_evaluator=MultiprocessEvaluator()
    )

    algorithm.observable.register(observer=VisualizerObserver())
    algorithm.observable.register(ProgressBarObserver(max=max_evaluations))

    algorithm.run()
    front = algorithm.get_result()

    # Plot front
    plot_front = Plot(plot_title='Pareto front approximation', axis_labels=['%SOP', '%TC'])
    plot_front.plot(front, label='NSGAIII-BB12010', filename='NSGAIII-BB12010')

    # Plot interactive front
    pareto_front = MSAPlot(plot_title='Pareto front approximation', axis_labels=['%SOP', '%TC'])
    pareto_front.plot(front, label='NSGAIII-BB12010', filename='NSGAIII-BB12010')

    print('Computing time: ' + str(algorithm.total_computing_time))
