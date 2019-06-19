from jmetal.algorithm.multiobjective.nsgaii import NSGAII
from jmetal.operator import BinaryTournamentSelection
from jmetal.util.comparator import RankingAndCrowdingDistanceComparator
from jmetal.util.observer import ProgressBarObserver, VisualizerObserver
from jmetal.util.solution_list.evaluator import SequentialEvaluator
from jmetal.util.termination_criterion import StoppingByEvaluations
from jmetal.util.visualization import Plot
from pymsa.core.score import SumOfPairs, PercentageOfTotallyConservedColumns

from sequoya.operator import SPXMSA, ShiftClosedGapGroups
from sequoya.problem import BAliBASE
from sequoya.util.visualization import MSAPlot

if __name__ == '__main__':
    # creates the problem
    problem = BAliBASE(instance='BB50011', path='../resources',
                       score_list=[SumOfPairs(), PercentageOfTotallyConservedColumns()])

    # creates the algorithm
    max_evaluations = 25000

    algorithm = NSGAII(
        problem=problem,
        population_size=100,
        offspring_population_size=100,
        mutation=ShiftClosedGapGroups(probability=0.3),
        crossover=SPXMSA(probability=0.7),
        selection=BinaryTournamentSelection(comparator=RankingAndCrowdingDistanceComparator()),
        termination_criterion=StoppingByEvaluations(max=max_evaluations),
        population_evaluator=SequentialEvaluator()
    )

    algorithm.observable.register(observer=ProgressBarObserver(max=max_evaluations))
    algorithm.observable.register(observer=VisualizerObserver())

    algorithm.run()
    front = algorithm.get_result()

    # plot front
    plot_front = Plot(plot_title='Pareto front approximation', axis_labels=['%SOP', '%TC'])
    plot_front.plot(front, label='NSGAIII-BB50011', filename='NSGAIII-BB50011')

    # plot interactive front
    pareto_front = MSAPlot(plot_title='Pareto front approximation', axis_labels=['%SOP', '%TC'])
    pareto_front.plot(front, label='NSGAIII-BB50011', filename='NSGAIII-BB50011')

    print('Computing time: ' + str(algorithm.total_computing_time))
