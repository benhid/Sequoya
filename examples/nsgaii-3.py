from jmetal.algorithm import NSGAII
from jmetal.component import SequentialEvaluator, ProgressBarObserver, RankingAndCrowdingDistanceComparator, \
    VisualizerObserver
from jmetal.operator import BinaryTournamentSelection
from pymsa.core.score import SumOfPairs, PercentageOfNonGaps, PercentageOfTotallyConservedColumns

from pym2sa.problem import BAliBASE
from pym2sa.operator import SPXMSA, TwoRandomAdjacentGapGroup
from pym2sa.util.graphic import MSAPlot


if __name__ == '__main__':
    # Creates the problem
    problem = BAliBASE(balibase_instance='BB12010', balibase_path='../resources',
                       score_list=[SumOfPairs(), PercentageOfTotallyConservedColumns(), PercentageOfNonGaps()])
    problem.obj_labels = ['SOP', '%TC', '%NonGaps']

    # Creates the algorithm
    algorithm = NSGAII(
        problem=problem,
        population_size=100,
        max_evaluations=10000,
        mutation=TwoRandomAdjacentGapGroup(probability=0.2),
        crossover=SPXMSA(probability=0.8),
        selection=BinaryTournamentSelection(comparator=RankingAndCrowdingDistanceComparator()),
        evaluator=SequentialEvaluator()
    )

    visualizer = VisualizerObserver()
    algorithm.observable.register(observer=visualizer)

    progress_bar = ProgressBarObserver(step=100, maximum=10000)
    algorithm.observable.register(progress_bar)

    algorithm.run()
    front = algorithm.get_result()

    # Plot the solution
    pareto_front = MSAPlot(plot_title='NSGAII for ' + problem.balibase_instance, axis_labels=problem.obj_labels)
    pareto_front.plot(front)
    pareto_front.to_html(filename='NSGAII-' + problem.balibase_instance)

    print('Computing time: ' + str(algorithm.total_computing_time))
