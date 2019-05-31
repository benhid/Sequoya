from math import sqrt

import matplotlib
from jmetal.problem import ZDT1

matplotlib.use('TkAgg')

from jmetal.algorithm.multiobjective.nsgaii import NSGAII
from jmetal.operator import BinaryTournamentSelection, PolynomialMutation, SBXCrossover
from jmetal.util.comparator import RankingAndCrowdingDistanceComparator
from jmetal.util.observer import ProgressBarObserver, VisualizerObserver
from jmetal.util.solution_list.evaluator import SequentialEvaluator
from jmetal.util.termination_criterion import StoppingByEvaluations
from jmetal.util.visualization import Plot

from sequoya.util.visualization import MSAPlot


class ZDT1Modified(ZDT1):
    """
    Problem ZDT1.

    .. note:: Version including a loop for increasing the computing time of the evaluation functions.
    """

    def evaluate(self, solution):
        g = self.__eval_g(solution)
        h = self.__eval_h(solution.variables[0], g)

        solution.objectives[0] = solution.variables[0]
        solution.objectives[1] = h * g

        s: float = 0.0
        for i in range(1000000000):
            s += i * 0.235 / 1.234

        return solution

    def __eval_g(self, solution):
        g = sum(solution.variables) - solution.variables[0]

        constant = 9.0 / (solution.number_of_variables - 1)
        g = constant * g
        g = g + 1.0

        return g

    def __eval_h(self, f: float, g: float) -> float:
        return 1.0 - sqrt(f / g)

    def get_name(self):
        return 'ZDT1m'


if __name__ == '__main__':
    # creates the problem
    problem = ZDT1Modified()

    # creates the algorithm
    max_evaluations = 1000

    algorithm = NSGAII(
        problem=problem,
        population_size=100,
        offspring_population_size=100,
        mutation=PolynomialMutation(probability=1.0 / problem.number_of_variables, distribution_index=20),
        crossover=SBXCrossover(probability=1.0, distribution_index=20),
        selection=BinaryTournamentSelection(comparator=RankingAndCrowdingDistanceComparator()),
        termination_criterion=StoppingByEvaluations(max=max_evaluations),
        population_evaluator=SequentialEvaluator()
    )

    algorithm.observable.register(observer=ProgressBarObserver(max=max_evaluations))
    algorithm.observable.register(observer=VisualizerObserver())

    algorithm.run()
    front = algorithm.get_result()

    # plot front
    plot_front = Plot(plot_title='Pareto front approximation')
    plot_front.plot(front, label='ZDT1', filename='ZDT1')

    # plot interactive front
    pareto_front = MSAPlot(plot_title='Pareto front approximation')
    pareto_front.plot(front, label='ZDT1', filename='ZDT1')

    print('Computing time: ' + str(algorithm.total_computing_time))
