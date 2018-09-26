from math import sqrt

from distributed import LocalCluster, Client
from jmetal.algorithm import NSGAII
from jmetal.component import SequentialEvaluator, ProgressBarObserver, RankingAndCrowdingDistanceComparator, \
    VisualizerObserver
from jmetal.core.problem import FloatProblem
from jmetal.core.solution import FloatSolution
from jmetal.operator import BinaryTournamentSelection, Polynomial, SBX
from pymsa.core.score import SumOfPairs, PercentageOfTotallyConservedColumns

from pym2sa.algorithm import NSGA2MSA
from pym2sa.algorithm.multiobjective.dnsgaii import dNSGAII
from pym2sa.component import ParallelEvaluator
from pym2sa.problem import BAliBASE
from pym2sa.operator import SPXMSA, ShiftClosedGapGroups
from pym2sa.util.graphic import MSAPlot

class ZDT1(FloatProblem):
    """ Problem ZDT1.

    .. note:: Bi-objective unconstrained problem. The default number of variables is 30.
    .. note:: Continuous problem having a convex Pareto front
    """

    def __init__(self, number_of_variables: int=30, rf_path: str=None):
        """ :param number_of_variables: Number of decision variables of the problem.
        :param rf_path: Path to the reference front file (if any). Default to None.
        """
        super(ZDT1, self).__init__(rf_path=rf_path)
        self.number_of_variables = number_of_variables
        self.number_of_objectives = 2
        self.number_of_constraints = 0

        self.obj_directions = [self.MINIMIZE, self.MINIMIZE]
        self.obj_labels = ['f(x)', 'f(y)']

        self.lower_bound = self.number_of_variables * [0.0]
        self.upper_bound = self.number_of_variables * [1.0]

    def evaluate(self, solution: FloatSolution) -> FloatSolution:
        g = self.__eval_g(solution)
        h = self.__eval_h(solution.variables[0], g)

        solution.objectives[0] = solution.variables[0]
        solution.objectives[1] = h * g

        sum: float = 0.0
        for i in range(2000000):
            sum += i * 0.235

        return solution

    def __eval_g(self, solution: FloatSolution):
        g = sum(solution.variables) - solution.variables[0]

        constant = 9.0 / (solution.number_of_variables - 1)
        g = constant * g
        g = g + 1.0

        return g

    def __eval_h(self, f: float, g: float) -> float:
        return 1.0 - sqrt(f / g)

    def get_name(self):
        return 'ZDT1'

def setup_distributed_client(address: str):
    new_client = Client(address)

    # This method will send (and import) the module up to all worker nodes in the cluster
    # Note: this file must be created by running `python setup.py install`
    new_client.upload_file('F:\Softw\Python\pyM2SA\dist\pym2sa.egg')

    # Also, each worker should install the dependencies
    def install_dependencies_on_workers():
        import os
        os.system('pip install pymsa jmetalpy')
    new_client.run(install_dependencies_on_workers)

    return new_client

if __name__ == '__main__':
    problem = ZDT1()

    # Setup Dask client
    client = setup_distributed_client('150.214.108.108:8786')

    # Creates the algorithm
    algorithm = dNSGAII(
        problem=problem,
        population_size=100,
        max_evaluations=1000,
        mutation=Polynomial(probability=1.0 / problem.number_of_variables, distribution_index=20),
        crossover=SBX(probability=1.0, distribution_index=20),
        selection=BinaryTournamentSelection(comparator=RankingAndCrowdingDistanceComparator()),
        number_of_cores=8,
        client=client
    )

    visualizer = VisualizerObserver()
    algorithm.observable.register(observer=visualizer)

    algorithm.run()
    front = algorithm.get_result()

    print('Computing time: ' + str(algorithm.total_computing_time))
