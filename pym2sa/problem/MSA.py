from pymsa.core.score import PercentageOfNonGaps, SumOfPairs, PercentageOfTotallyConservedColumns, Strike, Entropy
from jmetal.core.objective import Objective

from pym2sa.core.problem import MSAProblem
from pym2sa.core.solution import MSASolution


class MSA(MSAProblem):
    def __init__(self, number_of_variables: int) -> None:
        super(MSA, self).__init__()
        self.objectives = [self.Objective1(), self.Objective2()]
        self.number_of_objectives = len(self.objectives)
        self.number_of_variables = number_of_variables
        self.number_of_constraints = 0

    class Objective1(Objective):
        def compute(self, solution: MSASolution, problem: MSAProblem):
            return Entropy().compute(solution.decode_alignment())

        def is_a_minimization_objective(self) -> bool:
            return Entropy().is_minimization()

    class Objective2(Objective):
        def compute(self, solution: MSASolution, problem: MSAProblem):
            return SumOfPairs().compute(solution.decode_alignment())

        def is_a_minimization_objective(self) -> bool:
            return SumOfPairs().is_minimization()

    def create_solution(self) -> None:
        raise Exception("Not able to create any solution to MSA!")

    def get_name(self) -> str:
        return "Multiple Sequence Alignment (MSA) problem"
