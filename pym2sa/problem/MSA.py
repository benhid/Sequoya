from pymsa.core.score import PercentageOfNonGaps, SumOfPairs, PercentageOfTotallyConservedColumns, Strike, Entropy
from pymsa.core.substitution_matrix import PAM250

from pym2sa.core.problem import MSAProblem
from pym2sa.core.solution import MSASolution


class MSA(MSAProblem):

    def __init__(self, number_of_variables: int) -> None:
        super(MSA, self).__init__()
        self.number_of_variables = number_of_variables
        self.number_of_objectives = 2
        self.number_of_constraints = 0

    def evaluate(self, solution: MSASolution):
        solution.remove_full_of_gaps_columns()

        solution.objectives[0] = -1.0 * SumOfPairs(PAM250()).compute(solution.decode_alignment())
        solution.objectives[1] = -1.0 * PercentageOfTotallyConservedColumns().compute(solution.decode_alignment())

    def create_solution(self) -> None:
        raise Exception("Not able to create any solution to MSA!")

    def get_name(self) -> str:
        return "Multiple Sequence Alignment problem"
