from typing import List

from pymsa.core.score import Score

from pym2sa.core.problem import MSAProblem
from pym2sa.core.solution import MSASolution


class MSA(MSAProblem):

    def __init__(self, number_of_variables: int, score_list: List[Score]) -> None:
        super(MSA, self).__init__()
        self.score_list = score_list
        self.number_of_variables = number_of_variables
        self.number_of_objectives = len(self.score_list)
        self.number_of_constraints = 0

    def evaluate(self, solution: MSASolution):
        solution.remove_full_of_gaps_columns()

        for i, score in enumerate(self.score_list):
            if score.is_minimization():
                solution.objectives[i] = score.compute(solution.decode_alignment_as_list_of_sequences())
            else:
                solution.objectives[i] = -1.0 * score.compute(solution.decode_alignment_as_list_of_sequences())

    def get_name(self) -> str:
        return "Multiple Sequence Alignment problem"
