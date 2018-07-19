from typing import List

from pymsa.core.score import Score

from pym2sa.core.problem import MSAProblem
from pym2sa.core.solution import MSASolution


class MSA(MSAProblem):

    def __init__(self, score_list: List[Score]) -> None:
        """ Creates a new MSA problem.

        :param score_list: """
        super(MSA, self).__init__()
        self.score_list = score_list
        self.number_of_objectives = len(self.score_list)

        self.original_sequences: list = []
        self.sequences_names: list = []

    def create_solution(self) -> List[MSASolution]:
        raise NotImplementedError()

    def evaluate(self, solution: MSASolution):
        solution.remove_full_of_gaps_columns()

        for i, score in enumerate(self.score_list):
            if score.is_minimization():
                solution.objectives[i] = score.compute(solution.decode_alignment_as_list_of_sequences())
            else:
                solution.objectives[i] = -1.0 * score.compute(solution.decode_alignment_as_list_of_sequences())

    def get_name(self) -> str:
        return 'Multiple Sequence Alignment problem'
