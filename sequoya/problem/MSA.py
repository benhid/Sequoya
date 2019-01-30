from typing import List

from pymsa.core.score import Score

from sequoya.core.problem import MSAProblem
from sequoya.core.solution import MSASolution


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

    def evaluate(self, solution: MSASolution) -> MSASolution:
        solution.remove_full_of_gaps_columns()
        sequences = solution.decode_alignment_as_list_of_sequences()

        for i, score in enumerate(self.score_list):
            solution.objectives[i] = score.compute(sequences)

            if not score.is_minimization():
                # If score is maximization, invert sign
                solution.objectives[i] = -solution.objectives[i]

        return solution

    def get_name(self) -> str:
        return 'Multiple Sequence Alignment problem'
