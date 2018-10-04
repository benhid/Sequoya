from typing import List
import logging

from pym2sa.core.problem import MSAProblem
from pymsa.core.score import Score

from pym2sa.core.solution import MSASolution

logger = logging.getLogger('pyM2SA')


class MSA(MSAProblem):

    def __init__(self, score_list: List[Score], sequences_without_gaps: List[str], sequences_names: List[str]):
        """ Creates a new generic MSA problem.

        :param score_list: List of scores to evaluate MSAs.
        :param sequences_without_gaps: List of original sequences (without gaps).
        :param sequences_names: List of sequences names. """
        super(MSA, self).__init__()
        self.score_list = score_list
        self.sequences_without_gaps = sequences_without_gaps
        self.sequences_names = sequences_names

        self.number_of_objectives = len(self.score_list)
        self.number_of_variables = len(self.sequences_names)

        self.obj_directions = [0 for _ in range(len(self.score_list))]

        for i, score in enumerate(score_list):
            if score.is_minimization():
                self.obj_directions[i] = self.MINIMIZE
            else:
                self.obj_directions[i] = self.MAXIMIZE

    def create_solution(self) -> MSASolution:
        raise NotImplementedError()

    def evaluate(self, solution: MSASolution) -> MSASolution:
        """ Evaluate a multiple sequence alignment solution.

        :param solution: MSA to evaluate. """
        solution.remove_full_of_gaps_columns()

        for i, direction in enumerate(self.obj_directions):
            if direction == self.MAXIMIZE:
                solution.objectives[i] = self.score_list[i].compute(solution.decode_alignment_as_list_of_sequences())
            else:
                solution.objectives[i] = -1.0 * self.score_list[i].compute(solution.decode_alignment_as_list_of_sequences())

        return solution

    def get_name(self) -> str:
        return 'Multiple Sequence Alignment problem'
