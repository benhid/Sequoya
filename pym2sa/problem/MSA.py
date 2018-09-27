from typing import List
import logging

from pym2sa.core.problem import MSAProblem
from pymsa.core.score import Score

from pym2sa.core.solution import MSASolution

logger = logging.getLogger('pyM2SA')


class MSA(MSAProblem):

    def __init__(self, score_list: List[Score], original_sequences: [], sequences_names: []):
        """ Creates a new generic MSA problem.

        :param score_list: List of scores to evaluate MSAs.
        :param original_sequences: List of original sequences (without gaps).
        :param sequences_names: List of sequences names. """
        super(MSA, self).__init__()
        self.score_list = score_list
        self.original_sequences = original_sequences
        self.sequences_names = sequences_names

        self.number_of_constraints = 0
        self.number_of_objectives = len(self.score_list)
        self.number_of_variables = len(self.sequences_names)

    def create_solution(self) -> MSASolution:
        raise NotImplementedError()

    def evaluate(self, solution: MSASolution) -> MSASolution:
        """ Evaluate a multiple sequence alignment solution.

        :param solution: MSA to evaluate. """
        solution.remove_full_of_gaps_columns()

        for i, score in enumerate(self.score_list):
            if score.is_minimization():
                solution.objectives[i] = score.compute(solution.decode_alignment_as_list_of_sequences())
            else:
                solution.objectives[i] = -1.0 * score.compute(solution.decode_alignment_as_list_of_sequences())

        return solution

    def get_name(self) -> str:
        return 'Multiple Sequence Alignment problem'
