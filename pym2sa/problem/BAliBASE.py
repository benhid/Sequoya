import os
from typing import List
import logging
import random

from pymsa.core.score import Score
from pymsa.util.fasta import read_fasta_file_as_list_of_pairs

from pym2sa.core.solution import MSASolution
from pym2sa.operator import SPXMSA, TwoRandomAdjacentGapGroup
from pym2sa.problem.MSA import MSA

logger = logging.getLogger('pyM2SA')


class BAliBASE(MSA):

    DATA_FILES = ['tfa_clu', 'tfa_muscle', 'tfa_kalign', 'tfa_retalign',
                  'fasta_aln', 'tfa_probcons', 'tfa_mafft', 'tfa_fsa']

    def __init__(self, instance: str, balibase_path: str, score_list: List[Score]) -> None:
        """ Creates a new problem based on an instance of BAliBASE.

        :param instance: Instance name (e.g., BB12010).
        :param balibase_path: Path containing two directories: `bb_aligned`, with the pre-computed alignments and `bb_release`, with the original sequences.
        :param score_list: List of scores. """
        super(BAliBASE, self).__init__(score_list, [], [])
        self.balibase_path = balibase_path
        self.instance = instance

        self.instance_population = []
        self.instance_index = 0

        self.__import_instance()

    def create_solution(self) -> MSASolution:
        """ Read and import an instance of BAliBASE. """
        if self.instance_index < len(self.sequences_names):
            solution = self.instance_population[self.instance_index]
            self.instance_index += 1
        else:
            crossover_operator = SPXMSA(probability=1.0)
            mutation_operator = TwoRandomAdjacentGapGroup(probability=1.0)

            a = random.randint(0, len(self.sequences_names) - 1)
            b = random.randint(0, len(self.sequences_names) - 1)

            while a == b:
                b = random.randint(0, len(self.sequences_names) - 1)

            # We are only interested on one offspring
            solution = crossover_operator.execute([self.instance_population[a], self.instance_population[b]])[0]
            mutation_operator.execute(solution)

        return solution

    def __import_instance(self):
        original_sequences = self.__read_original_sequences()

        self.sequences_names = list(pair[0] for pair in original_sequences)
        self.original_sequences = list(pair[1] for pair in original_sequences)
        self.number_of_variables = len(self.sequences_names)

        aligned_sequences = self.__read_aligned_sequences()

        if len(aligned_sequences) < 2:
            raise Exception('More than one pre-computed alignment is required!')

        for msa in aligned_sequences:
            new_individual = MSASolution(self, msa)

            self.evaluate(new_individual)
            logger.info('Length of alignment: {0}, objective values: {1}'.format(
                new_individual.get_length_of_alignment(), new_individual.objectives)
            )

            self.instance_population.append(new_individual)

    def __read_original_sequences(self):
        bb3_release_path = self.__compute_path('bb_release')

        if os.path.isdir(bb3_release_path):
            fasta_file = read_fasta_file_as_list_of_pairs(self.instance + '.tfa', bb3_release_path)
        else:
            raise Exception('Instance not found at {}'.format(bb3_release_path))

        return fasta_file

    def __read_aligned_sequences(self) -> List[list]:
        bb3_aligned_path = self.__compute_path('bb_aligned')

        aligned_sequences = []
        if os.path.isdir(bb3_aligned_path):
            for file in os.listdir(bb3_aligned_path):
                if file.split('.')[0] == self.instance and file.split('.')[1] in self.DATA_FILES:
                    msa = read_fasta_file_as_list_of_pairs(file, bb3_aligned_path)
                    aligned_sequences.append(msa)
        else:
            raise Exception('Instance not found at {}'.format(bb3_aligned_path))

        return aligned_sequences

    def __compute_path(self, directory: str) -> str:
        return os.path.join(self.balibase_path, directory, 'RV' + self.instance[2:4] + '/')

    def get_name(self) -> str:
        return 'BAliBASE problem'
