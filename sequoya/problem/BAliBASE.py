import logging
import os
import random
from os import listdir
from typing import List

from pymsa.core.score import Score
from pymsa.util.fasta import read_fasta_file_as_list_of_pairs

from sequoya.core.solution import MSASolution
from sequoya.operator import SPXMSA, TwoRandomAdjacentGapGroup
from sequoya.problem.MSA import MSA

LOGGER = logging.getLogger('Sequoya')


class BAliBASE(MSA):
    DATA_FILES = ['tfa_clu', 'tfa_muscle', 'tfa_kalign', 'tfa_retalign', 'fasta_aln', 'tfa_probcons', 'tfa_mafft',
                  'tfa_fsa']

    def __init__(self, instance: str, path: str, score_list: List[Score], auto_import: bool = True) -> None:
        """
        Creates a new problem based on an instance of BAliBASE.

        :param instance: Instance name (e.g., BB12010).
        :param path: Path containing two directories: `bb3_aligned`, with the pre-computed alignments and
        `bb3_release`, with the original sequences.
        :param score_list: List of score functions.
        """
        super(BAliBASE, self).__init__(score_list)
        self.instance = instance
        self.path = path

        if auto_import:
            self.import_instance()

    def create_solution(self) -> MSASolution:
        crossover_operator = SPXMSA(probability=1.0)
        mutation_operator = TwoRandomAdjacentGapGroup(probability=1.0)

        a = random.randint(0, len(self.sequences) - 1)
        b = random.randint(0, len(self.sequences) - 1)

        while a == b:
            b = random.randint(0, len(self.sequences) - 1)

        offspring = crossover_operator.execute([self.sequences[a], self.sequences[b]])
        mutation_operator.execute(offspring[0])

        return offspring[0]

    def import_instance(self) -> List[MSASolution]:
        bb3_release_path = self._compute_path('bb3_release')
        assert os.path.isdir(bb3_release_path), 'Instance not found'

        msa = read_fasta_file_as_list_of_pairs(f'{bb3_release_path}/{self.instance}.tfa')
        self.identifiers = list(pair[0] for pair in msa)
        self.number_of_variables = len(self.identifiers)

        bb3_aligned_path = self._compute_path('bb3_aligned')
        assert os.path.isdir(bb3_aligned_path), 'Instance not found'

        multiple_alignments = []
        for file in listdir(bb3_aligned_path):
            name, fmt = file.split('.')

            if name == self.instance and fmt in self.DATA_FILES:
                msa = read_fasta_file_as_list_of_pairs(f'{bb3_aligned_path}/{file}')
                multiple_alignments.append(msa)

        if len(multiple_alignments) < 2:
            raise Exception('More than one pre-computed MSA is required')

        population = []
        for msa in multiple_alignments:
            new_individual = MSASolution(self, msa)
            population.append(new_individual)

        LOGGER.info('Instance imported')

        for index, individual in enumerate(population):
            self.evaluate(individual)
            LOGGER.info(f'ALN {index}, LEN: {individual.get_length_of_alignment()}, OBJ: {individual.objectives}')

        self.sequences = population

        return population

    def _compute_path(self, directory: str) -> str:
        return os.path.join(self.path, directory, 'RV' + self.instance[2:4] + '/')

    def get_name(self) -> str:
        return 'BAliBASE v3.0'
