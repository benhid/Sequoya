import os
from os import listdir
from typing import List
import logging
import random

from pymsa.core.score import Score
from pymsa.util.fasta import read_fasta_file_as_list_of_pairs

from sequoya.core.solution import MSASolution
from sequoya.operator import SPXMSA, TwoRandomAdjacentGapGroup
from sequoya.problem.MSA import MSA

LOGGER = logging.getLogger('Sequoya')


class BAliBASE(MSA):

    DATA_FILES = ['tfa_clu', 'tfa_muscle', 'tfa_kalign', 'tfa_retalign',
                  'fasta_aln', 'tfa_probcons', 'tfa_mafft', 'tfa_fsa']

    def __init__(self, balibase_instance: str, balibase_path: str, score_list: List[Score]) -> None:
        """ Creates a new problem based on an instance of BAliBASE.

        :param balibase_instance: Instance name (e.g., BB12010).
        :param balibase_path: Path containing two directories: `bb3_aligned`, with the pre-computed alignments and
        `bb3_release`, with the original sequences.
        :param score_list: List of pyMSA objects. """
        super(BAliBASE, self).__init__(score_list)
        self.balibase_instance = balibase_instance
        self.balibase_path = balibase_path
        self.instance = self.import_instance()
        self.counter = 0

    def create_solution(self) -> MSASolution:
        if self.counter < len(self.instance):
            offspring = self.instance[self.counter]
            self.counter += 1
        else:
            crossover_operator = SPXMSA(probability=1.0)
            mutation_operator = TwoRandomAdjacentGapGroup(probability=1.0)

            a = random.randint(0, len(self.instance) - 1)
            b = random.randint(0, len(self.instance) - 1)

            while a == b:
                b = random.randint(0, len(self.instance) - 1)

            offspring = crossover_operator.execute([self.instance[a], self.instance[b]])
            mutation_operator.execute(offspring[0])

            offspring = offspring[0]

        return offspring

    def import_instance(self) -> List[MSASolution]:
        self.read_original_sequences()

        bb3_aligned_path = self.compute_path('bb3_aligned')

        alignment_sequences = []
        if os.path.isdir(bb3_aligned_path):
            for file in listdir(bb3_aligned_path):
                if file.split('.')[0] == self.balibase_instance and file.split('.')[1] in self.DATA_FILES:
                    msa = read_fasta_file_as_list_of_pairs(bb3_aligned_path + '/' + file)
                    alignment_sequences.append(msa)
        else:
            raise Exception('Instance not found. Invalid path provided? {}'.format(bb3_aligned_path))

        if len(alignment_sequences) < 2:
            raise Exception('More than one pre-computed alignment is needed!')

        population = []

        for msa in alignment_sequences:
            new_individual = MSASolution(self, msa)
            population.append(new_individual)

        LOGGER.info('Instance imported')

        for index, individual in enumerate(population):
            self.evaluate(individual)
            LOGGER.info('Alignment {0} size: {1}, Objectives: {2}'.format(
                index, individual.get_length_of_alignment(), individual.objectives)
            )

        LOGGER.info('Number of pre-computed alignments: {0}'.format(len(population)))

        with open('PRECOMPUTED_ALIGNMENTS_{0}'.format(self.balibase_instance), 'w') as of:
            for solution in population:
                of.write(str(solution) + " ")
                of.write("\n")

        return population

    def read_original_sequences(self):
        bb3_release_path = self.compute_path('bb3_release')

        if os.path.isdir(bb3_release_path):
            fasta_file = read_fasta_file_as_list_of_pairs(bb3_release_path + '/' + self.balibase_instance + '.tfa', )
            self.sequences_names = list(pair[0] for pair in fasta_file)
            self.number_of_variables = len(self.sequences_names)
            self.original_sequences = fasta_file
        else:
            raise Exception('Instance not found. Invalid path provided? {}'.format(bb3_release_path))

    def compute_path(self, directory: str) -> str:
        return os.path.join(self.balibase_path, directory, 'RV' + self.balibase_instance[2:4] + '/')

    def get_name(self) -> str:
        return 'BAliBASE problem'
