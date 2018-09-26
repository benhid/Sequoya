import os
from os import listdir
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

    def create_solution(self) -> List[MSASolution]:
        raise NotImplementedError()

    def import_instance(self, population_size: int) -> List[MSASolution]:
        """ Read and import an instance of BAliBASE.

        :param population_size: If the instance has less pre-computed alignments than this value, the population will be increased until the population_size is met. """
        self.__read_original_sequences()
        aligned_sequences = self.__read_aligned_sequences()

        if len(aligned_sequences) < 2:
            raise Exception('More than one pre-computed alignment is needed!')

        population = []
        for msa in aligned_sequences:
            new_individual = MSASolution(self, msa)
            population.append(new_individual)

        logger.info('Instance imported')

        for index, individual in enumerate(population):
            self.evaluate(individual)
            logger.info('Alignment no. {0}, length of alignment: {1}, objective values: {2}'.format(
                index, individual.get_length_of_alignment(), individual.objectives)
            )

        logger.info('Number of pre-computed alignments: {0}'.format(len(population)))

        with open('pyM2SA-PRECOMPUTED_ALIGNMENTS_{}'.format(self.instance), 'w') as of:
            for solution in population:
                of.write(str(solution) + " ")
                of.write("\n")

        crossover_operator = SPXMSA(probability=1.0)
        mutation_operator = TwoRandomAdjacentGapGroup(probability=1.0)

        while len(population) < population_size:
            a = random.randint(0, len(population)-1)
            b = random.randint(0, len(population)-1)

            while a == b:
                b = random.randint(0, len(population) - 1)

            offspring = crossover_operator.execute([population[a], population[b]])
            mutation_operator.execute(offspring[0])
            mutation_operator.execute(offspring[1])

            population.append(offspring[0])
            logger.info('Population incremented, new population size: {}'.format(len(population)))

        logger.info('Final population size: {}'.format(len(population)))

        with open('pyM2SA-INITIAL_POPULATION_{}'.format(self.instance), 'w') as of:
            for solution in population:
                of.write(str(solution) + " ")
                of.write("\n")

        return population

    def __read_original_sequences(self):
        bb3_release_path = self.__compute_path('bb_release')

        if os.path.isdir(bb3_release_path):
            fasta_file = read_fasta_file_as_list_of_pairs(self.instance + '.tfa', bb3_release_path)

            self.sequences_names = list(pair[0] for pair in fasta_file)
            self.number_of_variables = len(self.sequences_names)
            self.original_sequences = fasta_file
        else:
            raise Exception('Instance not found at {}'.format(bb3_release_path))

    def __read_aligned_sequences(self):
        bb3_aligned_path = self.__compute_path('bb_aligned')

        alignment_sequences = []
        if os.path.isdir(bb3_aligned_path):
            for file in listdir(bb3_aligned_path):
                if file.split('.')[0] == self.instance and file.split('.')[1] in self.DATA_FILES:
                    msa = read_fasta_file_as_list_of_pairs(file, bb3_aligned_path)
                    alignment_sequences.append(msa)
        else:
            raise Exception('Instance not found at {}'.format(bb3_aligned_path))

        return alignment_sequences

    def __compute_path(self, directory: str) -> str:
        return os.path.join(self.balibase_path, directory, 'RV' + self.instance[2:4] + '/')

    def get_name(self) -> str:
        return 'BAliBASE problem'
