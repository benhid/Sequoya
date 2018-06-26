from os.path import dirname, join
from os import listdir
from typing import List
import logging
import random

from pymsa.util.fasta import read_fasta_file_as_list_of_pairs

from pym2sa.core.solution import MSASolution
from pym2sa.operator.crossover import SPXMSA
from pym2sa.operator.mutation import TwoRandomAdjacentGapGroup
from pym2sa.problem.MSA import MSA

logger = logging.getLogger(__name__)

BASE_PATH = dirname(join(dirname(__file__))) + '/problem/aligned'
DATA_FILES = ['tfa_clu', 'tfa_muscle', 'tfa_kalign', 'tfa_retalign',
              'fasta_aln', 'tfa_probcons', 'tfa_mafft', 'tfa_fsa']


class BAliBaseMSA(MSA):

    def __init__(self, instance: str) -> None:
        super(BAliBaseMSA, self).__init__(number_of_variables=0)
        self.problem_name = instance
        self.number_of_objectives = 2
        self.number_of_constraints = 0

    def create_solutions(self, population_size: int) -> List[MSASolution]:
        population = []

        r_v_file = '/RV{0}/'.format(self.problem_name[2:4])
        computed_path = BASE_PATH + r_v_file

        logger.info('Reading path ' + computed_path)

        try:
            for file in listdir(computed_path):
                if file.split('.')[0] == self.problem_name and file.split('.')[1] in DATA_FILES:
                    logger.info('...reading file ' + r_v_file + file)
                    fasta_file = read_fasta_file_as_list_of_pairs(file, computed_path)

                    msa = MSASolution(aligned_sequences=fasta_file, number_of_objectives=2)
                    population.append(msa)
        except FileNotFoundError:
            raise Exception('Instance not found. Invalid path provided? {0}'.format(computed_path))

        if len(population) < 2:
            raise Exception('More than one pre-computed alignment is needed!')

        logger.info('Initial population imported!')
        self.number_of_variables = population[0].number_of_variables

        for index, p in enumerate(population):
            self.evaluate(p)
            logger.info('Alignment {0} size: {1}, Objectives: {2}, {3}'.format(
                index, p.get_length_of_alignment(), p.objectives[0], p.objectives[1])
            )

        logger.info('Number of pre-computed alignments: {0}'.format(len(population)))

        with open('PRECOMPUTED_ALIGNMENTS', 'w') as of:
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
            population.append(offspring[1])
            logger.info('Population incremented by 2; new population {0}'.format(len(population)))

        logger.info('Population incremented to: {0}'.format(len(population)))

        with open('INITIAL_POPULATION', 'w') as of:
            for solution in population:
                of.write(str(solution) + " ")
                of.write("\n")

        return population

    def get_name(self) -> str:
        return "Multiple Sequence Alignment (MSA) BAliBASE problem"
