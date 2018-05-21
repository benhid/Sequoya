from os.path import dirname, join
from os import listdir
from typing import List
import logging
import random

from pymsa.core.score import PercentageOfNonGaps, SumOfPairs, PercentageOfTotallyConservedColumns, Strike, Entropy
from pymsa.util.fasta import read_fasta_file_as_list_of_pairs

from pym2sa.core.solution import MSASolution
from pym2sa.operators.crossover import SPXMSA
from pym2sa.problem.MSA import MSA

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

BASE_PATH = dirname(join(dirname(__file__)))
ALIGNED_SEQUENCES_PATH = '/problem/aligned'
DATA_FILES = ['tfa_clu', 'tfa_muscle', 'tfa_kalign', 'tfa_retalign',
              'fasta_aln', 'tfa_probcons', 'tfa_mafft', 'tfa_fsa']


class BalisebaseMSA(MSA):

    def __init__(self, instance: str, number_of_variables: int) -> None:
        super(BalisebaseMSA, self).__init__(number_of_variables)
        self.problem_name = instance
        self.number_of_objectives = 2
        self.number_of_variables = number_of_variables
        self.number_of_constraints = 0

    def create_solutions(self, population_size: int) -> List[MSASolution]:
        population = []

        r_v_file = '/RV' + self.problem_name[2:4] + '/'
        computed_path = BASE_PATH + ALIGNED_SEQUENCES_PATH + r_v_file
        logger.info('Reading path ' + computed_path)

        try:
            for file in listdir(computed_path):
                if file.split('.')[0] == self.problem_name and file.split('.')[1] in DATA_FILES:
                    logger.info('Reading file ' + r_v_file + file)
                    fasta_file = read_fasta_file_as_list_of_pairs(file, computed_path)

                    msa = MSASolution(aligned_sequences=fasta_file, number_of_objectives=2)
                    population.append(msa)
        except FileNotFoundError:
            raise Exception('Invalid path provided: {0}'.format(computed_path))

        if len(population) < 2:
            raise Exception('More than one pre-computed alignment is needed!')

        logger.info('Initial population imported!')

        for index, p in enumerate(population):
            logger.info('Alignment {0} size: {1}'.format(index, p.original_alignment_size))

        logger.info('Number of pre-computed alignments: {0}'.format(len(population)))

        crossover_operator = SPXMSA(probability=0.8)

        while len(population) < population_size:
            a, b = random.sample(range(len(population)-1), 2)
            offspring = crossover_operator.execute([population[a], population[b]])

            logger.info("Population incremented by 1; new population {0}".format(len(population)))
            population.append(offspring[0])

        logger.info('Population incremented to: {0}'.format(len(population)))

        return population

    def create_solution(self) -> Exception:
        raise NotImplemented

    def get_name(self) -> str:
        return "Multiple Sequence Alignment (MSA) BaliBASE problem"
