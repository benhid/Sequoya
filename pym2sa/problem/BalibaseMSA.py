import logging
import os
import random
from typing import List

from pymsa.core.score import PercentageOfNonGaps, SumOfPairs, PercentageOfTotallyConservedColumns, Strike, Entropy
from pymsa.util.fasta import read_fasta_file_as_list_of_pairs

from pym2sa.core.problem import MSAProblem
from pym2sa.core.solution import MSASolution
from pym2sa.operators.crossover import GapSequenceSolutionSinglePoint

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

BALIBASE_PATH = os.path.dirname(__file__)+'/bb3_release/'
ALIGNED_PATH = os.path.dirname(__file__)+'/aligned/'
DATA_FILES = ['tfa_clu', 'tfa_muscle', 'tfa_kalign', 'tfa_retalign',
              'fasta_aln', 'tfa_probcons', 'tfa_mafft', 'tfa_fsa']


class BalisebaseMSA(MSAProblem):

    def __init__(self, instance: str, number_of_variables: int) -> None:
        super(BalisebaseMSA, self).__init__()
        self.problem_name = instance
        self.number_of_objectives = 2
        self.number_of_variables = number_of_variables
        self.number_of_constraints = 0

    def evaluate(self, solution: MSASolution):
        solution.objectives[0] = PercentageOfTotallyConservedColumns().compute(solution.decode_alignment())
        solution.objectives[1] = -1.0 * SumOfPairs().compute(solution.decode_alignment())

    def create_solutions(self, population_size: int) -> List[MSASolution]:
        population = []

        r_v_path = 'RV' + self.problem_name[2:4]

        try:
            for file in os.listdir(ALIGNED_PATH+r_v_path):
                if file.split('.')[0] == self.problem_name and file.split('.')[1] in DATA_FILES:
                    logger.info('Reading file ' + file)

                    fasta_file = read_fasta_file_as_list_of_pairs(file, ALIGNED_PATH+r_v_path+'/')

                    msa = MSASolution(aligned_sequences=fasta_file, number_of_objectives=2)
                    population.append(msa)
        except FileNotFoundError:
            raise Exception('Invalid path provided: {0}'.format(ALIGNED_PATH+self.problem_name))

        if len(population) < 2:
            raise Exception('More than one pre-computed alignment is needed!')

        logger.info('Initial population imported!')

        for index, p in enumerate(population):
            logger.info('Alignment {0} size: {1}'.format(index, p.original_alignment_size))

        logger.info('Number of pre-computed alignments: {0}'.format(len(population)))

        crossover_operator = GapSequenceSolutionSinglePoint(probability=0.8)

        while len(population) < population_size:
            a, b = random.sample(range(len(population)-1), 2)
            offspring = crossover_operator.execute([population[a], population[b]])
            population.append(offspring[0])

        logger.info('Population incremented to: {0}'.format(len(population)))

        return population

    def create_solution(self) -> None:
        raise NotImplemented

    def get_name(self) -> str:
        return "Multiple Sequence Alignment (MSA) BaliBASE problem"
