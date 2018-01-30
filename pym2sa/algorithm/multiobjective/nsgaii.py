import logging
import os
from typing import List, TypeVar

from jmetal.algorithm.multiobjective.nsgaii import NSGAII
from pym2sa.core.solution import MSASolution
from pymsa.util import read_fasta_file_as_list_of_pairs

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

S = TypeVar('S')
R = TypeVar(List[S])


class NSGA2MSA(NSGAII[S, R]):

    def create_initial_population(self) -> List[MSASolution]:
        """ In this case, the initial population is not randomly generated but provided (precomputed alignments). """

        population = []
        directory = os.path.abspath(
            os.path.join(os.path.dirname(__file__), '..', '..', '..', 'pym2sa/runner/files')) + '/'

        try:
            for file in os.listdir(directory):
                logger.info('Reading file ' + directory + file)

                fasta_file = read_fasta_file_as_list_of_pairs(file, directory)

                msa = MSASolution(number_of_variables=len(fasta_file), number_of_objectives=2, number_of_constraints=0)
                msa.header = list(pair[0] for pair in fasta_file)
                msa.variables = list(pair[1] for pair in fasta_file)

                population.append(msa)
        except FileNotFoundError:
            raise Exception("Invalid path provided:", directory)

        if len(population) < 2:
            raise Exception("More than one precomputed alignment is needed!")

        return population
