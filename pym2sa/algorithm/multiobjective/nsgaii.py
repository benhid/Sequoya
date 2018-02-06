import logging
import os
from typing import List, TypeVar

from jmetal.algorithm.multiobjective.nsgaii import NSGAII
from jmetal.core.operator import Mutation, Selection, Crossover
from jmetal.core.problem import Problem

from pym2sa.core.solution import MSASolution
from pymsa.util.fasta import read_fasta_file_as_list_of_pairs

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

S = TypeVar('S')
R = TypeVar(List[S])


class NSGA2MSA(NSGAII[S, R]):
    def __init__(self, problem: Problem[S], population_size: int, initial_population_path: str, max_evaluations: int, mutation: Mutation[S],
                 crossover: Crossover[S, S], selection: Selection[List[S], S]):
        super().__init__(problem, population_size, max_evaluations, mutation, crossover, selection)
        self.initial_population_path = initial_population_path

    def create_initial_population(self) -> List[MSASolution]:
        """ In this case, the initial population is not randomly generated but provided
        (from precomputed alignments). """

        population = []

        try:
            for file in os.listdir(self.initial_population_path):
                logger.info('Reading file ' + self.initial_population_path + file)

                fasta_file = read_fasta_file_as_list_of_pairs(file, self.initial_population_path)

                msa = MSASolution(number_of_variables=len(fasta_file), number_of_objectives=2, number_of_constraints=0)
                msa.header = list(pair[0] for pair in fasta_file)
                msa.variables = list(pair[1] for pair in fasta_file)

                population.append(msa)
        except FileNotFoundError:
            raise Exception("Invalid path provided: {0}".format(self.initial_population_path))

        if len(population) < 2:
            raise Exception("More than one precomputed alignment is needed!")

        return population
