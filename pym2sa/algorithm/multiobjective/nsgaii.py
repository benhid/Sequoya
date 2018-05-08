import logging
from typing import List, TypeVar

from jmetal.algorithm.multiobjective.nsgaii import NSGAII
from jmetal.core.operator import Mutation, Selection, Crossover
from jmetal.core.problem import Problem

from pym2sa.core.solution import MSASolution

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

S = TypeVar('S')
R = TypeVar(List[S])


class NSGA2MSA(NSGAII[S, R]):

    def __init__(self, problem: Problem[S], population_size: int, max_evaluations: int,
                 mutation: Mutation[S], crossover: Crossover[S, S], selection: Selection[List[S], S]):
        super().__init__(problem, population_size, max_evaluations, mutation, crossover, selection)

    def create_initial_population(self) -> List[MSASolution]:
        population = self.problem.create_solutions(self.population_size)
        return population
