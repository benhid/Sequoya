from typing import List, TypeVar

from jmetal.algorithm.multiobjective.nsgaii import NSGAII
from jmetal.component.evaluator import SequentialEvaluator, Evaluator
from jmetal.core.operator import Mutation, Selection, Crossover
from jmetal.core.problem import Problem

from pym2sa.core.solution import MSASolution

S = TypeVar('S')
R = TypeVar(List[S])


class NSGA2MSA(NSGAII[S, R]):

    def __init__(self,
                 problem: Problem[S],
                 population_size: int,
                 max_evaluations: int,
                 mutation: Mutation[S],
                 crossover: Crossover[S, S],
                 selection: Selection[List[S], S],
                 evaluator: Evaluator = SequentialEvaluator()):
        super().__init__(problem, population_size, max_evaluations, mutation, crossover, selection, evaluator)

    def create_initial_population(self) -> List[MSASolution]:
        return self.problem.import_instance(self.population_size)
