import copy
import random
from typing import List

from jmetal.core.operator import Crossover

from pym2sa.core.solution import MSASolution


class SinglePointMSA(Crossover[MSASolution, MSASolution]):
    """ Implements a single point crossover for MSA representation. """

    def __init__(self, probability: float) -> None:
        if probability < 0 or probability > 1:
            raise Exception("Crossover probability value invalid: " + str(probability))

        super(SinglePointMSA, self).__init__(probability=probability)

    def execute(self, parents: List[MSASolution]) -> List[MSASolution]:
        if len(parents) != 2:
            raise Exception("The number of parents is not two: " + str(len(parents)))

        offspring = [copy.deepcopy(parents[0]), copy.deepcopy(parents[1])]

        if random.random() <= self.probability:
            cx_point = random.randint(1, offspring[0].get_lengh_of_alignment() - 1)

            left, right = [], []
            for solution in offspring:
                for i in range(0, solution.number_of_variables):
                    left.append(solution.variables[i][:cx_point])
                    right.append(solution.variables[i][cx_point:])

            for i in range(0, parents[0].number_of_variables):
                offspring[0].variables[i] = left[i] + right[i+parents[0].number_of_variables]
                offspring[1].variables[i] = left[i+parents[0].number_of_variables] + right[i]

        return offspring

    def get_number_of_parents(self) -> int:
        return 2