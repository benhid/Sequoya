import random

from jmetal.core.operator import Mutation

from pym2sa.core.solution import MSASolution


class RandomGapInsertion(Mutation[MSASolution]):
    """ Inserts a gap in a random position for each sequence. """

    __GAP_SYMBOL = '-'

    def __init__(self, probability: float) -> None:
        super(RandomGapInsertion, self).__init__(probability=probability)

    def execute(self, solution: MSASolution) -> MSASolution:
        if random.random() <= self.probability:
            for i in range(solution.number_of_variables):
                solution_as_list = list(solution.variables[i])
                solution_as_list.insert(random.randint(0, len(solution.variables[i])-1), self.__GAP_SYMBOL)
                solution.variables[i] = ''.join(solution_as_list)

        return solution
