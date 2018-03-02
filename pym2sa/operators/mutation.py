import random

from jmetal.core.operator import Mutation

from pym2sa.core.solution import MSASolution


class OneRandomGapInsertion(Mutation[MSASolution]):

    def __init__(self, probability: float, remove_gap_columns: bool = False) -> None:
        super(OneRandomGapInsertion, self).__init__(probability=probability)
        self.remove_gap_columns = remove_gap_columns

    def execute(self, solution: MSASolution) -> MSASolution:
        if solution is None:
            raise Exception("Solution is none")

        return self.do_mutation(solution)

    def do_mutation(self, solution: MSASolution) -> MSASolution:
        if random.random() <= self.probability:
            for seq_index in range(0, solution.number_of_variables):
                point = random.randint(0, solution.get_length_of_alignment() - 1)
                solution.add_gap_to_sequence(seq_index, point)

        return solution
