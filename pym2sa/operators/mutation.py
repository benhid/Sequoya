import random

from jmetal.core.operator import Mutation

from pym2sa.core.solution import MSASolution


class ShiftClosedGaps(Mutation[MSASolution]):
    """ Selects a random group and merges it with the adjacent gaps group. """

    def __init__(self, probability: float) -> None:
        super(ShiftClosedGaps, self).__init__(probability=probability)

    def execute(self, solution: MSASolution) -> MSASolution:
        if solution is None:
            raise Exception("Solution is none")

        return self.do_mutation(solution)

    def do_mutation(self, solution: MSASolution) -> MSASolution:
        if random.random() <= self.probability:
            # Select one random sequence from all
            """
            if solution.number_of_variables > 1:
                seq = random.randint(0, solution.number_of_variables - 1)
            else:
                seq = 0
            """
            for seq in range(solution.number_of_variables):

                gaps_group = solution.gaps_groups[seq]

                if len(gaps_group) >= 4:
                    random_gaps_group = random.randrange(0, len(gaps_group) - 2, 2)

                    to_add = gaps_group[random_gaps_group + 3] - gaps_group[random_gaps_group + 2] + 1
                    gaps_group[random_gaps_group + 1] += to_add

                    del gaps_group[random_gaps_group + 3]
                    del gaps_group[random_gaps_group + 2]

                # Sanity check: alignment is valid (same length for all sequences)
                if not solution.is_valid():
                    raise Exception("Mutated solution is not valid! {0}".format(solution.decode_alignment_as_list_of_pairs()))

        return solution


class TwoRandomAdjacentGapGroup(Mutation[MSASolution]):
    """ Selects a random group and merges it with the adjacent gaps group. """

    def __init__(self, probability: float, remove_gap_columns: bool = True) -> None:
        super(TwoRandomAdjacentGapGroup, self).__init__(probability=probability)
        self.remove_full_of_gap_columns = remove_gap_columns

    def execute(self, solution: MSASolution) -> MSASolution:
        if solution is None:
            raise Exception("Solution is none")

        return self.do_mutation(solution)

    def do_mutation(self, solution: MSASolution) -> MSASolution:
        if random.random() <= self.probability:
            # Select one random sequence from all
            if solution.number_of_variables > 1:
                seq = random.randint(0, solution.number_of_variables - 1)
            else:
                seq = 0

            gaps_group = solution.gaps_groups[seq]

            if len(gaps_group) >= 4:
                random_gaps_group = random.randrange(0, len(gaps_group) - 2, 2)

                right_is_closest = False

                if right_is_closest:
                    to_add = gaps_group[random_gaps_group + 3] - gaps_group[random_gaps_group + 2] + 1
                    gaps_group[random_gaps_group + 1] += to_add

                    del gaps_group[random_gaps_group + 3]
                    del gaps_group[random_gaps_group + 2]

            if self.remove_full_of_gap_columns:
                solution.remove_full_of_gaps_columns()

            # Sanity check: alignment is valid (same length for all sequences)
            if not solution.is_valid():
                raise Exception("Mutated solution is not valid! {0}".format(solution.decode_alignment_as_list_of_pairs()))

        return solution


class OneRandomGapInsertion(Mutation[MSASolution]):

    def __init__(self, probability: float, remove_gap_columns: bool = False) -> None:
        super(OneRandomGapInsertion, self).__init__(probability=probability)
        self.remove_full_of_gap_columns = remove_gap_columns

    def execute(self, solution: MSASolution) -> MSASolution:
        if solution is None:
            raise Exception("Solution is none")

        return self.do_mutation(solution)

    def do_mutation(self, solution: MSASolution) -> MSASolution:
        if random.random() <= self.probability:
            length_of_alignment = solution.get_length_of_alignment()
            for seq_index in range(solution.number_of_variables):
                point = random.randint(0, length_of_alignment - 1)
                solution.add_gap_to_sequence(seq_index, point)

            if self.remove_full_of_gap_columns:
                solution.remove_full_of_gaps_columns()

            # Sanity check: alignment is valid (same length for all sequences)
            if not solution.is_valid():
                raise Exception("Mutated solution is not valid! {0}".format(solution.decode_alignment_as_list_of_pairs()))

        return solution
