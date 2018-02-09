from typing import List, Tuple, Any
import copy
import random

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

        print("bef", offspring[0].decode_alignment_as_list_of_pairs())
        print("bef", offspring[1].decode_alignment_as_list_of_pairs())

        if random.random() <= self.probability:
            cx_point = random.randint(0, parents[0].get_length_of_alignment() - 1)
            self.split(cx_point, offspring[0])
            self.split(cx_point, offspring[1])

            print("after", cx_point, offspring[0].decode_alignment_as_list_of_pairs())
            print("after", cx_point, offspring[1].decode_alignment_as_list_of_pairs())

        return offspring

    def split(self, cx_point: int, parent: MSASolution):
        parent.split_gap_column(cx_point)

        parent_gaps_groups = parent.get_gaps_groups()
        parent_encoded_sequences = parent.variables

        for i in range(len(parent_gaps_groups)):
            gaps_group = parent_gaps_groups[i]

            if cx_point in gaps_group:
                idx_occurrence = gaps_group.index(cx_point)  # find index of first occurrence

                if idx_occurrence % 2 == 0:  # e.g. [ [0, 1, 2, 5] ], cx_point = 2
                    parent.gaps_groups[i] = gaps_group[:idx_occurrence]
                else:  # e.g. [ [1, 2, 2, 5] ], cx_point = 2
                    parent.gaps_groups[i] = gaps_group[:idx_occurrence + 1]
            else:
                # in this case, the cx_point is not in any gap group, so we'll have to delete every gaps group
                # after the cx_point; else, we will do nothing
                # e.g. if the gaps group is [3,7] and cx_point is 2, our new gaps group will be empty ( [] )
                #      else, if the gaps group is [1, 1, 3, 6] and cx_point is 8, our new gaps group will be the same
                parent.gaps_groups[i] = [] if gaps_group[1] > cx_point else gaps_group

            sequence = parent_encoded_sequences[i]
            number_of_gaps = parent.get_number_of_gaps_of_sequence(i)

            parent.variables[i] = sequence[:cx_point - number_of_gaps + 1]

    def get_number_of_parents(self) -> int:
        return 2
