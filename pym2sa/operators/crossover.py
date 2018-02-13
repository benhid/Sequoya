from typing import List
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
            raise Exception("The number of parents is not two (2) but " + str(len(parents)))

        offspring = [copy.deepcopy(parents[0]), copy.deepcopy(parents[1])]

        print("bef", offspring[0].decode_alignment_as_list_of_pairs())
        print("bef", offspring[1].decode_alignment_as_list_of_pairs())

        if random.random() <= self.probability:
            cx_point = random.randint(0, parents[0].get_length_of_alignment() - 1)
            self.split(cx_point, offspring[0])

            print("after", cx_point, offspring[0].decode_alignment_as_list_of_pairs())
            print("after", cx_point, offspring[1].decode_alignment_as_list_of_pairs())

        return offspring

    def split(self, cx_point: int, parent: MSASolution):
        new_parent = parent

        print('1. new parent start', new_parent.gaps_groups)

        new_parent.split_gap_column(cx_point)

        print('2. new parent after split', new_parent.gaps_groups)

        parent_gaps_groups = new_parent.gaps_groups
        parent_encoded_sequences = new_parent.variables

        for i in range(len(parent_gaps_groups)):
            gaps_group = parent_gaps_groups[i]

            print('3. gaps group', gaps_group)

            if cx_point in gaps_group:
                idx_occurrence = gaps_group.index(cx_point)  # find index of first occurrence
                print('3.1.1. cx_point (', cx_point, ') in gaps group, ocurrence', idx_occurrence, idx_occurrence % 2)
                if idx_occurrence % 2:
                    new_parent.gaps_groups[i] = gaps_group[:idx_occurrence + 1]
                else:
                    new_parent.gaps_groups[i] = gaps_group[:idx_occurrence]
                print('3.1.2. new gaps group', new_parent.gaps_groups[i])
            else:
                # in this case, the cx_point is not in any gap group, so we'll have to delete every gaps group
                # after the cx_point; else, we will do nothing
                # e.g. if the gaps group is [3,7] and cx_point is 2, our new gaps group will be empty ( [] )
                #      else, if the gaps group is [1, 1, 3, 6] and cx_point is 8, our new gaps group will be the same
                print('3.2.1. cx_point NOT in gaps group')
                inf_gaps_group = []

                for j in range(0, len(gaps_group) - 1, 2):
                    if gaps_group[j + 1] <= cx_point:
                        print(gaps_group[j + 1])
                        inf_gaps_group.append(gaps_group[j])
                        inf_gaps_group.append(gaps_group[j + 1])
                    else:
                        break

                new_parent.gaps_groups[i] = inf_gaps_group
                print('3.2.2. new gaps group', new_parent.gaps_groups[i])

            sequence = parent_encoded_sequences[i]
            number_of_gaps = [len(gaps_group) for gaps_group in new_parent.gaps_groups]

            print('4. new parent sequence before', new_parent.variables[i])
            seq_cx_point = len(sequence) - sum(number_of_gaps) + 1
            new_parent.variables[i] = sequence[:seq_cx_point] if seq_cx_point >= 0 else ''
            print('5. new parent sequence after', new_parent.variables[i])

        return new_parent

    def merge(self, parent_a: MSASolution, parent_b: MSASolution):
        pass

    def get_number_of_parents(self) -> int:
        return 2
