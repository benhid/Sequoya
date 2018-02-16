from typing import List
import copy
import random

from jmetal.core.operator import Crossover

from pym2sa.core.solution import MSASolution


class SinglePointMSA(Crossover[MSASolution, MSASolution]):
    """ Implements a single point crossover for MSA representation. """

    def __init__(self, probability: float) -> None:
        if not 0 < probability < 1:
            raise Exception("Crossover probability value invalid: " + str(probability))
        super(SinglePointMSA, self).__init__(probability=probability)

    def execute(self, parents: List[MSASolution]) -> List[MSASolution]:
        if len(parents) != 2:
            raise Exception("The number of parents is not two (2) but " + str(len(parents)))

        offspring = [copy.deepcopy(parents[0]), copy.deepcopy(parents[1])]

        if random.random() <= self.probability:
            # get random cut point
            cx_point = random.randint(0, parents[0].get_length_of_alignment() - 1)

            # split parent into two groups (first half, second half)
            parent_a_first_half, parent_a_second_half = self.split(cx_point, offspring[0])
            parent_b_first_half, parent_b_second_half = self.split(cx_point, offspring[1])

            # merge both parents
            offspring[0] = self.merge(parent_a_first_half, parent_b_second_half)
            offspring[1] = self.merge(parent_b_first_half, parent_a_second_half)

            if not offspring[0].is_valid() or not offspring[1].is_valid():
                raise Exception("Fatal error in crossover operator")

        return offspring

    def split(self, cx_point: int, parent: MSASolution):
        new_parent = copy.deepcopy(parent)
        new_parent.split_gap_column(cx_point)

        new_parent_first_half, new_parent_second_half = copy.deepcopy(new_parent), copy.deepcopy(new_parent)

        original_parent_gaps_groups = new_parent.gaps_groups
        original_parent_encoded_sequences = new_parent.variables

        # split gaps groups
        for i in range(len(original_parent_gaps_groups)):
            gaps_group = original_parent_gaps_groups[i]

            if cx_point in gaps_group:
                # find index of first occurrence
                idx_occurrence = gaps_group.index(cx_point)

                if idx_occurrence % 2:
                    new_parent_first_half.gaps_groups[i] = gaps_group[:idx_occurrence + 1]
                    new_parent_second_half.gaps_groups[i] = gaps_group[idx_occurrence + 1:]
                else:
                    if gaps_group[idx_occurrence] == gaps_group[idx_occurrence + 1] == cx_point:
                        new_parent_first_half.gaps_groups[i] = gaps_group[:idx_occurrence + 2]
                        new_parent_second_half.gaps_groups[i] = gaps_group[idx_occurrence + 2:]
                    else:
                        new_parent_first_half.gaps_groups[i] = gaps_group[:idx_occurrence]
                        new_parent_second_half.gaps_groups[i] = gaps_group[idx_occurrence:]
            else:
                inf_gaps_group, sup_gaps_group = [], []

                for j in range(0, len(gaps_group) - 1, 2):
                    if gaps_group[j + 1] <= cx_point:
                        inf_gaps_group.append(gaps_group[j])
                        inf_gaps_group.append(gaps_group[j + 1])
                    else:
                        sup_gaps_group.append(gaps_group[j])
                        sup_gaps_group.append(gaps_group[j + 1])

                    new_parent_first_half.gaps_groups[i] = inf_gaps_group
                    new_parent_second_half.gaps_groups[i] = sup_gaps_group

        # split sequences
        for i in range(len(new_parent_first_half.gaps_groups)):
            gaps_group = new_parent_first_half.gaps_groups[i]
            position_to_discount = 0

            # equivalent: sum([gaps_group[j+1]-gaps_group[j] for j in range(0, len(gaps_group) - 1, 2)])
            for j in range(0, len(gaps_group) - 1, 2):
                position_to_discount += gaps_group[j + 1] - gaps_group[j] + 1

            pos = cx_point - position_to_discount + 1
            sequence_cx_point = pos if pos > 0 else 0

            new_parent_first_half.variables[i] = original_parent_encoded_sequences[i][:sequence_cx_point]
            new_parent_second_half.variables[i] = original_parent_encoded_sequences[i][sequence_cx_point:]

        return [new_parent_first_half, new_parent_second_half]

    def merge(self, parent_a: MSASolution, parent_b: MSASolution) -> MSASolution:
        pass

    def get_number_of_parents(self) -> int:
        return 2
