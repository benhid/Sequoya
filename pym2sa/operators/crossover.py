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
            cx_point = random.randint(0, parents[0].get_length_of_original_alignment() - 1)

            # split parent into two groups (first half, second half)
            parent_a_first_half, parent_a_second_half = self.split(cx_point, offspring[0])
            parent_b_first_half, parent_b_second_half = self.split(cx_point, offspring[1])

            # merge both parents
            offspring = self.merge([parent_a_first_half, parent_b_second_half])

            for i in range(len(offspring)):
                if not offspring[i].is_valid():
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

    def merge(self, parents: List[MSASolution]) -> List[MSASolution]:
        offspring = [copy.deepcopy(parents[0]), copy.deepcopy(parents[1])]

        for i in range(parents[0].number_of_variables):
            length_of_gaps_parent_zero = parents[0].get_length_of_gaps(i)
            offset_parent_zero = len(offspring[0].variables[i]) + length_of_gaps_parent_zero

            # merge encoded sequences
            offspring[0].variables[i] = parents[0].variables[i] + parents[1].variables[i]

            # merge gaps groups
            offspring[0].gaps_groups[i] = parents[0].gaps_groups[i] + [gap + offset_parent_zero for gap in parents[1].gaps_groups[i]]

        for i in range(parents[1].number_of_variables):
            length_of_gaps_parent_one = parents[1].get_length_of_gaps(i)
            offset_parent_one = len(offspring[1].variables[i]) + length_of_gaps_parent_one

            # merge encoded sequences
            offspring[1].variables[i] = parents[1].variables[i] + parents[0].variables[i]

            # merge gaps groups
            offspring[1].gaps_groups[i] = parents[1].gaps_groups[i] + [gap + offset_parent_one for gap in
                                                                       parents[0].gaps_groups[i]]

        return offspring

    def get_number_of_parents(self) -> int:
        return 2


class SinglePointXMSA(Crossover[MSASolution, MSASolution]):
    """ Implements a single point crossover for MSA representation. """

    def __init__(self, probability: float) -> None:
        if not 0 <= probability <= 1:
            raise Exception("Crossover probability value invalid: " + str(probability))
        super(SinglePointXMSA, self).__init__(probability=probability)

    def execute(self, parents: List[MSASolution]) -> List[MSASolution]:
        if len(parents) != 2:
            raise Exception("The number of parents is not two (2) but " + str(len(parents)))

        offspring = [copy.deepcopy(parents[0]), copy.deepcopy(parents[1])]

        if random.random() <= self.probability:
            # get random cut point
            cx_point = random.randint(0, parents[0].get_length_of_original_alignment() - 1)

            # split parents
            a1, a2, b1, b2 = self.split_parents(cx_point, offspring[0], offspring[1])

            # merge to children
            offspring[0] = self.merge_parents(a1, b2)
            offspring[1] = self.merge_parents(b1, a2)

            # check if the alignment of the children are valid
            #for i in range(len(offspring)):
            #    if not offspring[i].is_valid():
            #       raise Exception("Fatal error in crossover operator: children don't have proper alignments")

        return offspring

    def split_parents(self, cx_point: int, parent_a: MSASolution, parent_b: MSASolution):
        # Copy both parents
        new_parent_a = copy.deepcopy(parent_a)
        new_parent_a.split_gap_column(cx_point)

        new_parent_b = copy.deepcopy(parent_b)
        new_parent_b.split_gap_column(cx_point)

        # the number of gaps groups in both parents should be the same at this point
        number_of_gaps_groups = len(new_parent_a.gaps_groups)

        a1, a2 = copy.deepcopy(new_parent_a), copy.deepcopy(new_parent_a)
        b1, b2 = copy.deepcopy(new_parent_b), copy.deepcopy(new_parent_b)

        # split gaps groups of parent A
        for i in range(number_of_gaps_groups):
            gaps_group = a1.gaps_groups[i]

            if cx_point in gaps_group:
                # find index of first occurrence
                idx_occurrence = gaps_group.index(cx_point)

                if idx_occurrence % 2:
                    a1.gaps_groups[i] = gaps_group[:idx_occurrence + 1]
                    a2.gaps_groups[i] = gaps_group[idx_occurrence + 1:]
                else:
                    if gaps_group[idx_occurrence] == gaps_group[idx_occurrence + 1] == cx_point:
                        a1.gaps_groups[i] = gaps_group[:idx_occurrence + 2]
                        a2.gaps_groups[i] = gaps_group[idx_occurrence + 2:]
                    else:
                        a1.gaps_groups[i] = gaps_group[:idx_occurrence]
                        a2.gaps_groups[i] = gaps_group[idx_occurrence:]
            else:
                inf_gaps_group, sup_gaps_group = [], []

                for j in range(0, len(gaps_group) - 1, 2):
                    if gaps_group[j + 1] <= cx_point:
                        inf_gaps_group.append(gaps_group[j])
                        inf_gaps_group.append(gaps_group[j + 1])
                    else:
                        sup_gaps_group.append(gaps_group[j])
                        sup_gaps_group.append(gaps_group[j + 1])

                    a1.gaps_groups[i] = inf_gaps_group
                    a2.gaps_groups[i] = sup_gaps_group

            a2.gaps_groups[i] = [element + 1 - a1.get_length_of_sequence(i) for element in a2.gaps_groups[i]]

        # split sequences of parent A
        for i in range(number_of_gaps_groups):
            sequence_cx_point = self.get_sequence_cx_point(a1.gaps_groups[i], cx_point)

            a1.variables[i] = new_parent_a.variables[i][:sequence_cx_point]
            a2.variables[i] = new_parent_a.variables[i][sequence_cx_point:]

        print('original a1', parent_a.variables, parent_a.gaps_groups)
        print('a1', a1.variables, a1.gaps_groups)
        print('a2', a2.variables, a2.gaps_groups)

        # split gaps groups of parent B based on parent A
        for i in range(number_of_gaps_groups):
            gaps_group = b1.gaps_groups[i]

            if cx_point in gaps_group:
                # find index of first occurrence
                idx_occurrence = gaps_group.index(cx_point)

                if idx_occurrence % 2:
                    b1.gaps_groups[i] = gaps_group[:idx_occurrence + 1]
                    b2.gaps_groups[i] = gaps_group[idx_occurrence + 1:]
                else:
                    if gaps_group[idx_occurrence] == gaps_group[idx_occurrence + 1] == cx_point:
                        b1.gaps_groups[i] = gaps_group[:idx_occurrence + 2]
                        b2.gaps_groups[i] = gaps_group[idx_occurrence + 2:]
                    else:
                        b1.gaps_groups[i] = gaps_group[:idx_occurrence]
                        b2.gaps_groups[i] = gaps_group[idx_occurrence:]
            else:
                inf_gaps_group, sup_gaps_group = [], []

                for j in range(0, len(gaps_group) - 1, 2):
                    if gaps_group[j + 1] <= cx_point:
                        inf_gaps_group.append(gaps_group[j])
                        inf_gaps_group.append(gaps_group[j + 1])
                    else:
                        sup_gaps_group.append(gaps_group[j])
                        sup_gaps_group.append(gaps_group[j + 1])

                    b1.gaps_groups[i] = inf_gaps_group
                    b2.gaps_groups[i] = sup_gaps_group

            b2.gaps_groups[i] = [element + 1 - b1.get_length_of_sequence(i) for element in b2.gaps_groups[i]]

        # split sequences of parent B
        for i in range(number_of_gaps_groups):
            # cut according to the same position of the parent A
            sequence_from_parent_a = a1.variables[i]
            number_of_elements_to_keep_from_first_parent = len(sequence_from_parent_a)

            b1.variables[i] = new_parent_b.variables[i][:number_of_elements_to_keep_from_first_parent]
            b2.variables[i] = new_parent_b.variables[i][number_of_elements_to_keep_from_first_parent:]

        print('original b1', parent_b.variables, parent_b.gaps_groups)
        print('b1', b1.variables, b1.gaps_groups)
        print('b2', b2.variables, b2.gaps_groups)

        return [a1, a2, b1, b2]

    def merge_parents(self, parent_a: MSASolution, parent_b: MSASolution) -> MSASolution:
        individual = copy.deepcopy(parent_a)

        for i in range(parent_a.number_of_variables):
            offset_parent_zero = parent_a.get_length_of_sequence(i)

            # merge encoded sequences
            individual.variables[i] = parent_a.variables[i] + parent_b.variables[i]

            # merge gaps groups
            individual.gaps_groups[i] = parent_a.gaps_groups[i] + \
                                        [gap + offset_parent_zero for gap in parent_b.gaps_groups[i]]

        return individual

    def get_sequence_cx_point(self, gaps_group, cx_point) -> int:
        position_to_discount = 0

        # equivalent to:
        #   sum([gaps_group[j+1]-gaps_group[j] for j in range(0, len(gaps_group) - 1, 2)])
        for j in range(0, len(gaps_group) - 1, 2):
            position_to_discount += gaps_group[j + 1] - gaps_group[j] + 1

        pos = cx_point - position_to_discount + 1
        sequence_cx_point = pos if pos > 0 else 0

        return sequence_cx_point

    def get_number_of_parents(self) -> int:
        return 2
