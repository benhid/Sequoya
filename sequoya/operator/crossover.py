import copy
import random
from typing import List

from jmetal.core.operator import Crossover

from sequoya.core.solution import MSASolution


class SPXMSA(Crossover[List[MSASolution], MSASolution]):
    """ Implements a single point crossover for MSA representation. """

    def __init__(self, probability: float, remove_gap_columns: bool = True) -> None:
        super(SPXMSA, self).__init__(probability=probability)

        self.remove_full_of_gap_columns = remove_gap_columns
        self.has_solution_been_crossed = None

    def execute(self, parents: List[MSASolution]) -> List[MSASolution]:
        if len(parents) != 2:
            raise Exception('The number of parents is not two (2) but {}'.format(len(parents)))

        return self.do_crossover(parents)

    def do_crossover(self, parents: List[MSASolution]) -> List[MSASolution]:
        if random.random() <= self.probability:
            cx_point = random.randint(1, parents[0].get_length_of_alignment() - 1)

            column_positions_in_first_parent = \
                self.find_original_positions_in_original_sequences(parents[0], cx_point)
            column_positions_in_second_parent = \
                self.__find_original_positions_in_aligned_sequences(parents[1], column_positions_in_first_parent)
            cutting_points_in_first_parent = \
                self.find_cutting_points_in_first_parent(parents[0], cx_point)

            offspring = self.cross_parents(
                cx_point, parents, cutting_points_in_first_parent, column_positions_in_second_parent)

            self.has_solution_been_crossed = True

            offspring[0].merge_gaps_groups()
            offspring[1].merge_gaps_groups()

            if self.remove_full_of_gap_columns:
                offspring[0].remove_full_of_gaps_columns()
                offspring[1].remove_full_of_gaps_columns()
        else:
            offspring = copy.deepcopy(parents)
            self.has_solution_been_crossed = False

        return offspring

    def cross_parents(self, cx_point: int, parents: List[MSASolution], cutting_points_in_first_parent: list,
                      column_positions_in_second_parent: list) -> List[MSASolution]:
        offspring_1 = copy.deepcopy(parents[0])
        offspring_2 = copy.deepcopy(parents[1])

        for i in range(offspring_1.number_of_variables):
            new_gap_group_list = []

            if cutting_points_in_first_parent[i] != -1:
                size_adjustment = cutting_points_in_first_parent[i] - column_positions_in_second_parent[i]
                cutting_point_passed = False

                # offspring 1: Add the gap groups of the first parent before the cutting point
                while not cutting_point_passed:
                    gaps_group = parents[0].gaps_groups[i]

                    for j in range(0, len(gaps_group) - 1, 2):
                        if gaps_group[j + 1] < cutting_points_in_first_parent[i]:
                            new_gap_group_list.append(gaps_group[j])
                            new_gap_group_list.append(gaps_group[j + 1])
                        else:
                            cutting_point_passed = True
                    cutting_point_passed = True

                # offspring 1: Add the gap groups on the first parent after the cutting point
                gaps_group = parents[1].gaps_groups[i]
                for j in range(0, len(gaps_group) - 1, 2):
                    if gaps_group[j + 1] < column_positions_in_second_parent[i]:
                        pass
                    else:
                        new_gap_group_list.append(gaps_group[j] + size_adjustment)
                        new_gap_group_list.append(gaps_group[j + 1] + size_adjustment)

                offspring_1.gaps_groups[i] = new_gap_group_list
                cutting_point_passed = False
                new_gap_group_list = []

                # offspring 2: Add the gap groups of the second parent before the cutting point
                while not cutting_point_passed:
                    gaps_group = parents[1].gaps_groups[i]

                    for j in range(0, len(gaps_group) - 1, 2):
                        if gaps_group[j + 1] < column_positions_in_second_parent[i]:
                            new_gap_group_list.append(gaps_group[j])
                            new_gap_group_list.append(gaps_group[j + 1])
                        else:
                            cutting_point_passed = True
                    cutting_point_passed = True

                # offspring 2: Add the gap groups on the second parent after the cutting point
                gaps_group = parents[0].gaps_groups[i]
                for j in range(0, len(gaps_group) - 1, 2):
                    if gaps_group[j + 1] < cutting_points_in_first_parent[i]:
                        pass
                    else:
                        new_gap_group_list.append(gaps_group[j] - size_adjustment)
                        new_gap_group_list.append(gaps_group[j + 1] - size_adjustment)

                offspring_2.gaps_groups[i] = new_gap_group_list

        max_sequence_length = self.find_length_of_the_largest_sequence(offspring_1)
        self.fill_sequences_with_gaps_to_reach_the_max_sequence_length(
            offspring_1, max_sequence_length, cutting_points_in_first_parent)

        max_sequence_length = self.find_length_of_the_largest_sequence(offspring_2)
        self.fill_sequences_with_gaps_to_reach_the_max_sequence_length(
            offspring_2, max_sequence_length, column_positions_in_second_parent)

        # Sanity check: alignment is valid (same length for all sequences)
        if not offspring_1.is_valid_msa():
            raise Exception('Offspring 1 is not valid! \n {0} \n {1} \n {2} \n {3} \n {4} \n {5} \n {6}'
                            .format(parents[0].decode_alignment_as_list_of_sequences(),
                                    parents[1].decode_alignment_as_list_of_sequences(), cx_point,
                                    offspring_1.decode_alignment_as_list_of_sequences(),
                                    offspring_2.decode_alignment_as_list_of_sequences(),
                                    cutting_points_in_first_parent, column_positions_in_second_parent))
        if not offspring_2.is_valid_msa():
            raise Exception('Offspring 2 is not valid! \n {0} \n {1} \n {2} \n {3} \n {4} \n {5} \n {6}'
                            .format(parents[0].decode_alignment_as_list_of_sequences(),
                                    parents[1].decode_alignment_as_list_of_sequences(), cx_point,
                                    offspring_1.decode_alignment_as_list_of_sequences(),
                                    offspring_2.decode_alignment_as_list_of_sequences(),
                                    cutting_points_in_first_parent, column_positions_in_second_parent))

        return [offspring_1, offspring_2]

    def find_original_positions_in_original_sequences(self, solution: MSASolution, column: int) -> list:
        """ Given a solution, find for each sequence the original positions of the symbol in the column
        in the original unaligned sequences """
        positions = [-1 for _ in range(solution.number_of_variables)]

        for i in range(solution.number_of_variables):
            positions[i] = self.__find_symbol_position_in_original_sequence(solution, i, column)

        return positions

    def __find_symbol_position_in_original_sequence(self, solution: MSASolution, seq_index: int, position: int):
        """ Given a symbol position, finds the corresponding position of the symbol in the original
        sequence if gaps are not taken into account. If the symbol is a gap the returned value is -1 """
        if position > solution.get_length_of_alignment():
            raise Exception('Position {0} is larger than the sequence size {1}'.format(position,
                                                                                       solution.get_length_of_alignment()))

        if not solution.is_gap_char_at_sequence(seq_index, position):
            symbol_position = solution.get_char_position_in_original_sequence(seq_index, position)
        else:
            position = solution.get_next_char_position_after_gap(seq_index, position)
            if position < 0:
                symbol_position = -1
            else:
                symbol_position = solution.get_char_position_in_original_sequence(seq_index, position)

        return symbol_position

    def __find_original_positions_in_aligned_sequences(self, solution: MSASolution,
                                                       column_positions_in_first_parent: list):
        positions = [-1 for _ in range(solution.number_of_variables)]

        for i in range(solution.number_of_variables):
            pos = column_positions_in_first_parent[i]
            positions[i] = solution.get_original_char_position_in_aligned_sequence(i, pos)

        return positions

    def find_cutting_points_in_first_parent(self, solution: MSASolution, position: int) -> list:
        """ Find the real cutting points in a solution. If the column is a gap then the next non-gap
        symbol must be found """
        positions = [-1 for _ in range(solution.number_of_variables)]

        for i in range(solution.number_of_variables):
            if solution.is_gap_char_at_sequence(i, position):
                positions[i] = solution.get_next_char_position_after_gap(i, position)
            else:
                positions[i] = position

        return positions

    def find_length_of_the_largest_sequence(self, solution: MSASolution):
        max_length = solution.get_length_of_sequence(0)

        for i in range(1, solution.number_of_variables):
            length_of_sequence_i = solution.get_length_of_sequence(i)
            if max_length < length_of_sequence_i:
                max_length = length_of_sequence_i

        return max_length

    def fill_sequences_with_gaps_to_reach_the_max_sequence_length(self, solution: MSASolution, max_length: int,
                                                                  cutting_points: list):
        for i in range(solution.number_of_variables):
            sequence_length = solution.get_length_of_sequence(i)
            if sequence_length != max_length:
                for j in range(sequence_length, max_length):
                    if cutting_points[i] == -1:
                        solution.add_gap_to_sequence_at_index(seq_index=i, gap_position=sequence_length - 1)
                    else:
                        solution.add_gap_to_sequence_at_index(seq_index=i, gap_position=cutting_points[i] + 1)

    def get_number_of_parents(self) -> int:
        return 2

    def get_number_of_children(self) -> int:
        return 2

    def get_name(self) -> str:
        return 'Simplex Crossover for MSA'
