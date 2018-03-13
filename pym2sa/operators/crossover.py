from typing import List
import copy
import random

from jmetal.core.operator import Crossover

from pym2sa.core.solution import MSASolution


class GapSequenceSolutionSinglePoint(Crossover[MSASolution, MSASolution]):
    """ Implements a single point crossover for MSA representation. """

    def __init__(self, probability: float, remove_gap_columns: bool = False) -> None:
        if not 0 <= probability <= 1:
            raise Exception("Crossover probability value invalid: " + str(probability))
        self.remove_gap_columns = remove_gap_columns
        self.has_solution_been_crossed = None
        super(GapSequenceSolutionSinglePoint, self).__init__(probability=probability)

    def execute(self, parents: List[MSASolution]) -> List[MSASolution]:
        if len(parents) != 2:
            raise Exception("The number of parents is not two (2) but " + str(len(parents)))

        return self.do_crossover(parents)

    def do_crossover(self, parents: List[MSASolution]) -> List[MSASolution]:
        if random.random() <= self.probability:
            cx_point = random.randint(0, parents[0].get_length_of_alignment() - 1)

            column_positions_in_first_parent = \
                self.__find_original_positions_in_original_sequences(parents[0], cx_point)
            column_positions_in_second_parent = \
                self.__find_original_positions_in_aligned_sequence(parents[1], column_positions_in_first_parent)

            cutting_points_in_first_parent =\
                self.__find_cutting_points_in_first_parent(parents[0], cx_point)

            offspring = self.cross_parents(parents, cutting_points_in_first_parent, column_positions_in_second_parent)

            #if self.remove_gap_columns:
            #    offspring[0].remove_full_of_gaps_columns()
            #    offspring[1].remove_full_of_gaps_columns()

            self.has_solution_been_crossed = True
        else:
            offspring = [copy.deepcopy(parents[0]), copy.deepcopy(parents[1])]
            self.has_solution_been_crossed = False

        return offspring

    def cross_parents(self, parents: List[MSASolution], cutting_points_in_first_parent: list,
                      column_positions_in_second_parent: list) -> List[MSASolution]:
        offspring_1 = copy.deepcopy(parents[0])
        offspring_2 = copy.deepcopy(parents[1])

        # Obtain first children
        for i in range(0, offspring_1.number_of_variables):
            parent_1_gaps_groups = parents[0].gaps_groups
            parent_2_gaps_groups = parents[1].gaps_groups

            if cutting_points_in_first_parent[i] != -1:
                # offspring 1: Add the gap groups of the first parent before the cutting point
                # offspring 1: Add the gap groups on the first parent after the cutting point
                # offspring 2: Add the gap groups of the second parent before the cutting point
                # offspring 2: Add the gap groups on the second parent after the cutting point
                pass

        max_sequence_length = self.__find_length_of_the_largest_sequence(offspring_1)
        self.__fill_sequences_with_gaps_to_reach_the_max_sequence_length(offspring_1, max_sequence_length,
                                                                         cutting_points_in_first_parent)
        self.__fill_sequences_with_gaps_to_reach_the_max_sequence_length(offspring_2, max_sequence_length,
                                                                         column_positions_in_second_parent)

        return [offspring_1, offspring_2]

    def __find_original_positions_in_original_sequences(self, solution: MSASolution, column: int) -> list:
        positions = []

        for i in range(solution.number_of_variables):
            positions[i] = self.__find_symbol_position_in_original_sequence(solution, i, column)

        return positions

    def __find_symbol_position_in_original_sequence(self, solution: MSASolution, seq_index: int, position: int):
        if position > solution.get_length_of_alignment():
            raise Exception("Position {0} is larger than the sequence size {1}".format(position, solution.get_length_of_alignment()))

        if not solution.is_gap_char_at_sequence(seq_index, position):
            symbol_position = solution.get_char_position_in_original_sequence(seq_index, position)
        else:
            position = solution.get_next_char_position_after_gap(seq_index, position)
            if position < 0:
                symbol_position = -1
            else:
                symbol_position = solution.get_char_position_in_original_sequence(seq_index, position)

        return symbol_position

    def __find_original_positions_in_aligned_sequence(self,  solution: MSASolution, column_positions_in_first_parent: list):
        positions = []

        for seq_index in range(0, solution.number_of_variables):
            pos = column_positions_in_first_parent[seq_index]
            positions[seq_index] = solution.get_original_char_position_in_aligned_sequence(seq_index, pos)

        return positions

    def __find_cutting_points_in_first_parent(self, solution: MSASolution, position: int) -> list:
        positions = []

        for seq_index in range(0, solution.number_of_variables):
            if solution.is_gap_char_at_sequence(seq_index, position):
                positions[seq_index] = solution.get_next_char_position_after_gap(seq_index, position)
            else:
                positions[seq_index] = position

        return positions

    def __find_length_of_the_largest_sequence(self, solution: MSASolution):
        max_length = solution.get_length_of_alignment()

        for i in range(0, solution.number_of_variables):
            if max_length < solution.get_length_of_alignment():
                max_length = solution.get_length_of_alignment()

        return max_length

    def __fill_sequences_with_gaps_to_reach_the_max_sequence_length(self, solution: MSASolution, max_length: int,
                                                                    positions: list):
        pass

    def get_number_of_parents(self) -> int:
        return 2
