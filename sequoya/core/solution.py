from jmetal.core.solution import Solution


class MSASolution(Solution[str]):
    """
    Class representing MSA solutions.
    """

    GAP_IDENTIFIER = '-'

    def __init__(self, problem, msa: list) -> None:
        super(MSASolution, self).__init__(number_of_variables=problem.number_of_variables,
                                          number_of_objectives=problem.number_of_objectives)

        self.sequences_names = problem.identifiers
        self.gaps_groups = [[] for _ in range(self.number_of_variables)]

        self.encode_alignment(list(pair[1] for pair in msa))

    def encode_alignment(self, aligned_sequences: list):
        # for each bb3_aligned sequence
        for index, seq in enumerate(aligned_sequences):
            # get gaps groups
            self.gaps_groups[index] = self.__get_gaps_group_of_sequence(seq)
            # get encoded sequence
            self.variables[index] = seq.replace(self.GAP_IDENTIFIER, "")

    def decode_sequence_at_index(self, seq_index: int):
        return self.__decode(self.variables[seq_index], self.gaps_groups[seq_index])

    def decode_alignment_as_list_of_sequences(self) -> list:
        aligned_sequences = []

        for i in range(self.number_of_variables):
            aligned_sequences.append(self.__decode(self.variables[i], self.gaps_groups[i]))

        return aligned_sequences

    def decode_alignment_as_list_of_pairs(self) -> list:
        list_of_pairs = []

        for i in range(self.number_of_variables):
            sequence = self.__decode(self.variables[i], self.gaps_groups[i])
            list_of_pairs.append((self.sequences_names[i], sequence))

        return list_of_pairs

    def __decode(self, encoded_sequence: str, gaps_group: list) -> str:
        aligned_sequence = list(encoded_sequence)

        # insert gap groups
        for i in range(0, len(gaps_group) - 1, 2):
            for j in range(gaps_group[i], gaps_group[i + 1] + 1):
                aligned_sequence.insert(j, '-')

        return ''.join(aligned_sequence)

    def merge_gaps_groups(self) -> None:
        for i in range(self.number_of_variables):
            gaps_group = self.gaps_groups[i]

            # for each gap on the gaps group
            j = 1
            while j <= len(gaps_group) - 2:
                if gaps_group[j] == gaps_group[j + 1]:
                    gaps_group[j] = gaps_group[j + 2]
                    del gaps_group[j + 1]
                    del gaps_group[j + 1]
                    j -= 2
                elif gaps_group[j] + 1 == gaps_group[j + 1]:
                    gaps_group[j] = gaps_group[j + 2]
                    del gaps_group[j + 1]
                    del gaps_group[j + 1]
                    j -= 2
                j += 2

    def add_gap_to_sequence_at_index(self, seq_index: int, gap_position: int):
        new_gaps_group = self.gaps_groups[seq_index]
        length_of_alignment = self.get_length_of_alignment()

        if gap_position == length_of_alignment:
            if self.is_gap_char_at_sequence(seq_index, gap_position - 1):
                new_gaps_group[gap_position - 1] += 1
            else:
                new_gaps_group.append(gap_position)
                new_gaps_group.append(gap_position)
        elif gap_position >= length_of_alignment:
            new_gaps_group.append(length_of_alignment)
            new_gaps_group.append(length_of_alignment)
        elif self.is_gap_char_at_sequence(seq_index, gap_position):
            # increments gaps group size
            gap_added = False
            for j in range(0, len(new_gaps_group) - 1, 2):
                if gap_added:
                    new_gaps_group[j] += 1
                    new_gaps_group[j + 1] += 1

                if new_gaps_group[j] == gap_position or new_gaps_group[j + 1] == gap_position:
                    new_gaps_group[j + 1] += 1
                    gap_added = True
                elif new_gaps_group[j] < gap_position < new_gaps_group[j + 1]:
                    new_gaps_group[j + 1] += 1
                    gap_added = True
        else:
            # add new group at position
            for j in range(0, len(new_gaps_group) - 1, 2):
                if new_gaps_group[j] > gap_position:
                    new_gaps_group[j] += 1
                    new_gaps_group[j + 1] += 1

            new_gaps_group.append(gap_position)
            new_gaps_group.append(gap_position)
            new_gaps_group.sort()

        self.gaps_groups[seq_index] = new_gaps_group

    def split_gap_column(self, column: int) -> None:
        # for each sequence
        for i in range(self.number_of_variables):
            # get gaps group
            gaps_group = self.gaps_groups[i]
            # for each gap on the gaps group
            for j in range(0, len(gaps_group) - 1, 2):
                if gaps_group[j] <= column < gaps_group[j + 1]:
                    gaps_group.insert(j + 1, column + 1)
                    gaps_group.insert(j + 1, column)

    def remove_full_of_gaps_columns(self) -> None:
        gap_columns = self.get_gap_columns_from_alignment()
        gap_columns.sort(reverse=True)

        for col in gap_columns:
            for seq_index in range(self.number_of_variables):
                self.remove_gap_from_sequence(seq_index, col)

    def remove_gap_column(self, column: int) -> None:
        if not self.is_gap_column(column):
            raise Exception("No gap group in position {0}".format(column))
        else:
            for i in range(self.number_of_variables):
                gaps_group = self.gaps_groups[i]

                for j in range(0, len(self.gaps_groups) - 1, 2):
                    if gaps_group[j] == column:
                        # e.g. index: 3, gaps_group: [[2, 3]] -> [[2, 2]]
                        gaps_group[j] += 1
                    elif gaps_group[j + 1] == column:
                        # e.g. index: 3, gaps_group: [[3, 5]] -> [[4, 5]]
                        gaps_group[j + 1] -= 1
                    else:
                        # e.g. index: 3, gaps_group: [[2, 5]] -> [[2, 2, 4, 5]]
                        gaps_group.insert(j + 1, column + 1)
                        gaps_group.insert(j + 1, column - 1)

    def remove_gap_group_from_sequence_at_column(self, seq_index: int, column_index: int) -> None:
        if not self.is_gap_char_at_sequence(seq_index, column_index):
            raise Exception("No gap group in position {0} at sequence {1}".format(column_index, seq_index))
        else:
            gaps_group = self.gaps_groups[seq_index]

            for j in range(0, len(gaps_group) - 1, 2):
                if gaps_group[j] == column_index or gaps_group[j + 1] == column_index \
                        or gaps_group[j] < column_index < gaps_group[j + 1]:
                    del gaps_group[j]
                    del gaps_group[j]

                    break

    def remove_gap_from_sequence(self, seq_index: int, position: int):
        new_gaps_group = self.gaps_groups[seq_index]

        if self.is_gap_char_at_sequence(seq_index, position):
            for j in range(0, len(new_gaps_group) - 1, 2):
                if new_gaps_group[j] == position or new_gaps_group[j + 1] == position:
                    if new_gaps_group[j] == new_gaps_group[j + 1]:
                        new_gaps_group[j + 2:] = [x - 1 for x in new_gaps_group[j + 2:]]
                        del new_gaps_group[j]
                        del new_gaps_group[j]
                        break
                    else:
                        new_gaps_group[j + 1] -= 1
                        new_gaps_group[j + 2:] = [x - 1 for x in new_gaps_group[j + 2:]]
                        break
                elif new_gaps_group[j] < position < new_gaps_group[j + 1]:
                    new_gaps_group[j + 1] -= 1
                    new_gaps_group[j + 2:] = [x - 1 for x in new_gaps_group[j + 2:]]
                    break

        self.gaps_groups[seq_index] = new_gaps_group

    def is_gap_column(self, column: int) -> bool:
        # check if the column index is in all gaps groups

        for i in range(self.number_of_variables):
            # for each sequence, get gaps group
            gaps_group = self.gaps_groups[i]
            gap_column = False

            # check if the column index in in its gaps group
            for j in range(0, len(gaps_group) - 1, 2):
                if gaps_group[j] <= column <= gaps_group[j + 1]:
                    gap_column = True
                    break

            # if gap_column is false for this sequence, there's no need to keep searching on the rest of them
            if not gap_column:
                return False

        return True

    def is_gap_char_at_sequence(self, seq_index: int, index: int) -> bool:
        assert seq_index <= self.number_of_variables - 1, "Sequence doesn't exist on this alignment"

        if index > self.get_length_of_sequence(seq_index) or index < 0:
            raise Exception(
                "Index out of sequence: index {0}, alignment length: {1}, sequence length: {2}".format(
                    index, self.get_length_of_alignment(), self.get_length_of_sequence(seq_index)))

        gaps_group = self.gaps_groups[seq_index]

        for a, b in zip(*[iter(gaps_group)] * 2):
            if a <= index <= b:
                return True

        return False

    def get_char_position_in_original_sequence(self, seq_index: int, position: int):
        gaps_group = self.gaps_groups[seq_index]
        gaps_on_the_left = []

        for j in range(0, len(gaps_group) - 1, 2):
            if gaps_group[j + 1] <= position:
                gaps_on_the_left.append(gaps_group[j])
                gaps_on_the_left.append(gaps_group[j + 1])

        if len(gaps_on_the_left) > 0:
            length_of_gaps = \
                sum([gaps_on_the_left[j + 1] - gaps_on_the_left[j] + 1 for j in range(0, len(gaps_on_the_left) - 1, 2)])

            position -= length_of_gaps

        return position

    def get_next_char_position_after_gap(self, seq_index: int, gap_position: int):
        if not self.is_gap_char_at_sequence(seq_index, gap_position):
            raise Exception("Symbol in position {0} is not a gap!".format(gap_position))

        position = -1
        gaps_group = self.gaps_groups[seq_index]

        for j in range(0, len(gaps_group) - 1, 2):
            if gaps_group[j] <= gap_position <= gaps_group[j + 1]:
                position = gaps_group[j + 1] + 1

        if position == self.get_length_of_sequence(seq_index):
            position = -1

        return position

    def get_original_char_position_in_aligned_sequence(self, seq_index: int, position: int):
        symbol_position = 0
        found_symbols = 0

        if position < 0:
            symbol_position = -1
        else:
            found = False
            while symbol_position < self.get_length_of_alignment() and not found:
                if self.is_gap_char_at_sequence(seq_index, symbol_position):
                    symbol_position += 1
                else:
                    if found_symbols == position:
                        found = True
                    else:
                        found_symbols += 1
                        symbol_position += 1

        return symbol_position

    def get_length_of_gaps(self, seq_index: int) -> int:
        if not self.gaps_groups[seq_index]:
            length_of_gaps = 0
        else:
            length_of_gaps = 0

            gaps_group = self.gaps_groups[seq_index]

            for a, b in zip(*[iter(gaps_group)] * 2):
                length_of_gaps += b - a + 1

        return length_of_gaps

    def __get_gaps_group_of_sequence(self, sequence: str) -> list:
        gaps_group = []
        gap_open = False
        start = 0

        for i in range(len(sequence)):
            if sequence[i] is '-':
                if not gap_open:
                    gap_open = True
                    start = i
            else:
                if gap_open:
                    gap_open = False

                    # e.g. [1, 2, 5, 6] means there are two gaps groups:
                    #   the first one from (1,2) and the second one from (5,6)
                    gaps_group.append(start)
                    gaps_group.append(i - 1)

        if gap_open:
            gaps_group.append(start)
            gaps_group.append(len(sequence) - 1)

        return gaps_group

    def get_number_of_gaps_groups_of_sequence(self, seq_index: int) -> float:
        return len(self.gaps_groups[seq_index]) / 2

    def get_gap_columns_from_alignment(self) -> list:
        gap_columns_index = []
        alignment_length = self.get_length_of_alignment()

        for j in range(alignment_length):
            if self.is_gap_column(j):
                gap_columns_index.append(j)

        return gap_columns_index

    def get_total_number_of_gaps(self) -> int:
        number_of_gaps = 0

        for i in range(self.number_of_variables):
            number_of_gaps += self.get_number_of_gaps_of_sequence_at_index(i)

        return number_of_gaps

    def get_number_of_gaps_of_sequence_at_index(self, seq_index: int):
        number_of_gaps = 0
        gaps_group = self.gaps_groups[seq_index]

        for i in range(0, len(gaps_group) - 1, 2):
            number_of_gaps += gaps_group[i + 1] - gaps_group[i] + 1

        return number_of_gaps

    def get_length_of_alignment(self) -> int:
        return len(self.variables[0]) + self.get_length_of_gaps(0)

    def get_length_of_sequence(self, seq_index: int) -> int:
        return len(self.variables[seq_index]) + self.get_length_of_gaps(seq_index)

    def is_valid_msa(self) -> bool:
        """ Check if all sequences have the same length """
        if not all(self.get_length_of_sequence(seq_index) == self.get_length_of_sequence(0)
                   for seq_index in range(1, self.number_of_variables)):
            return False
        return True

    def __str__(self):
        fasta = ""
        n = 200

        for index, seq in enumerate(self.sequences_names):
            fasta += '>' + seq + '\n'
            sequence = self.decode_sequence_at_index(index)
            sequence = [sequence[i:i + n] for i in range(0, len(sequence), n)]

            for wrap in sequence:
                fasta += wrap + "\n"

        fasta += "\n"

        return fasta
