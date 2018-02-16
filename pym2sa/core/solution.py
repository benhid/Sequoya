from jmetal.core.solution import Solution


class MSASolution(Solution[str]):
    """ Class representing MSA solutions """

    def __init__(self, aligned_sequences: list, number_of_objectives: int) -> None:
        super(MSASolution, self).__init__(number_of_variables=len(aligned_sequences),
                                          number_of_objectives=number_of_objectives)
        self.GAP_IDENTIFIER = '-'

        self.sequences_names = list(pair[0] for pair in aligned_sequences)
        self.original_sequences = list(pair[1] for pair in aligned_sequences)
        self.original_alignment_size = len(self.original_sequences[0])
        self.gaps_groups = []

        self.encode_alignment(self.original_sequences)

    def encode_alignment(self, aligned_sequences: list):
        # for each aligned sequence
        for i in range(len(aligned_sequences)):
            # get gaps groups
            self.gaps_groups.append(self.get_gaps_group(aligned_sequences[i]))
            # get encoded sequence
            self.variables[i] = aligned_sequences[i].replace(self.GAP_IDENTIFIER, "")

    def get_gaps_group(self, sequence: str) -> list:
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

    def decode_alignment(self):
        aligned_sequences = []

        for i in range(self.number_of_variables):
            aligned_sequences.append(self.decode(self.variables[i],
                                                 self.gaps_groups[i]))

        return aligned_sequences

    def decode(self, encoded_sequence: str, gaps_group: list) -> str:
        aligned_sequence = list(encoded_sequence)

        # insert gap groups
        for i in range(0, len(gaps_group) - 1, 2):
            for j in range(gaps_group[i], gaps_group[i + 1] + 1):
                aligned_sequence.insert(j, '-')

        return ''.join(aligned_sequence)

    def decode_alignment_as_list_of_pairs(self) -> list:
        list_of_pairs = []

        for i in range(self.number_of_variables):
            sequence = self.decode(self.variables[i], self.gaps_groups[i])
            list_of_pairs.append((self.sequences_names[i], sequence))

        return list_of_pairs

    def merge_gaps_groups(self) -> None:
        # e.g. [[2, 4, 4, 7], [2, 3], [5, 6, 6, 8]] = [[2, 7], [2, 3], [5, 8]]

        for i in range(self.number_of_variables):
            gaps_group = self.gaps_groups[i]

            # for each gap on the gaps group
            for j in range(0, len(gaps_group) - 1, 2):
                if gaps_group[i][j] == gaps_group[i][j+1]:
                    gaps_group.remove(j)
                    gaps_group.remove(j+1)

    def split_gap_column(self, column_index: int) -> None:
        # for each sequence
        for i in range(self.number_of_variables):
            # get gaps group
            gaps_group = self.gaps_groups[i]
            # for each gap on the gaps group
            for j in range(0, len(gaps_group) - 1, 2):
                if gaps_group[j] <= column_index < gaps_group[j + 1]:
                    gaps_group.insert(j + 1, column_index + 1)
                    gaps_group.insert(j + 1, column_index)

    def remove_gap_column(self, column_index: int) -> None:
        for i in range(self.number_of_variables):
            gaps_group = self.gaps_groups[i]

            for j in range(0, len(self.gaps_groups) - 1, 2):
                if gaps_group[j] <= column_index <= gaps_group[j + 1]:
                    if gaps_group[j] == column_index:
                        # e.g. index: 3, gaps_group: [[2, 3]] -> [[2, 2]]
                        gaps_group[j] += 1
                    elif gaps_group[j + 1] == column_index:
                        # e.g. index: 3, gaps_group: [[3, 5]] -> [[4, 5]]
                        gaps_group[j + 1] -= 1
                    else:
                        # e.g. index: 3, gaps_group: [[2, 5]] -> [[2, 2, 4, 5]]
                        gaps_group.insert(j + 1, column_index + 1)
                        gaps_group.insert(j + 1, column_index - 1)

    def is_gap_column(self, index: int) -> bool:
        # check if the column index is in all gaps groups

        for i in range(self.number_of_variables):
            # for each sequence, get gaps group
            gaps_group = self.gaps_groups[i]
            gap_column = False

            # check if the column index in in its gaps group
            for j in range(0, len(gaps_group) - 1, 2):
                if gaps_group[j] <= index <= gaps_group[j + 1]:
                    gap_column = True
                    break

            # if gap_column is false for this sequence, there's no need to keep searching on the rest of them
            if not gap_column:
                return False

        return True

    def get_gap_columns(self) -> list:
        gap_columns_index = []
        alignment_length = self.get_length_of_alignment()

        for j in range(alignment_length):
            if self.is_gap_column(j):
                if j not in gap_columns_index:
                    gap_columns_index.append(j)
            else:
                if j in gap_columns_index:
                    gap_columns_index.remove(j)

        return gap_columns_index

    def get_total_number_of_gaps(self) -> int:
        number_of_gaps = 0

        for i in range(self.number_of_variables):
            number_of_gaps += self.get_number_of_gaps_of_sequence(i)

        return number_of_gaps

    def get_number_of_gaps_of_sequence(self, seq_index: int):
        number_of_gaps = 0
        gaps_group = self.gaps_groups[seq_index]

        for i in range(0, len(gaps_group) - 1, 2):
            number_of_gaps += gaps_group[i + 1] - gaps_group[i] + 1

        return number_of_gaps

    def get_length_of_alignment(self) -> int:
        return self.original_alignment_size

    def is_valid(self) -> bool:
        # a solution is valid if the length of the alignment is equals to the length of the encoded sequences (
        # without gaps) + the number of gaps
        alignment_length = self.get_length_of_alignment()

        for i in range(self.number_of_variables):
            if alignment_length != len(self.variables[i]) + self.get_number_of_gaps_of_sequence(i):
                return False

        return True
