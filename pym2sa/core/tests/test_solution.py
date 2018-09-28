import unittest

from pym2sa.core.solution import MSASolution
from pym2sa.problem import MSA


class MSASolutionTestCases(unittest.TestCase):

    def test_should_return_original_sequences(self):
        # setup
        problem = MSA(score_list=[], original_sequences=['ACTGAC', 'ATCTC', 'AACTGC'], sequences_names=['seq1', 'seq2', 'seq3'])
        msa = MSASolution(problem, msa=[('seq1', 'AC---TGAC'), ('seq2', 'AT--CT--C'), ('seq3', 'AAC---TGC')])

        # check
        self.assertEqual(['AC---TGAC', 'AT--CT--C', 'AAC---TGC'], msa.decode_alignment_as_list_of_sequences())

    def test_should_return_original_alignment_size(self):
        # setup
        problem = MSA(score_list=[], original_sequences=['ACTGAC', 'ATCTC', 'AACTGC'], sequences_names=['seq1', 'seq2', 'seq3'])
        msa = MSASolution(problem, msa=[('seq1', 'AC---TGAC'), ('seq2', 'AT--CT--C'), ('seq3', 'AAC---TGC')])

        # check
        self.assertEqual(9, msa.get_length_of_alignment())

    def test_should_return_gaps_groups(self):
        # setup
        problem = MSA(score_list=[], original_sequences=['ACTGAC', 'ATCTC', 'AACTGC'], sequences_names=['seq1', 'seq2', 'seq3'])
        msa = MSASolution(problem, msa=[('seq1', 'AC---TGAC'), ('seq2', 'AT--CT--C'), ('seq3', 'AAC---TGC')])

        # check
        self.assertEqual([[2, 4], [2, 3, 6, 7], [3, 5]], msa.gaps_groups)

    def test_should_return_length_of_gaps_groups(self):
        # setup
        problem = MSA(score_list=[], original_sequences=['ACTGAC', 'ATCTC', 'AACTGC'], sequences_names=['seq1', 'seq2', 'seq3'])
        msa_1 = MSASolution(problem, msa=[('seq1', 'AC---TGAC'), ('seq2', 'AT--CT--C'), ('seq3', 'AAC---TGC')])

        problem = MSA(score_list=[], original_sequences=['GKGDPKKP', 'GKGDPKKP', 'MKKLKKHPDFP', 'MHI'], sequences_names=['seq1', 'seq2', 'seq3', 'seq4'])
        msa_2 = MSASolution(problem, msa=[('seq1', 'GKGD---PKKP'), ('seq2', 'GKGD---PKKP'),
                                          ('seq3', 'MKKLKKHPDFP'), ('seq4', 'M--------HI-')])

        # check
        self.assertEqual(3, msa_1.get_length_of_gaps(0))
        self.assertEqual(4, msa_1.get_length_of_gaps(1))
        self.assertEqual(3, msa_1.get_length_of_gaps(2))

        self.assertEqual(3, msa_2.get_length_of_gaps(0))
        self.assertEqual(3, msa_2.get_length_of_gaps(1))
        self.assertEqual(0, msa_2.get_length_of_gaps(2))
        self.assertEqual(9, msa_2.get_length_of_gaps(3))

    def test_should_return_number_of_gaps_of_all_sequences(self):
        # setup
        problem = MSA(score_list=[], original_sequences=['ACTGAC', 'ATCTC', 'AACTGC'], sequences_names=['seq1', 'seq2', 'seq3'])
        msa = MSASolution(problem, msa=[('seq1', 'AC---TGAC'), ('seq2', 'AT--CT--C'), ('seq3', 'AAC---TGC')])

        # check
        self.assertEqual(10, msa.get_total_number_of_gaps())

    def test_should_return_number_of_gaps_of_one_sequences(self):
        # setup
        problem = MSA(score_list=[], original_sequences=['ACTGAC', 'ATCTC', 'AACTGC'], sequences_names=['seq1', 'seq2', 'seq3'])
        msa = MSASolution(problem, msa=[('seq1', 'AC---TGAC'), ('seq2', 'AT--CT--C'), ('seq3', 'AAC---TGC')])

        # check
        self.assertEqual(3, msa.get_number_of_gaps_of_sequence_at_index(0))

    def test_should_return_if_msa_is_valid(self):
        # setup
        problem = MSA(score_list=[], original_sequences=['ACTGAC', 'ATCTC', 'AACTGC'], sequences_names=['seq1', 'seq2', 'seq3'])
        msa_valid = MSASolution(problem, msa=[('seq1', 'AC---TGAC'), ('seq2', 'AT--CT--C'), ('seq3', 'AAC---TGC')])
        msa_not_valid = MSASolution(problem, msa=[('seq1', 'A'), ('seq2', 'A'), ('seq3', 'AA')])

        # check
        self.assertTrue(msa_valid.is_valid_msa())
        self.assertFalse(msa_not_valid.is_valid_msa())

    def test_should_return_alignment_as_list_of_pairs(self):
        # setup
        problem = MSA(score_list=[], original_sequences=['ACTGAC', 'ATCTC', 'AACTGC'], sequences_names=['seq1', 'seq2', 'seq3'])
        msa = MSASolution(problem, msa=[('seq1', 'AC---TGAC'), ('seq2', 'AT--CT--C'), ('seq3', 'AAC---TGC')])

        # check
        self.assertEqual([('seq1', 'AC---TGAC'), ('seq2', 'AT--CT--C'), ('seq3', 'AAC---TGC')],
                         msa.decode_alignment_as_list_of_pairs())

    def test_should_merge_gaps_groups(self):
        # setup
        problem = MSA(score_list=[], original_sequences=['ACTGAC', 'ATCTC'], sequences_names=['seq1', 'seq2'])
        msa = MSASolution(problem, msa=[('seq1', 'ACTGAC'), ('seq2', 'ATCTC')])
        msa.gaps_groups[0] = [2, 4, 4, 5]
        msa.gaps_groups[1] = [2, 4, 5, 8]

        self.assertEqual(["AC-----TGAC", "AT-------CTC"], msa.decode_alignment_as_list_of_sequences())

        # run
        msa.merge_gaps_groups()

        # check
        self.assertEqual([2, 5], msa.gaps_groups[0])
        self.assertEqual([2, 8], msa.gaps_groups[1])

    def test_should_merge_gaps_groups_case_b(self):
        # setup
        problem = MSA(score_list=[], original_sequences=['ACTGAC'], sequences_names=['seq1'])
        msa = MSASolution(problem, msa=[('seq1', 'ACTGAC')])
        msa.gaps_groups[0] = [2, 4, 4, 8, 8, 10]

        self.assertEqual(["AC-----------TGAC"], msa.decode_alignment_as_list_of_sequences())

        # run
        msa.merge_gaps_groups()

        # check
        self.assertEqual([2, 10], msa.gaps_groups[0])

    def test_should_merge_gaps_groups_case_c(self):
        # setup
        problem = MSA(score_list=[], original_sequences=['ACTGAC'], sequences_names=['seq1'])
        msa = MSASolution(problem, msa=[('seq1', 'ACTGAC')])
        msa.gaps_groups[0] = [2, 4, 4, 8, 8, 10, 10, 10]

        self.assertEqual(["AC------------TGAC"], msa.decode_alignment_as_list_of_sequences())

        # run
        msa.merge_gaps_groups()

        # check
        self.assertEqual([2, 10], msa.gaps_groups[0])

    def test_should_return_is_gap_column(self):
        # setup
        problem = MSA(score_list=[], original_sequences=['ACTGAC', 'ATCTC', 'AACTGC'], sequences_names=['seq1', 'seq2', 'seq3'])
        msa = MSASolution(problem, msa=[('seq1', 'AC---TGAC'), ('seq2', 'AT--CT--C'), ('seq3', 'AAC---TGC')])

        # check
        self.assertTrue(msa.is_gap_column(3))
        self.assertFalse(msa.is_gap_column(4))

    def test_should_split_gap_column(self):
        # setup
        problem = MSA(score_list=[], original_sequences=['AC', 'TC', 'AA'], sequences_names=['seq1', 'seq2', 'seq3'])
        msa_1 = MSASolution(problem, msa=[('seq1', '----AC'), ('seq2', 'T----C'), ('seq3', '--A-A-')])
        msa_2 = MSASolution(problem, msa=[('seq1', '----AC'), ('seq2', 'T----C'), ('seq3', '--A-A-')])
        msa_3 = MSASolution(problem, msa=[('seq1', '----AC'), ('seq2', 'T----C'), ('seq3', '-A---A')])

        # actual = [[0, 2], [1, 3], [0, 1, 3, 3, 5, 5]]
        msa_1.split_gap_column(1)
        msa_2.split_gap_column(2)

        # actual = [[0, 2], [1, 3], [0, 0, 2, 4]]
        msa_3.split_gap_column(3)

        # check
        self.assertEqual([[0, 1, 2, 3], [1, 1, 2, 4], [0, 1, 3, 3, 5, 5]], msa_1.gaps_groups)
        self.assertEqual([[0, 2, 3, 3], [1, 2, 3, 4], [0, 1, 3, 3, 5, 5]], msa_2.gaps_groups)
        self.assertEqual([[0, 3], [1, 3, 4, 4], [0, 0, 2, 3, 4, 4]], msa_3.gaps_groups)

    def test_should_remove_gap_column(self):
        # setup
        problem = MSA(score_list=[], original_sequences=['ACTGAC', 'ATCTC', 'AACTGC'], sequences_names=['seq1', 'seq2', 'seq3'])
        msa = MSASolution(problem, msa=[('seq1', 'AC---TGAC'), ('seq2', 'AT--CT--C'), ('seq3', 'AAC---TGC')])

        msa.remove_gap_column(3)

        # check
        self.assertEqual([[2, 2, 4, 4], [2, 2, 6, 7], [4, 5]], msa.gaps_groups)

    def test_should_remove_gap(self):
        # setup
        problem = MSA(score_list=[], original_sequences=['ACTGAC', 'ACTGAC', 'ACTGAC'], sequences_names=['seq1', 'seq2', 'seq3'])
        msa = MSASolution(problem, msa=[('seq1', 'AC---TGAC'), ('seq2', 'AC---TGAC'), ('seq3', 'AC---TGAC')])

        msa.remove_gap_from_sequence(0, 2)
        msa.remove_gap_from_sequence(1, 2)
        msa.remove_gap_from_sequence(2, 2)

        # check
        self.assertEqual(['AC--TGAC', 'AC--TGAC', 'AC--TGAC'], msa.decode_alignment_as_list_of_sequences())

    def test_should_remove_gap_case_b(self):
        # setup
        problem = MSA(score_list=[], original_sequences=['ABCDE', 'ABCDE'], sequences_names=['seq1', 'seq2'])
        msa = MSASolution(problem, msa=[('seq1', 'AB--CD-E'), ('seq2', 'AB--CDE-')])

        msa.remove_gap_from_sequence(0, 3)
        self.assertEqual(['AB-CD-E', 'AB--CDE-'], msa.decode_alignment_as_list_of_sequences())

        msa.remove_gap_from_sequence(0, 2)
        self.assertEqual(['ABCD-E', 'AB--CDE-'], msa.decode_alignment_as_list_of_sequences())

        msa.remove_gap_from_sequence(1, 3)
        self.assertEqual(['ABCD-E', 'AB-CDE-'], msa.decode_alignment_as_list_of_sequences())

        msa.remove_gap_from_sequence(1, 2)
        self.assertEqual(['ABCD-E', 'ABCDE-'], msa.decode_alignment_as_list_of_sequences())

    def test_should_remove_gap_case_c(self):
        # setup
        problem = MSA(score_list=[], original_sequences=['ABCDE'], sequences_names=['seq1'])
        msa = MSASolution(problem, msa=[('seq1', 'AB----CD-E-')])

        msa.remove_gap_from_sequence(0, 3)
        self.assertEqual(['AB---CD-E-'], msa.decode_alignment_as_list_of_sequences())

        msa.remove_gap_from_sequence(0, 3)
        self.assertEqual(['AB--CD-E-'], msa.decode_alignment_as_list_of_sequences())

        msa.remove_gap_from_sequence(0, 8)
        self.assertEqual(['AB--CD-E'], msa.decode_alignment_as_list_of_sequences())

    def test_should_remove_all_gap_columns(self):
        # setup
        problem = MSA(score_list=[], original_sequences=['ACTGAC', 'ACTGAC', 'ACTGAC'], sequences_names=['seq1', 'seq2', 'seq3'])
        msa = MSASolution(problem, msa=[('seq1', 'AC---TGAC'), ('seq2', 'AC---TGAC'), ('seq3', 'AC---TGAC')])

        msa.remove_full_of_gaps_columns()

        # check
        self.assertEqual(['ACTGAC', 'ACTGAC', 'ACTGAC'], msa.decode_alignment_as_list_of_sequences())

    def test_should_remove_all_gap_columns_case_b(self):
        # setup
        problem = MSA(score_list=[], original_sequences=['ACTGC', 'ACAC', 'ACAC'], sequences_names=['seq1', 'seq2', 'seq3'])
        msa = MSASolution(problem, msa=[('seq1', 'AC--T--GC'), ('seq2', 'AC-----AC'), ('seq3', 'A---C--AC')])

        msa.remove_full_of_gaps_columns()

        # check
        self.assertEqual(['ACTGC', 'AC-AC', 'A-CAC'], msa.decode_alignment_as_list_of_sequences())

    def test_should_remove_all_gap_columns_case_c(self):
        # setup
        problem = MSA(score_list=[], original_sequences=['', '', 'AA'], sequences_names=['seq1', 'seq2', 'seq3'])
        msa = MSASolution(problem, msa=[('seq1', '----'), ('seq2', '----'), ('seq3', '-AA-')])

        msa.remove_full_of_gaps_columns()

        # check
        self.assertEqual(['--', '--', 'AA'], msa.decode_alignment_as_list_of_sequences())

    def test_should_remove_all_gap_columns_case_d(self):
        # setup
        problem = MSA(score_list=[], original_sequences=['ABCDE', 'ABCDE'], sequences_names=['seq1', 'seq2'])
        msa = MSASolution(problem, msa=[('seq1', 'AB--CDE-'), ('seq2', 'AB--CD-E')])

        msa.remove_full_of_gaps_columns()

        # check
        self.assertEqual(['ABCDE-', 'ABCD-E'], msa.decode_alignment_as_list_of_sequences())

    def test_should_return_gap_columns(self):
        # setup
        problem = MSA(score_list=[], original_sequences=['AA', 'AA', 'AA'], sequences_names=['seq1', 'seq2', 'seq3'])
        msa = MSASolution(problem, msa=[('seq1', '--AA-'), ('seq2', '--AA-'), ('seq3', '--AA-')])

        # check
        self.assertEqual([0, 1, 4], msa.get_gap_columns_from_alignment())

    def test_should_is_gap_at_char_sequence_raise_exception_if_position_is_negative(self):
        # setup
        problem = MSA(score_list=[], original_sequences=['AA', 'AA', 'AA'], sequences_names=['seq1', 'seq2', 'seq3'])
        msa = MSASolution(problem, msa=[('seq1', '--AA-'), ('seq2', '--AA-'), ('seq3', '--AA-')])

        # check
        self.assertTrue(msa.is_gap_char_at_sequence(0, 0))
        self.assertTrue(msa.is_gap_char_at_sequence(0, 1))
        self.assertTrue(msa.is_gap_char_at_sequence(1, 4))
        self.assertTrue(msa.is_gap_char_at_sequence(2, 1))

        with self.assertRaises(Exception):
            msa.is_gap_char_at_sequence(0, -1)

    def test_remove_gap_group_at_column(self):
        # setup
        problem = MSA(score_list=[], original_sequences=['AA', 'AA', 'AA'], sequences_names=['seq1', 'seq2', 'seq3'])
        msa = MSASolution(problem, msa=[('seq1', '--AA-'), ('seq2', '---AA'), ('seq3', '--AA-')])

        msa.remove_gap_group_from_sequence_at_column(0, 0)
        msa.remove_gap_group_from_sequence_at_column(1, 0)
        msa.remove_gap_group_from_sequence_at_column(2, 4)

        # check
        self.assertEqual([[4, 4], [], [0, 1]], msa.gaps_groups)

    def test_should_get_the_right_char_position_in_the_original_sequence(self):
        # setup
        problem = MSA(score_list=[], original_sequences=['ABC'], sequences_names=['seq1'])
        msa_1 = MSASolution(problem, msa=[('seq1', '-ABC')])
        msa_2 = MSASolution(problem, msa=[('seq1', 'A--B-C')])

        # check
        self.assertEqual(0, msa_1.get_char_position_in_original_sequence(seq_index=0, position=1))
        self.assertEqual(1, msa_1.get_char_position_in_original_sequence(0, 2))
        self.assertEqual(2, msa_1.get_char_position_in_original_sequence(0, 3))

        self.assertEqual(0, msa_2.get_char_position_in_original_sequence(0, 0))
        self.assertEqual(1, msa_2.get_char_position_in_original_sequence(0, 3))
        self.assertEqual(2, msa_2.get_char_position_in_original_sequence(0, 5))

    def test_should_get_next_char_position_after_gap(self):
        # setup
        problem = MSA(score_list=[], original_sequences=['ABC', 'ABC'], sequences_names=['seq1', 'seq2'])
        msa_1 = MSASolution(problem, msa=[('seq1', '-ABC'), ('seq2', 'AB-C')])

        problem = MSA(score_list=[], original_sequences=['ABC'], sequences_names=['seq1'])
        msa_2 = MSASolution(problem, msa=[('seq1', 'A--BC')])
        msa_extra_gaps = MSASolution(problem, msa=[('seq1', '---A---BC------')])

        # check
        self.assertEqual(1, msa_1.get_next_char_position_after_gap(seq_index=0, gap_position=0))
        self.assertEqual(3, msa_1.get_next_char_position_after_gap(seq_index=1, gap_position=2))

        self.assertEqual(3, msa_2.get_next_char_position_after_gap(seq_index=0, gap_position=1))
        self.assertEqual(3, msa_2.get_next_char_position_after_gap(seq_index=0, gap_position=2))

    def test_should_get_original_char_position_in_aligned_sequence(self):
        # setup
        problem = MSA(score_list=[], original_sequences=['ABC', 'ABCD', 'AB'], sequences_names=['seq1', 'seq2', 'seq3'])
        msa = MSASolution(problem, msa=[('seq1', '-ABC'), ('seq2', 'ABCD'), ('seq3', '--AB')])

        # check
        self.assertEqual(1, msa.get_original_char_position_in_aligned_sequence(seq_index=0, position=0))
        self.assertEqual(2, msa.get_original_char_position_in_aligned_sequence(seq_index=0, position=1))
        self.assertEqual(3, msa.get_original_char_position_in_aligned_sequence(seq_index=0, position=2))

        self.assertEqual(0, msa.get_original_char_position_in_aligned_sequence(seq_index=1, position=0))
        self.assertEqual(1, msa.get_original_char_position_in_aligned_sequence(seq_index=1, position=1))
        self.assertEqual(2, msa.get_original_char_position_in_aligned_sequence(seq_index=1, position=2))

        self.assertEqual(2, msa.get_original_char_position_in_aligned_sequence(seq_index=2, position=0))
        self.assertEqual(3, msa.get_original_char_position_in_aligned_sequence(seq_index=2, position=1))

    def test_should_increments_gaps_group(self):
        # setup
        problem = MSA(score_list=[], original_sequences=['A'], sequences_names=['seq1'])
        msa_1 = MSASolution(problem, msa=[('seq1', 'A-')])
        msa_2 = MSASolution(problem, msa=[('seq1', '-A')])

        problem = MSA(score_list=[], original_sequences=['AC'], sequences_names=['seq1'])
        msa_3 = MSASolution(problem, msa=[('seq1', 'A-C')])
        msa_4 = MSASolution(problem, msa=[('seq1', 'A--C')])
        msa_5 = MSASolution(problem, msa=[('seq1', 'A--C-')])

        # run
        msa_1.add_gap_to_sequence_at_index(seq_index=0, gap_position=1)
        msa_2.add_gap_to_sequence_at_index(seq_index=0, gap_position=0)
        msa_3.add_gap_to_sequence_at_index(seq_index=0, gap_position=1)
        msa_4.add_gap_to_sequence_at_index(seq_index=0, gap_position=1)
        msa_5.add_gap_to_sequence_at_index(seq_index=0, gap_position=2)
        msa_5.add_gap_to_sequence_at_index(seq_index=0, gap_position=5)
        msa_5.add_gap_to_sequence_at_index(seq_index=0, gap_position=50)

        # check
        self.assertEqual([('seq1', 'A--')], msa_1.decode_alignment_as_list_of_pairs())
        self.assertEqual([('seq1', '--A')], msa_2.decode_alignment_as_list_of_pairs())
        self.assertEqual([('seq1', 'A--C')], msa_3.decode_alignment_as_list_of_pairs())
        self.assertEqual([('seq1', 'A---C')], msa_4.decode_alignment_as_list_of_pairs())
        self.assertEqual([('seq1', 'A---C---')], msa_5.decode_alignment_as_list_of_pairs())

    def test_should_create_new_gaps_group(self):
        # setup
        problem = MSA(score_list=[], original_sequences=['A'], sequences_names=['seq1'])
        msa_1 = MSASolution(problem, msa=[('seq1', 'A-')])
        msa_2 = MSASolution(problem, msa=[('seq1', '-A')])

        problem = MSA(score_list=[], original_sequences=['AC'], sequences_names=['seq1'])
        msa_3 = MSASolution(problem, msa=[('seq1', 'A-C')])

        problem = MSA(score_list=[], original_sequences=['AAA'], sequences_names=['seq1'])
        msa_4 = MSASolution(problem, msa=[('seq1', 'AAA')])

        # run
        msa_1.add_gap_to_sequence_at_index(seq_index=0, gap_position=0)
        msa_2.add_gap_to_sequence_at_index(seq_index=0, gap_position=2)
        msa_3.add_gap_to_sequence_at_index(seq_index=0, gap_position=3)
        msa_4.add_gap_to_sequence_at_index(seq_index=0, gap_position=1)

        # check
        self.assertEqual([('seq1', '-A-')], msa_1.decode_alignment_as_list_of_pairs())
        self.assertEqual([('seq1', '-A-')], msa_2.decode_alignment_as_list_of_pairs())
        self.assertEqual([('seq1', 'A-C-')], msa_3.decode_alignment_as_list_of_pairs())
        self.assertEqual([('seq1', 'A-AA')], msa_4.decode_alignment_as_list_of_pairs())


if __name__ == "__main__":
    unittest.main()
