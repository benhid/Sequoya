import unittest

from pym2sa.core.solution import MSASolution


class MSASolutionTestCases(unittest.TestCase):
    def setUp(self):
        print("setUp: RUNNING MSASolutionTestCases")

    def tearDown(self):
        print("tearDown: TEST ENDED")

    def test_should_return_original_sequences(self):
        # setup
        msa = MSASolution(aligned_sequences=[('seq1', 'AC---TGAC'), ('seq2', 'AT--CT--C'), ('seq3', 'AAC---TGC')],
                          number_of_objectives=2)

        # check
        self.assertEqual(['AC---TGAC', 'AT--CT--C', 'AAC---TGC'], msa.original_sequences)

    def test_should_return_original_alignment_size(self):
        # setup
        msa = MSASolution(aligned_sequences=[('seq1', 'AC---TGAC'), ('seq2', 'AT--CT--C'), ('seq3', 'AAC---TGC')],
                          number_of_objectives=2)

        # check
        self.assertEqual(9, msa.original_alignment_size)

    def test_should_return_gaps_groups(self):
        # setup
        msa = MSASolution(aligned_sequences=[('seq1', 'AC---TGAC'), ('seq2', 'AT--CT--C'), ('seq3', 'AAC---TGC')],
                          number_of_objectives=2)

        # check
        self.assertEqual([[2, 4], [2, 3, 6, 7], [3, 5]], msa.gaps_groups)

    def test_should_decode_sequences(self):
        # setup
        msa = MSASolution(aligned_sequences=[('seq1', 'AC---TGAC'), ('seq2', 'AT--CT--C'), ('seq3', 'AAC---TGC')],
                          number_of_objectives=2)

        # check
        self.assertEqual(['AC---TGAC', 'AT--CT--C', 'AAC---TGC'], msa.decode_alignment())

    def test_should_return_number_of_gaps_of_all_sequences(self):
        # setup
        msa = MSASolution(aligned_sequences=[('seq1', 'AC---TGAC'), ('seq2', 'AT--CT--C'), ('seq3', 'AAC---TGC')],
                          number_of_objectives=2)

        # check
        self.assertEqual(10, msa.get_total_number_of_gaps())

    def test_should_return_number_of_gaps_of_one_sequences(self):
        # setup
        msa = MSASolution(aligned_sequences=[('seq1', 'AC---TGAC'), ('seq2', 'AT--CT--C'), ('seq3', 'AAC---TGC')],
                          number_of_objectives=2)

        # check
        self.assertEqual(3, msa.get_number_of_gaps_of_sequence(0))

    def test_should_return_is_valid(self):
        # setup
        msa_valid = MSASolution(aligned_sequences=[('seq1', 'AC---TGAC'), ('seq2', 'AT--CT--C'), ('seq3', 'AAC---TGC')],
                                number_of_objectives=2)
        msa_not_valid = MSASolution(aligned_sequences=[('seq1', 'AC----TGAC'), ('seq2', 'AT--CT--C'), ('seq3', 'AAC---TGC')],
                                    number_of_objectives=2)

        # check
        self.assertTrue(msa_valid.is_valid())
        self.assertFalse(msa_not_valid.is_valid())

    def test_should_return_alignment_as_list_of_pairs(self):
        # setup
        aln_seq = [('seq1', 'AC---TGAC'), ('seq2', 'AT--CT--C'), ('seq3', 'AAC---TGC')]
        msa = MSASolution(aligned_sequences=aln_seq,
                          number_of_objectives=2)

        # check
        self.assertEqual(aln_seq, msa.decode_alignment_as_list_of_pairs())

    def test_should_merge_gaps_groups(self):
        # todo
        pass

    def test_should_return_is_gap_column(self):
        # setup
        msa = MSASolution(aligned_sequences=[('seq1', 'AC---TGAC'), ('seq2', 'AT--CT--C'), ('seq3', 'AAC---TGC')],
                          number_of_objectives=2)

        # check
        self.assertTrue(msa.is_gap_column(3))
        self.assertFalse(msa.is_gap_column(4))

    def test_should_split_gap_column(self):
        # setup
        msa_1 = MSASolution(aligned_sequences=[('seq1', '---AC'), ('seq2', 'T---C'), ('seq3', '--A-A-')],
                            number_of_objectives=2)
        msa_2 = MSASolution(aligned_sequences=[('seq1', '---AC'), ('seq2', 'T---C'), ('seq3', '--A-A-')],
                            number_of_objectives=2)
        msa_3 = MSASolution(aligned_sequences=[('seq1', '---AC'), ('seq2', 'T---C'), ('seq3', '-A----')],
                            number_of_objectives=2)

        # actual = [[0, 2], [1, 3], [0, 1, 3, 3, 5, 5]]
        msa_1.split_gap_column(1)
        msa_2.split_gap_column(2)

        # actual = [[0, 2], [1, 3], [0, 0, 2, 5]]
        msa_3.split_gap_column(4)

        # check
        self.assertEqual([[0, 1, 1, 2], [1, 3], [0, 1, 3, 3, 5, 5]], msa_1.gaps_groups)
        self.assertEqual([[0, 2], [1, 2, 2, 3], [0, 1, 3, 3, 5, 5]], msa_2.gaps_groups)
        self.assertEqual([[0, 2], [1, 3], [0, 0, 2, 4, 4, 5]], msa_3.gaps_groups)

    def test_should_remove_gap_column(self):
        # setup
        msa = MSASolution(aligned_sequences=[('seq1', 'AC---TGAC'), ('seq2', 'AT--CT--C'), ('seq3', 'AAC---TGC')],
                          number_of_objectives=2)
        msa.remove_gap_column(3)

        # check
        self.assertEqual([[2, 2, 4, 4], [2, 2, 6, 7], [4, 5]], msa.gaps_groups)

    def test_should_return_gap_columns(self):
        # setup
        msa = MSASolution(aligned_sequences=[('seq1', '--AA-'), ('seq2', '--AA-'), ('seq3', '--AA-')],
                          number_of_objectives=2)

        # check
        self.assertEqual([0, 1, 4], msa.get_gap_columns())


if __name__ == "__main__":
    unittest.main()
