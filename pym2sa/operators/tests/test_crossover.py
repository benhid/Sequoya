import unittest

from pym2sa.core.solution import MSASolution
from pym2sa.operators.crossover import GapSequenceSolutionSinglePoint

'''
class GapSequenceSolutionSinglePointTestCases(unittest.TestCase):
    def setUp(self):
        print("setUp: RUNNING GapSequenceSolutionSinglePointTestCases")

    def tearDown(self):
        print("tearDown: TEST ENDED")

    def test_should_the_solution_remain_unchanged_if_the_probability_is_zero(self):
        # setup
        msa_1 = MSASolution(aligned_sequences=[('seq1', 'ACTC'), ('seq2', 'A-TC'), ('seq3', 'A--C')],
                            number_of_objectives=3)
        msa_2 = MSASolution(aligned_sequences=[('seq1', 'CT-G'), ('seq2', '-T-G'), ('seq3', '-ATG')],
                            number_of_objectives=3)
        crossover = GapSequenceSolutionSinglePoint(probability=0.0)

        # run
        offspring = crossover.execute([msa_1, msa_2])

        # check
        self.assertEqual([('seq1', 'ACTC'), ('seq2', 'A-TC'), ('seq3', 'A--C')], offspring[0].decode_alignment_as_list_of_pairs())
        self.assertEqual([('seq1', 'CT-G'), ('seq2', '-T-G'), ('seq3', '-ATG')], offspring[1].decode_alignment_as_list_of_pairs())

    def test_should_split_two_parents_with_no_gaps_in_half(self):
        # setup
        msa_1 = MSASolution(aligned_sequences=[('seq1', 'AB'), ('seq2', 'CD')], number_of_objectives=2)
        msa_2 = MSASolution(aligned_sequences=[('seq1', 'EF'), ('seq2', 'GH')], number_of_objectives=2)

        # run
        crossover = GapSequenceSolutionSinglePoint(probability=1.0)
        a1, a2, b1, b2 = crossover.split_parents(cx_point=0, parent_a=msa_1, parent_b=msa_2)

        # check
        self.assertEqual([('seq1', 'A'), ('seq2', 'C')], a1.decode_alignment_as_list_of_pairs())
        self.assertEqual([('seq1', 'B'), ('seq2', 'D')], a2.decode_alignment_as_list_of_pairs())
        self.assertEqual([('seq1', 'E'), ('seq2', 'G')], b1.decode_alignment_as_list_of_pairs())
        self.assertEqual([('seq1', 'F'), ('seq2', 'H')], b2.decode_alignment_as_list_of_pairs())

    def test_should_split_two_parents_with_one_gap_in_half(self):
        # setup
        msa_1 = MSASolution(aligned_sequences=[('seq1', 'A-JB'), ('seq2', '-C-D')], number_of_objectives=2)
        msa_2 = MSASolution(aligned_sequences=[('seq1', 'E-F-'), ('seq2', '-G-H')], number_of_objectives=2)

        # run
        crossover = GapSequenceSolutionSinglePoint(probability=1.0)
        a1, a2, b1, b2 = crossover.split_parents(cx_point=1, parent_a=msa_1, parent_b=msa_2)

        # check
        self.assertEqual([('seq1', 'A-'), ('seq2', '-C')], a1.decode_alignment_as_list_of_pairs())
        self.assertEqual([('seq1', 'JB'), ('seq2', '-D')], a2.decode_alignment_as_list_of_pairs())
        self.assertEqual([('seq1', 'E-'), ('seq2', '-G')], b1.decode_alignment_as_list_of_pairs())
        self.assertEqual([('seq1', 'F-'), ('seq2', '-H')], b2.decode_alignment_as_list_of_pairs())

    def test_should_split_two_parents_with_multiple_gaps_in_half(self):
        # setup
        msa_1 = MSASolution(aligned_sequences=[('seq1', 'AB--CD-E')], number_of_objectives=2)
        msa_2 = MSASolution(aligned_sequences=[('seq1', 'AB-CDE-')], number_of_objectives=2)

        # run
        crossover = GapSequenceSolutionSinglePoint(probability=1.0)
        a1, a2, b1, b2 = crossover.split_parents(cx_point=4, parent_a=msa_1, parent_b=msa_2)

        # check
        self.assertEqual([('seq1', 'AB--C')], a1.decode_alignment_as_list_of_pairs())
        self.assertEqual([('seq1', 'D-E')], a2.decode_alignment_as_list_of_pairs())
        self.assertEqual([('seq1', 'AB-C')], b1.decode_alignment_as_list_of_pairs())
        self.assertEqual([('seq1', 'DE-')], b2.decode_alignment_as_list_of_pairs())

    def test_should_merge_two_parents_with_no_gaps(self):
        # setup
        msa_1 = MSASolution(aligned_sequences=[('seq1', 'BC'), ('seq2', 'AB')], number_of_objectives=2)
        msa_2 = MSASolution(aligned_sequences=[('seq1', 'BC'), ('seq2', 'AB')], number_of_objectives=2)

        # run
        crossover = GapSequenceSolutionSinglePoint(probability=1.0)
        individual = crossover.merge_parents(msa_1, msa_2)

        # check
        self.assertEqual([('seq1', 'BCBC'), ('seq2', 'ABAB')], individual.decode_alignment_as_list_of_pairs())

    def test_should_merge_two_parents_with_one_gaps(self):
        # setup
        msa_1 = MSASolution(aligned_sequences=[('seq1', 'B-C'), ('seq2', 'AB-')], number_of_objectives=2)
        msa_2 = MSASolution(aligned_sequences=[('seq1', 'BC-'), ('seq2', 'A-B')], number_of_objectives=2)

        # run
        crossover = GapSequenceSolutionSinglePoint(probability=1.0)
        individual = crossover.merge_parents(msa_1, msa_2)

        # check
        self.assertEqual([('seq1', 'B-CBC-'), ('seq2', 'AB-A-B')], individual.decode_alignment_as_list_of_pairs())

    def test_should_execute_if_parents_have_no_gaps(self):
        # setup
        msa_1 = MSASolution(aligned_sequences=[('seq1', 'BCDE'), ('seq2', 'ABCE')], number_of_objectives=2)
        msa_2 = MSASolution(aligned_sequences=[('seq1', 'BCDE'), ('seq2', 'ABCE')], number_of_objectives=2)

        # run
        crossover = GapSequenceSolutionSinglePoint(probability=1.0)
        offspring = crossover.execute([msa_1, msa_2])

        # check
        self.assertEqual([('seq1', 'BCDE'), ('seq2', 'ABCE')], offspring[0].decode_alignment_as_list_of_pairs())
        self.assertEqual([('seq1', 'BCDE'), ('seq2', 'ABCE')], offspring[1].decode_alignment_as_list_of_pairs())
'''


class GapSequenceSolutionSinglePointTestCases(unittest.TestCase):
    def setUp(self):
        print("setUp: RUNNING GapSequenceSolutionSinglePointTestCases")

    def tearDown(self):
        print("tearDown: TEST ENDED")

    def test_should_the_solution_remain_unchanged_if_the_probability_is_zero(self):
        # setup
        msa_1 = MSASolution(aligned_sequences=[('seq1', 'ACTC'), ('seq2', 'A-TC'), ('seq3', 'A--C')],
                            number_of_objectives=3)
        msa_2 = MSASolution(aligned_sequences=[('seq1', 'CT-G'), ('seq2', '-T-G'), ('seq3', '-ATG')],
                            number_of_objectives=3)
        crossover = GapSequenceSolutionSinglePoint(probability=0.0)

        # run
        offspring = crossover.execute([msa_1, msa_2])

        # check
        self.assertEqual([('seq1', 'ACTC'), ('seq2', 'A-TC'), ('seq3', 'A--C')], offspring[0].decode_alignment_as_list_of_pairs())
        self.assertEqual([('seq1', 'CT-G'), ('seq2', '-T-G'), ('seq3', '-ATG')], offspring[1].decode_alignment_as_list_of_pairs())



if __name__ == "__main__":
    unittest.main()
