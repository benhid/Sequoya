import unittest
from unittest import mock

from pym2sa.core.solution import MSASolution
from pym2sa.operators.crossover import SinglePointMSA, SinglePointXMSA

"""
class SinglePointMSATestCases(unittest.TestCase):
    def setUp(self):
        print("setUp: RUNNING SinglePointMSATestCases")

    def tearDown(self):
        print("tearDown: TEST ENDED")

    def test_should_the_solution_remain_unchanged_if_the_probability_is_almost_zero(self):
        # setup
        msa_1 = MSASolution(aligned_sequences=[('seq1', 'ACTC'), ('seq2', 'A-TC'), ('seq3', 'A--C')],
                            number_of_objectives=2)
        msa_2 = MSASolution(aligned_sequences=[('seq1', 'CT-G'), ('seq2', '-T-G'), ('seq3', '-ATG')],
                            number_of_objectives=2)
        crossover = SinglePointMSA(probability=0.001)

        # run
        offspring = crossover.execute([msa_1, msa_2])

        # check
        self.assertEqual([('seq1', 'ACTC'), ('seq2', 'A-TC'), ('seq3', 'A--C')], offspring[0].decode_alignment_as_list_of_pairs())
        self.assertEqual([('seq1', 'CT-G'), ('seq2', '-T-G'), ('seq3', '-ATG')], offspring[1].decode_alignment_as_list_of_pairs())

    def test_should_split_first_half_of_parent_until_third_character(self):
        # setup
        msa_1 = MSASolution(aligned_sequences=[('seq1', 'ACT---C'), ('seq2', '---A-TC'), ('seq3', 'A--C---')],
                            number_of_objectives=2)
        msa_2 = MSASolution(aligned_sequences=[('seq1', 'C-T---G'), ('seq2', '-T---G-'), ('seq3', '----ATG')],
                            number_of_objectives=2)

        # actual msa_1 = [[3, 5], [0, 2, 4, 4], [1, 2, 4, 6]]
        # actual msa_2 = [[1, 1, 3, 5], [0, 0, 2, 4, 6, 6], [0, 3]]

        # run
        crossover = SinglePointMSA(probability=0.999)
        new_msa_1_first_half, new_msa_1_second_half = crossover.split(cx_point=2, parent=msa_1)
        new_msa_2_first_half, new_msa_2_second_half = crossover.split(cx_point=2, parent=msa_2)

        # check
        self.assertEqual([[], [0, 2], [1, 2]], new_msa_1_first_half.gaps_groups)
        self.assertEqual([[3, 5], [4, 4], [4, 6]], new_msa_1_second_half.gaps_groups)
        self.assertEqual(['ACT', '', 'A'], new_msa_1_first_half.variables)
        self.assertEqual(['C', 'ATC', 'C'], new_msa_1_second_half.variables)

        self.assertEqual([[1, 1], [0, 0, 2, 2], [0, 2]], new_msa_2_first_half.gaps_groups)
        self.assertEqual([[3, 5], [3, 4, 6, 6], [3, 3]], new_msa_2_second_half.gaps_groups)
        self.assertEqual(['CT', 'T', ''], new_msa_2_first_half.variables)
        self.assertEqual(['G', 'G', 'ATG'], new_msa_2_second_half.variables)

    def test_should_split_first_half_of_parent_until_fourth_character(self):
        # setup
        msa_1 = MSASolution(aligned_sequences=[('seq1', 'A-CT--C'), ('seq2', '-A-A-TC'), ('seq3', 'A--C---')],
                            number_of_objectives=2)
        msa_2 = MSASolution(aligned_sequences=[('seq1', 'A--CT-C'), ('seq2', '-AA--TC'), ('seq3', 'A--C---')],
                            number_of_objectives=2)

        # run
        crossover = SinglePointMSA(probability=0.999)
        new_msa_1_first_half, new_msa_1_second_half = crossover.split(cx_point=2, parent=msa_1)
        new_msa_2_first_half, new_msa_2_second_half = crossover.split(cx_point=2, parent=msa_2)

        # check
        self.assertEqual([[], [0, 2], [1, 2]], new_msa_1_first_half.gaps_groups)
        self.assertEqual([[3, 5], [4, 4], [4, 6]], new_msa_1_second_half.gaps_groups)

        self.assertEqual([[1, 1], [0, 0, 2, 2], [0, 2]], new_msa_2_first_half.gaps_groups)
        self.assertEqual([[3, 5], [3, 4, 6, 6], [3, 3]], new_msa_2_second_half.gaps_groups)
        self.assertEqual(['CT', 'T', ''], new_msa_2_first_half.variables)
        self.assertEqual(['G', 'G', 'ATG'], new_msa_2_second_half.variables)

    def test_should_merge_two_parents(self):
        # setup
        msa_1 = MSASolution(aligned_sequences=[('seq1', 'ACT'), ('seq2', '---'), ('seq3', 'A--')],
                            number_of_objectives=2)
        msa_2 = MSASolution(aligned_sequences=[('seq1', '---C'), ('seq2', 'A-TC'), ('seq3', 'C---')],
                            number_of_objectives=2)

        # run
        crossover = SinglePointMSA(probability=0.999)
        new_msa = crossover.merge(parents=[msa_1, msa_2])

        # check
        self.assertEqual([[3, 5], [0, 2, 4, 4], [1, 2, 4, 6]], new_msa[0].gaps_groups)
        self.assertEqual(['ACTC', 'ATC', 'AC'], new_msa[0].variables)
        self.assertEqual([[0, 2], [1, 1, 4, 6], [1, 3, 5, 6]], new_msa[1].gaps_groups)
        self.assertEqual(['CACT', 'ATC', 'CA'], new_msa[1].variables)

    @mock.patch('random.randint')
    def test_should_the_operator_work_if_the_third_bit_is_selected(self, random_call):
        # setup
        msa_1 = MSASolution(aligned_sequences=[('seq1', 'A-CT--C'), ('seq2', '-A-A-TC'), ('seq3', 'A--C---')],
                            number_of_objectives=2)
        msa_2 = MSASolution(aligned_sequences=[('seq1', 'A--CT-C'), ('seq2', '-AA--TC'), ('seq3', 'A--C---')],
                            number_of_objectives=2)

        # run
        crossover = SinglePointMSA(probability=0.999)
        random_call.return_value = 3

        offspring = crossover.execute([msa_1, msa_2])

        # check
        self.assertEqual([('seq1', 'A-CTTC-'), ('seq2', '-A-ATC-'), ('seq3', 'A--C---')],
                         offspring[0].decode_alignment_as_list_of_pairs())
        self.assertEqual([('seq1', 'TCAC-T-'), ('seq2', 'TCA-A--'), ('seq3', 'AC-----')],
                         offspring[1].decode_alignment_as_list_of_pairs())
"""

class SinglePointXMSATestCases(unittest.TestCase):
    def setUp(self):
        print("setUp: RUNNING SinglePointXMSATestCases")

    def tearDown(self):
        print("tearDown: TEST ENDED")

    def test_should_the_solution_remain_unchanged_if_the_probability_is_zero(self):
        # setup
        msa_1 = MSASolution(aligned_sequences=[('seq1', 'ACTC'), ('seq2', 'A-TC'), ('seq3', 'A--C')],
                            number_of_objectives=3)
        msa_2 = MSASolution(aligned_sequences=[('seq1', 'CT-G'), ('seq2', '-T-G'), ('seq3', '-ATG')],
                            number_of_objectives=3)
        crossover = SinglePointXMSA(probability=0.0)

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
        crossover = SinglePointXMSA(probability=1.0)
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
        crossover = SinglePointXMSA(probability=1.0)
        a1, a2, b1, b2 = crossover.split_parents(cx_point=1, parent_a=msa_1, parent_b=msa_2)

        # check
        self.assertEqual([('seq1', 'A-'), ('seq2', '-C')], a1.decode_alignment_as_list_of_pairs())
        self.assertEqual([('seq1', 'JB'), ('seq2', '-D')], a2.decode_alignment_as_list_of_pairs())
        self.assertEqual([('seq1', 'E-'), ('seq2', '-G')], b1.decode_alignment_as_list_of_pairs())
        self.assertEqual([('seq1', 'F-'), ('seq2', '-H')], b2.decode_alignment_as_list_of_pairs())

    def test_should_merge_two_parents_with_no_gaps(self):
        # setup
        msa_1 = MSASolution(aligned_sequences=[('seq1', 'BC'), ('seq2', 'AB')], number_of_objectives=2)
        msa_2 = MSASolution(aligned_sequences=[('seq1', 'BC'), ('seq2', 'AB')], number_of_objectives=2)

        # run
        crossover = SinglePointXMSA(probability=1.0)
        individual = crossover.merge_parents(msa_1, msa_2)

        # check
        self.assertEqual([('seq1', 'BCBC'), ('seq2', 'ABAB')], individual.decode_alignment_as_list_of_pairs())

    def test_should_merge_two_parents_with_one_gaps(self):
        # setup
        msa_1 = MSASolution(aligned_sequences=[('seq1', 'B-C'), ('seq2', 'AB-')], number_of_objectives=2)
        msa_2 = MSASolution(aligned_sequences=[('seq1', 'BC-'), ('seq2', 'A-B')], number_of_objectives=2)

        # run
        crossover = SinglePointXMSA(probability=1.0)
        individual = crossover.merge_parents(msa_1, msa_2)

        # check
        self.assertEqual([('seq1', 'B-CBC-'), ('seq2', 'AB-A-B')], individual.decode_alignment_as_list_of_pairs())

    def test_should_execute_if_parents_have_no_gaps(self):
        # setup
        msa_1 = MSASolution(aligned_sequences=[('seq1', 'BCDE'), ('seq2', 'ABCE')], number_of_objectives=2)
        msa_2 = MSASolution(aligned_sequences=[('seq1', 'BCDE'), ('seq2', 'ABCE')], number_of_objectives=2)

        # run
        crossover = SinglePointXMSA(probability=1.0)
        offspring = crossover.execute([msa_1, msa_2])

        # check
        self.assertEqual([('seq1', 'BCDE'), ('seq2', 'ABCE')], offspring[0].decode_alignment_as_list_of_pairs())
        self.assertEqual([('seq1', 'BCDE'), ('seq2', 'ABCE')], offspring[1].decode_alignment_as_list_of_pairs())


if __name__ == "__main__":
    unittest.main()
