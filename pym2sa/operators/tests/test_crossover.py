import unittest
from unittest import mock

from pym2sa.core.solution import MSASolution
from pym2sa.operators.crossover import SinglePointMSA


class SinglePointMSATestCases(unittest.TestCase):
    def setUp(self):
        print("setUp: RUNNING SinglePointMSATestCases")

    def tearDown(self):
        print("tearDown: TEST ENDED")

    def test_should_the_solution_remain_unchanged_if_the_probability_is_zero(self):
        # setup
        msa_1 = MSASolution(aligned_sequences=[('seq1', 'ACTC'), ('seq2', 'A-TC'), ('seq3', 'A--C')],
                            number_of_objectives=2)
        msa_2 = MSASolution(aligned_sequences=[('seq1', 'CT-G'), ('seq2', '-T-G'), ('seq3', '-ATG')],
                            number_of_objectives=2)
        crossover = SinglePointMSA(probability=0)

        # run
        offspring = crossover.execute([msa_1, msa_2])

        # check
        self.assertEqual([('seq1', 'ACTC'), ('seq2', 'A-TC'), ('seq3', 'A--C')], offspring[0].decode_alignment_as_list_of_pairs())
        self.assertEqual([('seq1', 'CT-G'), ('seq2', '-T-G'), ('seq3', '-ATG')], offspring[1].decode_alignment_as_list_of_pairs())

    def test_should_split_first_half_of_parent(self):
        # setup
        msa_1 = MSASolution(aligned_sequences=[('seq1', 'ACT---C'), ('seq2', '---A-TC'), ('seq3', 'A--C---')],
                            number_of_objectives=2)
        msa_2 = MSASolution(aligned_sequences=[('seq1', 'C-T---G'), ('seq2', '-T---G-'), ('seq3', '----ATG')],
                            number_of_objectives=2)

        # actual msa_1 = [[3, 5], [0, 2, 4, 4], [1, 2, 4, 6]]
        # actual msa_2 = [[1, 1, 3, 5], [0, 0, 2, 4, 6, 6], [0, 3]]

        # run
        crossover = SinglePointMSA(probability=0.999)
        new_msa_1 = crossover.split(cx_point=2, parent=msa_1)
        new_msa_2 = crossover.split(cx_point=2, parent=msa_2)

        # check
        self.assertEqual([[], [0, 2], [1, 2]], new_msa_1.gaps_groups)
        self.assertEqual(['ACT', '', 'A'], new_msa_1.variables)

        self.assertEqual([[1, 1], [0, 0], [0, 2]], new_msa_2.gaps_groups)
        self.assertEqual(['CT', 'T', ''], new_msa_2.variables)

    @mock.patch('random.randint')
    def test_should_the_operator_work_if_the_second_bit_is_selected(self, random_call):
        # setup
        msa_1 = MSASolution(aligned_sequences=[('seq1', 'ACT---C'), ('seq2', '---A-TC'), ('seq3', 'A--C---')],
                            number_of_objectives=2)
        msa_2 = MSASolution(aligned_sequences=[('seq1', 'C-T---G'), ('seq2', '-T---G-'), ('seq3', '----ATG')],
                            number_of_objectives=2)

        # actual msa_1 = [[3, 5], [0, 2, 4, 4], [1, 2, 4, 6]]
        # actual msa_2 = [[1, 1, 3, 5], [0, 0, 2, 4, 6, 6], [0, 3]]

        # run
        crossover = SinglePointMSA(probability=0.999)
        random_call.return_value = 2  # first 3 characters (0, 1, 2)

        offspring = crossover.execute([msa_1, msa_2])

        # check
        self.assertEqual([[], [0, 2], [1, 2]], offspring[0].gaps_groups)
        self.assertEqual(['ACT', '', 'A'], offspring[0].variables)

        self.assertEqual([[1, 1], [0, 0], []], offspring[1].gaps_groups)
        self.assertEqual(['ACT', '', 'A'], offspring[1].variables)


if __name__ == "__main__":
    unittest.main()
