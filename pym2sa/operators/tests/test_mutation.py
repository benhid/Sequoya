import unittest
from unittest import mock

from pym2sa.core.solution import MSASolution
from pym2sa.operators.mutation import OneRandomGapInsertion


class OneRandomGapInsertionTestCases(unittest.TestCase):

    def setUp(self):
        print("setUp: RUNNING OneRandomGapInsertionTestCases")

    def tearDown(self):
        print("tearDown: TEST ENDED")

    @mock.patch('random.randint')
    def test_should_execute_mutation_case_a(self, random_call):
        # setup
        msa = MSASolution(aligned_sequences=[('seq1', 'AB-D')], number_of_objectives=2)
        mutation = OneRandomGapInsertion(probability=1.0)

        # run
        random_call.return_value = 1
        solution = mutation.execute(msa)

        # check
        self.assertEqual([('seq1', 'A-B-D')], solution.decode_alignment_as_list_of_pairs())

    @mock.patch('random.randint')
    def test_should_execute_mutation_case_b(self, random_call):
        # setup
        msa = MSASolution(aligned_sequences=[('seq1', 'AB---D'), ('seq1', 'A--B-D')], number_of_objectives=2)
        mutation = OneRandomGapInsertion(probability=1.0)

        # run
        random_call.return_value = 1
        solution = mutation.execute(msa)

        # check
        self.assertEqual([('seq1', 'A-B---D'), ('seq1', 'A---B-D')], solution.decode_alignment_as_list_of_pairs())

    @mock.patch('random.randint')
    def test_should_execute_mutation_case_c(self, random_call):
        # setup
        msa = MSASolution(aligned_sequences=[('seq1', 'ABBB'), ('seq1', '----')], number_of_objectives=2)
        mutation = OneRandomGapInsertion(probability=1.0)

        # run
        random_call.return_value = 1
        solution = mutation.execute(msa)

        # check
        self.assertEqual([('seq1', 'A-BBB'), ('seq1', '-----')], solution.decode_alignment_as_list_of_pairs())

    @mock.patch('random.randint')
    def test_should_execute_mutation_case_d(self, random_call):
        # setup
        msa = MSASolution(aligned_sequences=[('seq1', 'A-G-A-GA-G-G-T-'), ('seq1', '--CAC-A--GGT--G')], number_of_objectives=2)
        mutation = OneRandomGapInsertion(probability=1.0)

        # run
        random_call.return_value = 8
        solution = mutation.execute(msa)

        # check
        self.assertEqual([('seq1', 'A-G-A-GA--G-G-T-'), ('seq1', '--CAC-A---GGT--G')], solution.decode_alignment_as_list_of_pairs())


if __name__ == "__main__":
    unittest.main()
