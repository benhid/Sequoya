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
    def test_should_execute_mutation(self, random_call):
        # setup
        msa = MSASolution(aligned_sequences=[('seq1', 'AB-D')], number_of_objectives=2)
        mutation = OneRandomGapInsertion(probability=1.0)

        # run
        random_call.return_value = 1
        solution = mutation.execute(msa)

        # check
        self.assertEqual([('seq1', 'A-B-D')], solution.decode_alignment_as_list_of_pairs())


if __name__ == "__main__":
    unittest.main()
