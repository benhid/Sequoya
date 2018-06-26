import unittest
from unittest import mock

from pym2sa.core.solution import MSASolution
from pym2sa.operator.mutation import OneRandomGapInsertion, TwoRandomAdjacentGapGroup, ShiftGapGroup


class ShiftGapGroupTestCases(unittest.TestCase):

    def setUp(self):
        print("setUp: RUNNING ShiftGapGroupTestCases")

    def tearDown(self):
        print("tearDown: TEST ENDED")

    @mock.patch('random.randrange')
    @mock.patch('random.randint')
    def test_should_execute_mutation_case_a(self, random_shift, random_group):
        # setup
        msa = MSASolution(aligned_sequences=[('seq1', '--AB--CD--')], number_of_objectives=2)
        mutation = ShiftGapGroup(probability=1.0, remove_gap_columns=False)

        # run
        random_group.return_value = 2
        random_shift.return_value = 0
        solution = mutation.execute(msa)

        # check
        self.assertEqual([[0, 1, 3, 4, 8, 9]], solution.gaps_groups)
        self.assertEqual([('seq1', '--A--BCD--')], solution.decode_alignment_as_list_of_pairs())

    @mock.patch('random.randrange')
    @mock.patch('random.randint')
    def test_should_execute_mutation_case_a(self, random_shift, random_group):
        # setup
        msa = MSASolution(aligned_sequences=[('seq1', '--AB--CD--')], number_of_objectives=2)
        mutation = ShiftGapGroup(probability=1.0, remove_gap_columns=False)

        # run
        random_group.return_value = 2
        random_shift.return_value = 1
        solution = mutation.execute(msa)

        # check
        self.assertEqual([[0, 1, 5, 6, 8, 9]], solution.gaps_groups)
        self.assertEqual([('seq1', '--ABC--D--')], solution.decode_alignment_as_list_of_pairs())

    @mock.patch('random.randrange')
    @mock.patch('random.randint')
    def test_should_execute_mutation_case_a(self, random_shift, random_group):
        # setup
        msa = MSASolution(aligned_sequences=[('seq1', '--AB--CD--')], number_of_objectives=2)
        mutation = ShiftGapGroup(probability=1.0, remove_gap_columns=False)

        # run
        random_group.return_value = 0
        random_shift.return_value = 1
        solution = mutation.execute(msa)

        # check
        self.assertEqual([('seq1', 'A--B--CD--')], solution.decode_alignment_as_list_of_pairs())

    @mock.patch('random.randrange')
    @mock.patch('random.randint')
    def test_should_execute_mutation_case_a(self, random_shift, random_group):
        # setup
        msa = MSASolution(aligned_sequences=[('seq1', '--AB--C--D')], number_of_objectives=2)
        mutation = ShiftGapGroup(probability=1.0, remove_gap_columns=False)

        # run
        random_group.return_value = 2
        random_shift.return_value = 1
        solution = mutation.execute(msa)

        # check
        self.assertEqual([('seq1', '--ABC----D')], solution.decode_alignment_as_list_of_pairs())


class TwoRandomAdjacentGapGroupTestCases(unittest.TestCase):

    def setUp(self):
        print("setUp: RUNNING TwoRandomAdjacentGapGroupTestCases")

    def tearDown(self):
        print("tearDown: TEST ENDED")

    def test_should_execute_mutation_case_a(self):
        # setup
        msa = MSASolution(aligned_sequences=[('seq1', 'AAA--B--D')], number_of_objectives=2)
        mutation = TwoRandomAdjacentGapGroup(probability=1.0, remove_gap_columns=False)

        # run
        solution = mutation.execute(msa)

        # check
        self.assertEqual([('seq1', 'AAA----BD')], solution.decode_alignment_as_list_of_pairs())

    def test_should_execute_mutation_case_b(self):
        # setup
        msa = MSASolution(aligned_sequences=[('seq1', '--CD--')], number_of_objectives=2)
        mutation = TwoRandomAdjacentGapGroup(probability=1.0, remove_gap_columns=False)

        # run
        solution = mutation.execute(msa)

        # check
        self.assertEqual([('seq1', '----CD')], solution.decode_alignment_as_list_of_pairs())

    @mock.patch('random.randrange')
    def test_should_execute_mutation_case_c(self, random_call):
        # setup
        msa = MSASolution(aligned_sequences=[('seq1', '-CD--D-')], number_of_objectives=2)
        mutation = TwoRandomAdjacentGapGroup(probability=1.0, remove_gap_columns=False)

        # run
        random_call.return_value = 0
        solution = mutation.execute(msa)

        # check
        self.assertEqual([('seq1', '---CDD-')], solution.decode_alignment_as_list_of_pairs())

    @mock.patch('random.randrange')
    def test_should_execute_mutation_case_d(self, random_call):
        # setup
        msa = MSASolution(aligned_sequences=[('seq1', '-A-B-C-')], number_of_objectives=2)
        mutation = TwoRandomAdjacentGapGroup(probability=1.0, remove_gap_columns=False)

        # run
        random_call.return_value = 2
        solution = mutation.execute(msa)

        # check
        self.assertEqual([('seq1', '-A--BC-')], solution.decode_alignment_as_list_of_pairs())


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
