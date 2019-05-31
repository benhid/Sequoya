import unittest
from unittest import mock

from sequoya.core.solution import MSASolution
from sequoya.operator import OneRandomGapInsertion, TwoRandomAdjacentGapGroup, ShiftGapGroup, ShiftClosedGapGroups
from sequoya.problem import MSA


class ShiftClosedGapGroupsTestCases(unittest.TestCase):

    def setUp(self):
        self.problem = MSA(score_list=[])
        self.problem.identifiers = ['seq1']
        self.problem.number_of_variables = 1

    @mock.patch('random.randrange')
    def test_should_execute_mutation_case_a(self, random_group):
        # setup
        msa = MSASolution(self.problem, msa=[('seq1', '--AB--CD--')])
        mutation = ShiftClosedGapGroups(probability=1.0, remove_gap_columns=False)

        # run
        random_group.return_value = 2
        solution = mutation.execute(msa)

        # check
        self.assertEqual([('seq1', '--AB--CD--')], solution.decode_alignment_as_list_of_pairs())

    @mock.patch('random.randrange')
    def test_should_execute_mutation_case_b(self, random_group):
        # setup
        msa = MSASolution(self.problem, msa=[('seq1', '--AB--CD---E')])
        mutation = ShiftClosedGapGroups(probability=1.0, remove_gap_columns=False)

        # run
        random_group.return_value = 2
        solution = mutation.execute(msa)

        # check
        self.assertEqual([('seq1', '--AB---CD--E')], solution.decode_alignment_as_list_of_pairs())

    @mock.patch('random.randrange')
    def test_should_execute_mutation_case_c(self, random_group):
        # setup
        msa = MSASolution(self.problem, msa=[('seq1', '--AB-----CD---E')])
        mutation = ShiftClosedGapGroups(probability=1.0, remove_gap_columns=False)

        # run
        random_group.return_value = 2
        solution = mutation.execute(msa)

        # check
        self.assertEqual([('seq1', '--AB---CD-----E')], solution.decode_alignment_as_list_of_pairs())


    @mock.patch('random.randrange')
    def test_should_execute_mutation_case_d(self, random_group):
        # setup
        msa = MSASolution(self.problem, msa=[('seq1', '--A--')])
        mutation = ShiftClosedGapGroups(probability=1.0, remove_gap_columns=False)

        # run
        random_group.return_value = 0
        solution = mutation.execute(msa)

        # check
        self.assertEqual([('seq1', '--A--')], solution.decode_alignment_as_list_of_pairs())

    @mock.patch('random.randrange')
    def test_should_execute_mutation_case_e(self, random_group):
        # setup
        problem = MSA(score_list=[])
        problem.identifiers = ['seq1', 'seq2']
        problem.number_of_variables = 2

        msa = MSASolution(problem, msa=[('seq1', '---B--AA--'), ('seq2', '--B---AA--')])
        mutation = ShiftClosedGapGroups(probability=1.0, remove_gap_columns=False)

        # run
        random_group.return_value = 0
        solution = mutation.execute(msa)

        # check
        self.assertEqual([('seq1', '--B---AA--'), ('seq2', '---B--AA--')], solution.decode_alignment_as_list_of_pairs())


class ShiftGapGroupTestCases(unittest.TestCase):

    def setUp(self):
        self.problem = MSA(score_list=[])
        self.problem.identifiers = ['seq1']
        self.problem.number_of_variables = 1

    @mock.patch('random.randrange')
    @mock.patch('random.randint')
    def test_should_execute_mutation_case_a(self, random_shift, random_group):
        # setup
        msa = MSASolution(self.problem, msa=[('seq1', '--AB--CD--')])
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
    def test_should_execute_mutation_case_b(self, random_shift, random_group):
        # setup
        msa = MSASolution(self.problem, msa=[('seq1', '--AB--CD--')])
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
    def test_should_execute_mutation_case_c(self, random_shift, random_group):
        # setup
        msa = MSASolution(self.problem, msa=[('seq1', '--AB--CD--')])
        mutation = ShiftGapGroup(probability=1.0, remove_gap_columns=False)

        # run
        random_group.return_value = 0
        random_shift.return_value = 1
        solution = mutation.execute(msa)

        # check
        self.assertEqual([('seq1', 'A--B--CD--')], solution.decode_alignment_as_list_of_pairs())

    @mock.patch('random.randrange')
    @mock.patch('random.randint')
    def test_should_execute_mutation_case_d(self, random_shift, random_group):
        # setup
        msa = MSASolution(self.problem, msa=[('seq1', '--AB--C--D')])
        mutation = ShiftGapGroup(probability=1.0, remove_gap_columns=False)

        # run
        random_group.return_value = 2
        random_shift.return_value = 1
        solution = mutation.execute(msa)

        # check
        self.assertEqual([('seq1', '--ABC----D')], solution.decode_alignment_as_list_of_pairs())


class TwoRandomAdjacentGapGroupTestCases(unittest.TestCase):

    def setUp(self):
        self.problem = MSA(score_list=[])
        self.problem.identifiers = ['seq1']
        self.problem.number_of_variables = 1

    def test_should_execute_mutation_case_a(self):
        # setup
        msa = MSASolution(self.problem, msa=[('seq1', 'AAA--B--D')])
        mutation = TwoRandomAdjacentGapGroup(probability=1.0, remove_gap_columns=False)

        # run
        solution = mutation.execute(msa)

        # check
        self.assertEqual([('seq1', 'AAA----BD')], solution.decode_alignment_as_list_of_pairs())

    def test_should_execute_mutation_case_b(self):
        # setup
        msa = MSASolution(self.problem, msa=[('seq1', '--CD--')])
        mutation = TwoRandomAdjacentGapGroup(probability=1.0, remove_gap_columns=False)

        # run
        solution = mutation.execute(msa)

        # check
        self.assertEqual([('seq1', '----CD')], solution.decode_alignment_as_list_of_pairs())

    @mock.patch('random.randrange')
    def test_should_execute_mutation_case_c(self, random_call):
        # setup
        msa = MSASolution(self.problem, msa=[('seq1', '-CD--D-')])
        mutation = TwoRandomAdjacentGapGroup(probability=1.0, remove_gap_columns=False)

        # run
        random_call.return_value = 0
        solution = mutation.execute(msa)

        # check
        self.assertEqual([('seq1', '---CDD-')], solution.decode_alignment_as_list_of_pairs())

    @mock.patch('random.randrange')
    def test_should_execute_mutation_case_d(self, random_call):
        # setup
        msa = MSASolution(self.problem, msa=[('seq1', '-A-B-C-')])
        mutation = TwoRandomAdjacentGapGroup(probability=1.0, remove_gap_columns=False)

        # run
        random_call.return_value = 2
        solution = mutation.execute(msa)

        # check
        self.assertEqual([('seq1', '-A--BC-')], solution.decode_alignment_as_list_of_pairs())


class OneRandomGapInsertionTestCases(unittest.TestCase):

    def setUp(self):
        self.problem = MSA(score_list=[])
        self.problem.identifiers = ['seq1', 'seq2']
        self.problem.number_of_variables = 2

    @mock.patch('random.randint')
    def test_should_execute_mutation_case_a(self, random_call):
        # setup
        problem = MSA(score_list=[])
        problem.identifiers = ['seq1']
        problem.number_of_variables = 1
        msa = MSASolution(problem, msa=[('seq1', 'AB-D')])

        mutation = OneRandomGapInsertion(probability=1.0)

        # run
        random_call.return_value = 1
        solution = mutation.execute(msa)

        # check
        self.assertEqual([('seq1', 'A-B-D')], solution.decode_alignment_as_list_of_pairs())

    @mock.patch('random.randint')
    def test_should_execute_mutation_case_b(self, random_call):
        # setup
        msa = MSASolution(self.problem, msa=[('seq1', 'AB---D'), ('seq2', 'A--B-D')])
        mutation = OneRandomGapInsertion(probability=1.0)

        # run
        random_call.return_value = 1
        solution = mutation.execute(msa)

        # check
        self.assertEqual([('seq1', 'A-B---D'), ('seq2', 'A---B-D')], solution.decode_alignment_as_list_of_pairs())

    @mock.patch('random.randint')
    def test_should_execute_mutation_case_c(self, random_call):
        # setup
        msa = MSASolution(self.problem, msa=[('seq1', 'ABBB'), ('seq2', '----')])
        mutation = OneRandomGapInsertion(probability=1.0)

        # run
        random_call.return_value = 1
        solution = mutation.execute(msa)

        # check
        self.assertEqual([('seq1', 'A-BBB'), ('seq2', '-----')], solution.decode_alignment_as_list_of_pairs())

    @mock.patch('random.randint')
    def test_should_execute_mutation_case_d(self, random_call):
        # setup
        msa = MSASolution(self.problem, msa=[('seq1', 'A-G-A-GA-G-G-T-'), ('seq2', '--CAC-A--GGT--G')])
        mutation = OneRandomGapInsertion(probability=1.0)

        # run
        random_call.return_value = 8
        solution = mutation.execute(msa)

        # check
        self.assertEqual([('seq1', 'A-G-A-GA--G-G-T-'), ('seq2', '--CAC-A---GGT--G')], solution.decode_alignment_as_list_of_pairs())


if __name__ == "__main__":
    unittest.main()
