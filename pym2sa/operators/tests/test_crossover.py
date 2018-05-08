import unittest
from unittest import mock

from pym2sa.core.solution import MSASolution
from pym2sa.operators.crossover import GapSequenceSolutionSinglePoint


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
        crossover = GapSequenceSolutionSinglePoint(probability=1.0, remove_gap_columns=False)

        # run
        offspring = crossover.execute([msa_1, msa_2])

        # check
        self.assertEqual([('seq1', 'ACTC'), ('seq2', 'A-TC'), ('seq3', 'A--C')], offspring[0].decode_alignment_as_list_of_pairs())
        self.assertEqual([('seq1', 'CT-G'), ('seq2', '-T-G'), ('seq3', '-ATG')], offspring[1].decode_alignment_as_list_of_pairs())

    def test_should_find_the_cutting_points_in_the_first_parent_return_the_column_position_if_it_is_occupied_by_non_gap(self):
        # setup
        msa = MSASolution(aligned_sequences=[('seq1', 'BCDE'), ('seq2', 'ABCE')], number_of_objectives=2)
        crossover = GapSequenceSolutionSinglePoint(probability=1.0, remove_gap_columns=False)

        # run
        cutting_points = crossover.find_cutting_points_in_first_parent(msa, 1)

        # check
        self.assertEqual(1, cutting_points[0])
        self.assertEqual(1, cutting_points[1])

    def test_should_find_the_cutting_points_in_the_first_parent_return_the_column_position_if_it_is_occupied_by_gap(self):
        # setup
        msa = MSASolution(aligned_sequences=[('seq1', 'BC-DE'), ('seq2', 'ABC-E')], number_of_objectives=2)
        crossover = GapSequenceSolutionSinglePoint(probability=1.0, remove_gap_columns=False)

        # run
        cutting_points = crossover.find_cutting_points_in_first_parent(msa, 2)

        # check
        self.assertEqual(3, cutting_points[0])
        self.assertEqual(2, cutting_points[1])

    def test_should_find_the_cutting_points_in_the_first_parent_return_minus_one_if_the_point_is_in_a_gap_group_ending_the_sequence(self):
        # setup
        msa = MSASolution(aligned_sequences=[('seq1', 'BC-D-E--'),
                                             ('seq2', 'ABC-E---')], number_of_objectives=2)
        crossover = GapSequenceSolutionSinglePoint(probability=1.0, remove_gap_columns=False)

        # run
        cutting_points = crossover.find_cutting_points_in_first_parent(msa, 6)

        # check
        self.assertEqual([-1, -1], cutting_points)

    def test_should_find_original_positions_in_solution_with_gaps(self):
        # setup
        msa = MSASolution(aligned_sequences=[('seq1', 'BC-D-E---'), ('seq2', '--C--E---')], number_of_objectives=2)
        crossover = GapSequenceSolutionSinglePoint(probability=1.0, remove_gap_columns=False)

        # run
        cutting_points = crossover.find_original_positions_in_original_sequences(msa, 5)

        # check
        self.assertEqual(3, cutting_points[0])
        self.assertEqual(1, cutting_points[1])

    def test_should_find_original_positions_in_solution_with_no_gaps(self):
        # setup
        msa = MSASolution(aligned_sequences=[('seq1', 'ABCD'), ('seq2', 'DCBA')], number_of_objectives=2)
        crossover = GapSequenceSolutionSinglePoint(probability=1.0, remove_gap_columns=False)

        # run
        cutting_points = crossover.find_original_positions_in_original_sequences(msa, 2)

        # check
        self.assertEqual(2, cutting_points[0])
        self.assertEqual(2, cutting_points[1])

    @mock.patch('random.randint')
    def test_should_single_point_crossover_work_properly_case_a(self, random_call):
        """ AB--C|D-E, AB-C|DE- => AB--CDE-, AB-CD-E- """
        # setup
        msa_1 = MSASolution(aligned_sequences=[('seq1', 'AB--CD-E')], number_of_objectives=2)
        msa_2 = MSASolution(aligned_sequences=[('seq1', 'AB--CDE-')], number_of_objectives=2)
        crossover = GapSequenceSolutionSinglePoint(probability=1.0, remove_gap_columns=False)

        # run
        random_call.return_value = 4
        children = crossover.execute([msa_1, msa_2])

        # check
        self.assertEqual(["AB--CDE-"], children[0].decode_alignment())
        self.assertEqual(["AB--CD-E"], children[1].decode_alignment())

    @mock.patch('random.randint')
    def test_should_single_point_crossover_work_properly_case_a_with_remove_gap_columns(self, random_call):
        # setup
        msa_1 = MSASolution(aligned_sequences=[('seq1', 'AB--CD-E')], number_of_objectives=2)
        msa_2 = MSASolution(aligned_sequences=[('seq1', 'AB--CDE-')], number_of_objectives=2)
        crossover = GapSequenceSolutionSinglePoint(probability=1.0, remove_gap_columns=True)

        # run
        random_call.return_value = 4
        children = crossover.execute([msa_1, msa_2])

        # check
        self.assertEqual(["ABCDE-"], children[0].decode_alignment())
        self.assertEqual(["ABCD-E"], children[1].decode_alignment())


    @mock.patch('random.randint')
    def test_should_single_point_crossover_work_properly_case_b(self, random_call):
        """ A-BC|D-E-, A-B-C|DE-   => A-BCDE-, A-B-CD-E- """
        # setup
        msa_1 = MSASolution(aligned_sequences=[('seq1', 'A-BCD-E')], number_of_objectives=2)
        msa_2 = MSASolution(aligned_sequences=[('seq1', 'A-B-CDE-')], number_of_objectives=2)
        crossover = GapSequenceSolutionSinglePoint(probability=1.0, remove_gap_columns=False)

        # run
        random_call.return_value = 3
        children = crossover.execute([msa_1, msa_2])

        # check
        self.assertEqual(["A-BCDE-"], children[0].decode_alignment())
        self.assertEqual(["A-B-CD-E"], children[1].decode_alignment())

    @mock.patch('random.randint')
    def test_should_single_point_crossover_work_properly_case_c(self, random_call):
        """ A|B-CD-EF, ---A|BCD-EF => ABCD-EF, ---AB-CD-EF """
        # setup
        msa_1 = MSASolution(aligned_sequences=[('seq1', 'AB-CD-EF')], number_of_objectives=2)
        msa_2 = MSASolution(aligned_sequences=[('seq1', '---ABCD-EF')], number_of_objectives=2)
        crossover = GapSequenceSolutionSinglePoint(probability=1.0, remove_gap_columns=False)

        # run
        random_call.return_value = 0
        children = crossover.execute([msa_1, msa_2])

        # check
        self.assertEqual(["ABCD-EF"], children[0].decode_alignment())
        self.assertEqual(["---AB-CD-EF"], children[1].decode_alignment())

    @mock.patch('random.randint')
    def test_should_single_point_crossover_work_properly_case_d(self, random_call):
        """ GKGD---P|KK,  GKGD-P|KK  => GKGD---PKK, GKGD-PKK
            M------Q|DR,  --M--Q|DR  => M------QDR, --M--QDR """
        # setup
        msa_1 = MSASolution(aligned_sequences=[('seq1', 'GKGD---PKK'), ('seq2', 'M------QDR')], number_of_objectives=2)
        msa_2 = MSASolution(aligned_sequences=[('seq1', 'GKGD-PKK'),   ('seq1', '--M--QDR')], number_of_objectives=2)
        crossover = GapSequenceSolutionSinglePoint(probability=1.0, remove_gap_columns=False)

        # run
        random_call.return_value = 7
        children = crossover.execute([msa_1, msa_2])

        # check
        self.assertEqual(["GKGD---PKK", "M------QDR"], children[0].decode_alignment())
        self.assertEqual(["GKGD-PKK",  "--M--QDR"], children[1].decode_alignment())

    @mock.patch('random.randint')
    def test_should_single_point_crossover_work_properly_case_e(self, random_call):
        """ GKGD---PK|KP, GKGD-PK|KP   => GKGD---PK-KP, GKGD-PK--KP
            M------QD|RV, --M--QD|RV   => M------QD-RV, --M--QD--RV
            MKKLKKHPD|FP, MKKLKKHPD|FP => MKKLKKHPD-FP, MKKLKKHPDFP
            M--------|HI, ---M--H|I-   => M--------HI-, ---M--H---I """
        # setup
        msa_1 = MSASolution(aligned_sequences=[('seq1', 'GKGD---PKKP'),
                                               ('seq2', 'M------QDRV'),
                                               ('seq3', 'MKKLKKHPDFP'),
                                               ('seq4', 'M--------HI')], number_of_objectives=2)
        msa_2 = MSASolution(aligned_sequences=[('seq1', 'GKGD-PKKP'),
                                               ('seq2', '--M--QDRV'),
                                               ('seq3', 'MKKLKKHPDFP'),
                                               ('seq4', '---M--HI-')], number_of_objectives=2)
        crossover = GapSequenceSolutionSinglePoint(probability=1.0, remove_gap_columns=False)

        # run
        random_call.return_value = 8
        children = crossover.execute([msa_1, msa_2])

        # check
        self.assertEqual(["GKGD---PK-KP", "M------QD-RV", "MKKLKKHPD-FP", "M--------HI-"], children[0].decode_alignment())
        self.assertEqual(["GKGD-PK--KP", "--M--QD--RV", "MKKLKKHPDFP", "---M--H---I"], children[1].decode_alignment())


if __name__ == "__main__":
    unittest.main()
