import unittest
from unittest import mock

from pymsa import SumOfPairs, PercentageOfTotallyConservedColumns

from pym2sa.algorithm.multiobjective.dnsgaii import reproduction
from pym2sa.core.solution import MSASolution
from pym2sa.operator import SPXMSA, OneRandomGapInsertion
from pym2sa.problem import MSA


class dNSGAIITestCases(unittest.TestCase):

    @mock.patch('random.randint')
    def test_should_reproduction_return_best_offspring(self, random_call):
        # setup
        problem = MSA(score_list=[PercentageOfTotallyConservedColumns(), SumOfPairs()],
                      original_sequences=['ABC', 'ACB'],
                      sequences_names=['seq1', 'seq2'])
        problem.obj_labels = ['%TC', 'SOP']

        msa_1 = MSASolution(problem, msa=[('seq1', '---ABC'), ('seq2', 'ACB---')])
        msa_2 = MSASolution(problem, msa=[('seq1', 'ABC-'), ('seq2', 'ACB-')])

        crossover = SPXMSA(probability=1.0, remove_gap_columns=False)
        mutation = OneRandomGapInsertion(probability=0.0, remove_gap_columns=False)

        # run
        random_call.return_value = 2
        offspring = reproduction([msa_1, msa_2], problem, crossover, mutation)

        # check
        self.assertEqual([('seq1', 'A--BC'), ('seq2', 'ACB--')], offspring.decode_alignment_as_list_of_pairs())


if __name__ == "__main__":
    unittest.main()
