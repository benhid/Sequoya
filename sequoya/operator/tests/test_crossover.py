import unittest
from unittest import mock

from sequoya.core.solution import MSASolution
from sequoya.operator.crossover import SPXMSA
from sequoya.problem import MSA


class GapSequenceSolutionSinglePointTestCases(unittest.TestCase):

    def setUp(self):
        self.problem = MSA(score_list=[])
        self.problem.identifiers = ['seq1']
        self.problem.number_of_variables = 1

    def test_should_the_solution_remain_unchanged_if_the_probability_is_zero(self):
        # setup
        problem = MSA(score_list=[])
        problem.identifiers = ['seq1', 'seq2', 'seq3']
        problem.number_of_variables = 3
        msa_1 = MSASolution(problem, msa=[('seq1', 'ACTC'), ('seq2', 'A-TC'), ('seq3', 'A--C')])
        msa_2 = MSASolution(problem, msa=[('seq1', 'CT-G'), ('seq2', '-T-G'), ('seq3', '-ATG')])

        crossover = SPXMSA(probability=0.0, remove_gap_columns=False)

        # run
        offspring = crossover.execute([msa_1, msa_2])

        # check
        self.assertEqual([('seq1', 'ACTC'), ('seq2', 'A-TC'), ('seq3', 'A--C')],
                         offspring[0].decode_alignment_as_list_of_pairs())
        self.assertEqual([('seq1', 'CT-G'), ('seq2', '-T-G'), ('seq3', '-ATG')],
                         offspring[1].decode_alignment_as_list_of_pairs())

    def test_should_find_the_cutting_points_in_the_first_parent_return_the_column_position_if_it_is_occupied_by_non_gap(
            self):
        # setup
        problem = MSA(score_list=[])
        problem.identifiers = ['seq1', 'seq2']
        problem.number_of_variables = 2
        msa = MSASolution(problem, msa=[('seq1', 'BCDE'), ('seq2', 'ABCE')])

        crossover = SPXMSA(probability=1.0, remove_gap_columns=False)

        # run
        cutting_points = crossover.find_cutting_points_in_first_parent(msa, 1)

        # check
        self.assertEqual(1, cutting_points[0])
        self.assertEqual(1, cutting_points[1])

    def test_should_find_the_cutting_points_in_the_first_parent_return_the_column_position_if_it_is_occupied_by_gap(
            self):
        # setup
        problem = MSA(score_list=[])
        problem.identifiers = ['seq1', 'seq2']
        problem.number_of_variables = 2
        msa = MSASolution(problem, msa=[('seq1', 'BC-DE'), ('seq2', 'ABC-E')])

        crossover = SPXMSA(probability=1.0, remove_gap_columns=False)

        # run
        cutting_points = crossover.find_cutting_points_in_first_parent(msa, 2)

        # check
        self.assertEqual(3, cutting_points[0])
        self.assertEqual(2, cutting_points[1])

    def test_should_find_the_cutting_points_in_the_first_parent_return_minus_one_if_the_point_is_in_a_gap_group_ending_the_sequence(
            self):
        # setup
        problem = MSA(score_list=[])
        problem.identifiers = ['seq1', 'seq2']
        problem.number_of_variables = 2
        msa = MSASolution(problem, msa=[('seq1', 'BC-D-E--'), ('seq2', 'ABC-E---')])

        crossover = SPXMSA(probability=1.0, remove_gap_columns=False)

        # run
        cutting_points = crossover.find_cutting_points_in_first_parent(msa, 6)

        # check
        self.assertEqual([-1, -1], cutting_points)

    def test_should_find_original_positions_in_solution_with_gaps(self):
        # setup
        problem = MSA(score_list=[])
        problem.identifiers = ['seq1', 'seq2']
        problem.number_of_variables = 2
        msa = MSASolution(problem, msa=[('seq1', 'BC-D-E---'), ('seq2', '--C--E---')])

        crossover = SPXMSA(probability=1.0, remove_gap_columns=False)

        # run
        cutting_points = crossover.find_original_positions_in_original_sequences(msa, 5)

        # check
        self.assertEqual(3, cutting_points[0])
        self.assertEqual(1, cutting_points[1])

    def test_should_find_original_positions_in_solution_with_no_gaps(self):
        # setup
        problem = MSA(score_list=[])
        problem.identifiers = ['seq1', 'seq2']
        problem.number_of_variables = 2
        msa = MSASolution(problem, msa=[('seq1', 'ABCD'), ('seq2', 'DCBA')])

        crossover = SPXMSA(probability=1.0, remove_gap_columns=False)

        # run
        cutting_points = crossover.find_original_positions_in_original_sequences(msa, 2)

        # check
        self.assertEqual(2, cutting_points[0])
        self.assertEqual(2, cutting_points[1])

    @mock.patch('random.randint')
    def test_should_single_point_crossover_work_properly_case_a(self, random_call):
        """ AB--C|D-E, AB-C|DE- => AB--CDE-, AB-CD-E- """
        # setup
        problem = MSA(score_list=[])
        problem.identifiers = ['seq1']
        problem.number_of_variables = 1
        msa_1 = MSASolution(problem, msa=[('seq1', 'AB--CD-E')])
        msa_2 = MSASolution(problem, msa=[('seq1', 'AB--CDE-')])

        crossover = SPXMSA(probability=1.0, remove_gap_columns=False)

        # run
        random_call.return_value = 4
        children = crossover.execute([msa_1, msa_2])

        # check
        self.assertEqual(["AB--CDE-"], children[0].decode_alignment_as_list_of_sequences())
        self.assertEqual(["AB--CD-E"], children[1].decode_alignment_as_list_of_sequences())

    @mock.patch('random.randint')
    def test_should_single_point_crossover_work_properly_case_a_with_remove_gap_columns(self, random_call):
        # setup
        problem = MSA(score_list=[])
        problem.identifiers = ['seq1']
        problem.number_of_variables = 1
        msa_1 = MSASolution(problem, msa=[('seq1', 'AB--CD-E')])
        msa_2 = MSASolution(problem, msa=[('seq1', 'AB--CDE-')])

        crossover = SPXMSA(probability=1.0, remove_gap_columns=False)
        crossover_remove_full = SPXMSA(probability=1.0, remove_gap_columns=True)

        # run
        random_call.return_value = 4
        children_1 = crossover.execute([msa_1, msa_2])
        children_2 = crossover_remove_full.execute([msa_1, msa_2])

        # check
        self.assertEqual(["AB--CDE-"], children_1[0].decode_alignment_as_list_of_sequences())
        self.assertEqual(["AB--CD-E"], children_1[1].decode_alignment_as_list_of_sequences())
        self.assertEqual(["ABCDE"], children_2[0].decode_alignment_as_list_of_sequences())
        self.assertEqual(["ABCDE"], children_2[1].decode_alignment_as_list_of_sequences())

    @mock.patch('random.randint')
    def test_should_single_point_crossover_work_properly_case_b(self, random_call):
        """ A-BC|D-E-, A-B-C|DE-   => A-BCDE-, A-B-CD-E- """
        # setup
        problem = MSA(score_list=[])
        problem.identifiers = ['seq1']
        problem.number_of_variables = 1
        msa_1 = MSASolution(problem, msa=[('seq1', 'A-BCD-E')])
        msa_2 = MSASolution(problem, msa=[('seq1', 'A-B-CDE-')])

        crossover = SPXMSA(probability=1.0, remove_gap_columns=False)

        # run
        random_call.return_value = 3
        children = crossover.execute([msa_1, msa_2])

        # check
        self.assertEqual(["A-BCDE-"], children[0].decode_alignment_as_list_of_sequences())
        self.assertEqual(["A-B-CD-E"], children[1].decode_alignment_as_list_of_sequences())

    @mock.patch('random.randint')
    def test_should_single_point_crossover_work_properly_case_c(self, random_call):
        """ A|B-CD-EF, ---A|BCD-EF => ABCD-EF, ---AB-CD-EF """
        # setup
        problem = MSA(score_list=[])
        problem.identifiers = ['seq1']
        problem.number_of_variables = 1
        msa_1 = MSASolution(problem, msa=[('seq1', 'AB-CD-EF')])
        msa_2 = MSASolution(problem, msa=[('seq1', '---ABCD-EF')])

        crossover = SPXMSA(probability=1.0, remove_gap_columns=False)

        # run
        random_call.return_value = 0
        children = crossover.execute([msa_1, msa_2])

        # check
        self.assertEqual(["ABCD-EF"], children[0].decode_alignment_as_list_of_sequences())
        self.assertEqual(["---AB-CD-EF"], children[1].decode_alignment_as_list_of_sequences())

    @mock.patch('random.randint')
    def test_should_single_point_crossover_work_properly_case_d(self, random_call):
        """ GKGD---P|KK,  GKGD-P|KK  => GKGD---PKK, GKGD-PKK
            M------Q|DR,  --M--Q|DR  => M------QDR, --M--QDR """
        # setup
        problem = MSA(score_list=[])
        problem.identifiers = ['seq1', 'seq2']
        problem.number_of_variables = 2
        msa_1 = MSASolution(problem, msa=[('seq1', 'GKGD---PKK'), ('seq2', 'M------QDR')])
        msa_2 = MSASolution(problem, msa=[('seq1', 'GKGD-PKK'), ('seq1', '--M--QDR')])

        crossover = SPXMSA(probability=1.0, remove_gap_columns=False)

        # run
        random_call.return_value = 7
        children = crossover.execute([msa_1, msa_2])

        # check
        self.assertEqual(["GKGD---PKK", "M------QDR"], children[0].decode_alignment_as_list_of_sequences())
        self.assertEqual(["GKGD-PKK", "--M--QDR"], children[1].decode_alignment_as_list_of_sequences())

    @mock.patch('random.randint')
    def test_should_single_point_crossover_work_properly_case_g(self, random_call):
        """ Example of MSA in Ortuño's paper
            GKGD---PK|KP, GKGD-PK|KP   => GKGD---PK-KP, GKGD-PK--KP
            M------QD|RV, --M--QD|RV   => M------QD-RV, --M--QD--RV
            MKKLKKHPD|FP, MKKLKKHPD|FP => MKKLKKHPD-FP, MKKLKKHPDFP
            M--------|HI, ---M--H|I-   => M--------HI-, ---M--H---I """
        # setup
        problem = MSA(score_list=[])
        problem.identifiers = ['seq1', 'seq2', 'seq3', 'seq4']
        problem.number_of_variables = 4
        msa_1 = MSASolution(problem, msa=[('seq1', 'GKGD---PKKP'), ('seq2', 'M------QDRV'),
                                          ('seq3', 'MKKLKKHPDFP'), ('seq4', 'M--------HI')])
        msa_2 = MSASolution(problem, msa=[('seq1', 'GKGD-PKKP'), ('seq2', '--M--QDRV'),
                                          ('seq3', 'MKKLKKHPDFP'), ('seq4', '---M--HI-')])

        crossover = SPXMSA(probability=1.0, remove_gap_columns=False)

        # run
        random_call.return_value = 8
        children = crossover.execute([msa_1, msa_2])

        # check
        self.assertEqual(["GKGD---PK-KP", "M------QD-RV", "MKKLKKHPD-FP", "M--------HI-"],
                         children[0].decode_alignment_as_list_of_sequences())
        self.assertEqual(["GKGD-PK--KP", "--M--QD--RV", "MKKLKKHPDFP", "---M--H---I"],
                         children[1].decode_alignment_as_list_of_sequences())

    @mock.patch('random.randint')
    def test_should_single_point_crossover_work_properly_case_f(self, random_call):
        """ GKGD---P|KK, GKGD-P|KK   => GKGD---PKK, GKGD-P-KK
            M------Q|DR-, --M--Q|DR  => M------QDR, --M--QDR- """
        # setup
        problem = MSA(score_list=[])
        problem.identifiers = ['seq1', 'seq2']
        problem.number_of_variables = 2
        msa_1 = MSASolution(problem, msa=[('seq1', 'GKGD---PKK'), ('seq2', 'M------QDR-')])
        msa_2 = MSASolution(problem, msa=[('seq1', 'GKGD-PKK'), ('seq2', '--M--QDR')])

        crossover = SPXMSA(probability=1.0, remove_gap_columns=False)

        # run
        random_call.return_value = 7
        children = crossover.execute([msa_1, msa_2])

        # check
        self.assertEqual(["GKGD---PKK", "M------QDR"], children[0].decode_alignment_as_list_of_sequences())
        self.assertEqual(["GKGD-P-KK", "--M--QDR-"], children[1].decode_alignment_as_list_of_sequences())

    @mock.patch('random.randint')
    def test_should_single_point_crossover_work_properly_case_h(self, random_call):
        """ MSA with no crossover in the first sequence
            -----------|-M, --M|------  =>  ------------M------, --M """
        # setup
        problem = MSA(score_list=[])
        problem.identifiers = ['seq1']
        problem.number_of_variables = 1
        msa_1 = MSASolution(problem, msa=[('seq1', '------------M')])
        msa_2 = MSASolution(problem, msa=[('seq1', '--M------')])
        crossover = SPXMSA(probability=1.0, remove_gap_columns=False)

        # run
        random_call.return_value = 10
        children = crossover.execute([msa_1, msa_2])

        # check
        self.assertEqual(["------------M------"], children[0].decode_alignment_as_list_of_sequences())
        self.assertEqual(["--M"], children[1].decode_alignment_as_list_of_sequences())

    @mock.patch('random.randint')
    def test_should_single_point_crossover_work_properly_case_i(self, random_call):
        """ MSA with no crossover in the first sequence
            GKGD---PKKP|--, GKGD-PKKP|  =>  GKGD---PKKP--------, GKGD-PKKP--
            -----------|-M, --M|------  =>  ------------M------, --M-------- """
        # setup
        problem = MSA(score_list=[])
        problem.identifiers = ['seq1', 'seq2']
        problem.number_of_variables = 2
        msa_1 = MSASolution(problem, msa=[('seq1', 'GKGD---PKKP--'), ('seq2', '------------M')])
        msa_2 = MSASolution(problem, msa=[('seq1', 'GKGD-PKKP'), ('seq2', '--M------')])

        crossover = SPXMSA(probability=1.0, remove_gap_columns=False)

        # run
        random_call.return_value = 10
        children = crossover.execute([msa_1, msa_2])

        # check
        self.assertEqual(["GKGD---PKKP--------", "------------M------"],
                         children[0].decode_alignment_as_list_of_sequences())
        self.assertEqual(19, children[0].get_length_of_alignment())
        self.assertEqual(["GKGD-PKKP--", "--M--------"], children[1].decode_alignment_as_list_of_sequences())
        self.assertEqual(11, children[1].get_length_of_alignment())

    @mock.patch('random.randint')
    def test_should_single_point_crossover_work_properly_case_j(self, random_call):
        # setup
        problem = MSA(score_list=[])
        problem.identifiers = ['seq1', 'seq2']
        problem.number_of_variables = 2
        msa_1 = MSASolution(problem, msa=[('seq1', 'MIKMIM-IK'), ('seq2', 'A-B-CDEF-')])
        msa_2 = MSASolution(problem, msa=[('seq1', '--MIKMIMIK'), ('seq2', 'ABC-D-E-F-')])

        crossover = SPXMSA(probability=1.0, remove_gap_columns=True)

        # run
        random_call.return_value = 2
        children = crossover.execute([msa_1, msa_2])

        # check
        self.assertEqual(["MIK--MIMIK", "A-BCD-E-F-"], children[0].decode_alignment_as_list_of_sequences())
        self.assertEqual(["--MIKMIM-IK", "AB----CDEF-"], children[1].decode_alignment_as_list_of_sequences())

    @mock.patch('random.randint')
    def test_should_single_point_crossover_work_properly_real_case(self, random_call):
        # setup
        problem = MSA(score_list=[])
        problem.identifiers = ['1bbt_ac', '1al2_ad', '1b35_C', '1bbt_ab', '1mec_aa', '1bbt_aa', '1al2_ab',
                                   '1al2_ac']
        problem.number_of_variables = 8
        msa_1 = MSASolution(problem, msa=[
            ('1bbt_ac',
             '------GIFPVACSDGYGGLVTTDPKTAD---PVYGKVFNPPRNQLPGRFTNLLDVAEACP--------TFLRFEGGVPYVTTKTDSDRVLAQFDMSL----AAKHMSNTFLAG---------------------LAQYYTQYSGT-----INLHFMFTGPTDAKA-------RYMVAY----APPGMEPPKTPEAAAH---------------CIHAEWDTGLNSKF---------TFSIPYLSAADYT----YTASDVAETTNV--------QGWVCLFQ--------ITHGKADG-------DALVVLASAGKDF-----------------------ELRLPVDARAE----'),
            ('1al2_ad',
             '-------GLPVMNTPGSNQYLTADNFQSP---CALPEFDVTPPIDIPGEVKNMMELAEIDTMIPFDL--SATKKNTMEMYRVRLSDKPHTDDPILCLSLSPASDPRLSHTMLGE---------------------ILNYYTHWAGS-----LKFTFLFCGSMMATG-------KLLVSY----APPGADPPKKRKEAML---------------GTHVIWDIGLQSSC---------TMVVPWISNTT------YRQTIDDSFTE---------GGYISVFYQTRIV---VPLSTPRE-------MDILGFVSACNDF-----------------------SVRLLRDTTHIEQKA'),
            ('1b35_C',
             'SKPTVQGKIGECKLRGQGRMANFDGMDMSHKMALSSTNEIETNEGLAGTSLDVMDLSRVLSIPNYWDRFTWKTSDVINTVLWDNYVSPFKVKPYSATI-----TDRFRCTHMGK---------------------VANAFTYWRGS-----MVYTFKFVKTQYHSG---RLRISFIPYYYNTTISTGTPDVSRTQKI---------------------VVDLRTSTAV---------SFTVPYIGSRPWLYCIRPESSWLSKDNTDGALMYNCVSGIVRVEVLNQLVAAQNVFSEIDVICEVNGGPDLEFAGPTCPRY----------VPYAGDFTLADTRKIEAERTQEYSNNED'),
            ('1bbt_ab',
             '-------LLEDRILTTRNGHTTSTTQSS----VGVTYGYATAEDFVSGPNTSGLETRVV----------QAERFFKTHLFDWVTSDSFGRCHLLELPT---------DHKGVYGS--------------------LTDSYAYMRNG-----WDVEVTAVGNQFNGG-------CLLVAM----VPELCSIQKRELYQLT--------------LFPHQFINPRTNMTA---------HITVPFVGVNR------YDQYKVHKP-----------WTLVVMVVAPLTV---NTEGAPQI-------KVYANIAPTNVHV-----------------------AGEFPSKE-------'),
            ('1mec_aa',
             '------------------GVENAEKGVTEN--TDATADFVAQPVYLPENQTKVAFFYDRSSPIGRFAVKSGSLESGFAPFSNKACPNSVILTPGPQFDPAYDQLRPQRLTEIWGNGNEETSEVFPLKTKQDYSFCLFSPFVYYKCD-----LEVTLSPHTSGAHGL---------LVRW----CPTGTPTKPTTQVLHEVSSLSEGRT------PQVYSAGPGTSNQI---------SFVVPYNSPLSVLPAVWYNGHKRFDNTGD--------LGIAPNSDFGTLF---FAGTKPDI-------KFTVYLRYKNMRVFCPRP--TVFFPWPT----SGDKIDMTPRAGVL-----'),
            ('1bbt_aa',
             '---------------------TTSAGESADPVTTTVENYGGETQIQRRQHTDVSFI--------------------MDRFVKVTPQNQINILDLMQVP---------SHTLVGG---------------------LLRASTYYFSD-----LEIAVK------HEG---------DLTW----VPNGAPEK---------------------------ALDNTTNPTAYHKAPLT--RLALPYTAPHRVLATV-YNGECRTLPTSFN-------YGAIKATRVTELL---YRMKRAETYCP----RPLLAIHPTEARH---------------------KQKIVAP----------'),
            ('1al2_ab',
             '------AATSRDALPNTEASGPTHSKEIP---ALTAVETGATNPLVPSDTVQTRHVVQH----------RSRSESSIESFFARGACVTIMTVDNPAST-----TNKDKLFAVWKITYKDTVQLRR----------KLEFFTYSRFD-----MELTFVVTANFTETNNGHALNQVYQIMY----IPPGAPVP----EKWD-----------------DYTWQTSSNPSIFYTYGTAPARISVPYVGISN-AYSHFYDGFSKVPLKDQSAALGDSLYGAASLNDFGILAVRVVNDHNPTKVT----SKIRVYLKPKHIRVWCPRPPRAVAYYGPGVDYKDGTLTPLSTKDLTTY----'),
            ('1al2_ac',
             '----EACGYSDRVLQLTLGNSTITTQEA----ANSVVAYGRWPEYLRDSEANPVDQPTEPDV-------AACRFYTLDTVSWTKESRGWWWKLPDALRDMGLFGQNMYYHYLGRSGYTVHVQCNASKFHQGALGVFAVPEMCLAGDSNTTTMHTSYQNANPGEKGG-------TFTGTF----TPDNNQTSPARRFCPVDYLLGNGTLLGNAFVFPHQIINLRTNNCA---------TLVLPYVNSLS------IDSMVKHNN-----------WGIAILPLAPLNF---ASESSPEI-------PITLTIAPMCCEF-------------------NGLRNITLPRLQ-------'),
        ])
        msa_2 = MSASolution(problem, msa=[
            ('1bbt_ac',
             '------GIFPVACSDGYGGLVTTDPKTAD---PVYGKVFNPPRNQLPGRFTNLLDVAEACP--------TFLRFEGGVPYVTTKTDSDRVLAQFDMSL----AAKHMSNTFLAG---------------------LAQYYTQYSGT-----INLHFMFTGPTDAKA-------RYMVAY----APPGMEPPKTPEAAAH---------------CIHAEWDTGLNSKF---------TFSIPYLSAADYT----YTASDVAETTNV--------QGWVCLFQ--------ITHGKADG-------DALVVLASAGKDF-----------------------ELRLPVDARAE----'),
            ('1al2_ad',
             '-------GLPVMNTPGSNQYLTADNFQSP---CALPEFDVTPPIDIPGEVKNMMELAEIDTMIPFDL--SATKKNTMEMYRVRLSDKPHTDDPILCLSLSPASDPRLSHTMLGE---------------------ILNYYTHWAGS-----LKFTFLFCGSMMATG-------KLLVSY----APPGADPPKKRKEAML---------------GTHVIWDIGLQSSC---------TMVVPWISNTT------YRQTIDDSFTE---------GGYISVFYQTRIV---VPLSTPRE-------MDILGFVSACNDF-----------------------SVRLLRDTTHIEQKA'),
            ('1b35_C',
             'SKPTVQGKIGECKLRGQGRMANFDGMDMSHKMALSSTNEIETNEGLAGTSLDVMDLSRVLSIPNYWDRFTWKTSDVINTVLWDNYVSPFKVKPYSATI-----TDRFRCTHMGK---------------------VANAFTYWRGS-----MVYTFKFVKTQYHSG---RLRISFIPYYYNTTISTGTPDVSRTQKI---------------------VVDLRTSTAV---------SFTVPYIGSRPWLYCIRPESSWLSKDNTDGALMYNCVSGIVRVEVLNQLVAAQNVFSEIDVICEVNGGPDLEFAGPTCPRY----------VPYAGDFTLADTRKIEAERTQEYSNNED'),
            ('1bbt_ab',
             '-------LLEDRILTTRNGHTTSTTQSS----VGVTYGYATAEDFVSGPNTSGLETRVV----------QAERFFKTHLFDWVTSDSFGRCHLLELPT---------DHKGVYGS--------------------LTDSYAYMRNG-----WDVEVTAVGNQFNGG-------CLLVAM----VPELCSIQKRELYQLT--------------LFPHQFINPRTNMTA---------HITVPFVGVNR------YDQYKVHKP-----------WTLVVMVVAPLTV---NTEGAPQI-------KVYANIAPTNVHV-----------------------AGEFPSKE-------'),
            ('1mec_aa',
             '------------------GVENAEKGVTEN--TDATADFVAQPVYLPENQTKVAFFYDRSSPIGRFAVKSGSLESGFAPFSNKACPNSVILTPGPQFDPAYDQLRPQRLTEIWGNGNEETSEVFPLKTKQDYSFCLFSPFVYYKCD-----LEVTLSPHTSGAHGL---------LVRW----CPTGTPTKPTTQVLHEVSSLSEGRT------PQVYSAGPGTSNQI---------SFVVPYNSPLSVLPAVWYNGHKRFDNTGD--------LGIAPNSDFGTLF---FAGTKPDI-------KFTVYLRYKNMRVFCPRP--TVFFPWPT----SGDKIDMTPRAGVL-----'),
            ('1bbt_aa',
             '---------------------TTSAGESADPVTTTVENYGGETQIQRRQHTDVSFI--------------------MDRFVKVTPQNQINILDLMQVP---------SHTLVGG---------------------LLRASTYYFSD-----LEIAVK------HEG---------DLTW----VPNGAPEK---------------------------ALDNTTNPTAYHKAPLT--RLALPYTAPHRVLATV-YNGECRTLPTSFN-------YGAIKATRVTELL---YRMKRAETYCP----RPLLAIHPTEARH---------------------KQKIVAP----------'),
            ('1al2_ab',
             '------AATSRDALPNTEASGPTHSKEIP---ALTAVETGATNPLVPSDTVQTRHVVQH----------RSRSESSIESFFARGACVTIMTVDNPAST-----TNKDKLFAVWKITYKDTVQLRR----------KLEFFTYSRFD-----MELTFVVTANFTETNNGHALNQVYQIMY----IPPGAPVP----EKWD-----------------DYTWQTSSNPSIFYTYGTAPARISVPYVGISN-AYSHFYDGFSKVPLKDQSAALGDSLYGAASLNDFGILAVRVVNDHNPTKVT----SKIRVYLKPKHIRVWCPRPPRAVAYYGPGVDYKDGTLTPLSTKDLTTY----'),
            ('1al2_ac',
             '----EACGYSDRVLQLTLGNSTITTQEA----ANSVVAYGRWPEYLRDSEANPVDQPTEPDV-------AACRFYTLDTVSWTKESRGWWWKLPDALRDMGLFGQNMYYHYLGRSGYTVHVQCNASKFHQGALGVFAVPEMCLAGDSNTTTMHTSYQNANPGEKGG-------TFTGTF----TPDNNQTSPARRFCPVDYLLGNGTLLGNAFVFPHQIINLRTNNCA---------TLVLPYVNSLS------IDSMVKHNN-----------WGIAILPLAPLNF---ASESSPEI-------PITLTIAPMCCEF-------------------NGLRNITLPRLQ-------'),
        ])

        crossover = SPXMSA(probability=1.0, remove_gap_columns=False)

        # run
        random_call.return_value = 176
        children = crossover.execute([msa_1, msa_2])

        # check
        self.assertTrue(children[0].is_valid_msa())
        self.assertTrue(children[1].is_valid_msa())

    def test_should_single_point_crossover_work_properly_dummy_case(self):
        # setup
        problem = MSA(score_list=[])
        problem.identifiers = ['a', 'b', 'c', 'd']
        problem.number_of_variables = 4
        msa_1 = MSASolution(problem, msa=[
            ('a', '----GKGDPKKPRGKMSSYAFFVQTSREEHKKKHPDASVNFSEFSKKCSERWKTMSAKEKGKFEDMAKADKARYEREMKTYIPPK----------GE'),
            ('b', '-------MQDRVKRPMNAFIVWSRDQRRKMALENPRMRN--SEISKQLGYQWKMLTEAEKWPFFQEAQKLQAMHREKYPNYKYRP---RRKAKMLPK'),
            ('c', 'MKKLK---KHPDFPKKPLTPYFRFFMEKRAKYAKLHPEMSNLDLTKILSKKYKELPEKKKMKYIQDFQREKQEFERNLARFREDH---PDLIQNAKK'),
            ('d', '---------MHIKKPLNAFMLYMKEMRANVVAESTLKES--AAINQILGRRWHALSREEQAKYYELARKERQLHMQLYPGWSARDNYGKKKKRKREK')
        ])
        msa_2 = MSASolution(problem, msa=[
            ('a', '----GKGDPKKPRGKMSSYAFFVQTSREEHKKKHPDASVNFSEFSKKCSERWKTMSAKEKGKFEDMAKADKARYEREMKTYIPPK---GE-------'),
            ('b', '----M---QDRVKRPMNAFIVWSRDQRRKMALENPRMRN--SEISKQLGYQWKMLTEAEKWPFFQEAQKLQAMHREKYPNYKYRP---RRKAKMLPK'),
            ('c', 'MKKLK-KHPDFPKKPLTPYFRFFMEKRAKYAKLHPEMSN--LDLTKILSKKYKELPEKKKMKYIQDFQREKQEFERNLARFREDH---PDLIQNAKK'),
            ('d', '-------MH--IKKPLNAFMLYMKEMRANVVAESTLKES--AAINQILGRRWHALSREEQAKYYELARKERQLHMQLYPGWSARDNYGKKKKRKREK')
        ])

        crossover = SPXMSA(probability=1.0, remove_gap_columns=False)

        # run
        children = crossover.cross_parents(10, [msa_1, msa_2], [10, 10, 10, 10], [10, 10, 8, 8])

        # check
        self.assertTrue(children[0].is_valid_msa())
        self.assertTrue(children[1].is_valid_msa())

    def test_should_single_point_crossover_work_properly_real_case(self):
        # setup
        problem = MSA(score_list=[])
        problem.identifiers = ['a', 'b', 'c', 'd']
        problem.number_of_variables = 4
        msa_1 = MSASolution(problem, msa=[
            ('a', '----GKGDPKKPRGKMSSYAFFVQTSREEHKKKHPDASVNFSEFSKKCSERWKTMSAKEKGKFEDMAKADKARYEREMKTYIPPK----------GE'),
            ('b', '-------MQDRVKRPMNAFIVWSRDQRRKMALENPRMRN--SEISKQLGYQWKMLTEAEKWPFFQEAQKLQAMHREKYPNYKYRP---RRKAKMLPK'),
            ('c', 'MKKLK---KHPDFPKKPLTPYFRFFMEKRAKYAKLHPEMSNLDLTKILSKKYKELPEKKKMKYIQDFQREKQEFERNLARFREDH---PDLIQNAKK'),
            ('d', '---------MHIKKPLNAFMLYMKEMRANVVAESTLKES--AAINQILGRRWHALSREEQAKYYELARKERQLHMQLYPGWSARDNYGKKKKRKREK')
        ])
        msa_2 = MSASolution(problem, msa=[
            ('a', '----GKGDPKKPRGKMSSYAFFVQTSREEHKKKHPDASVNFSEFSKKCSERWKTMSAKEKGKFEDMAKADKARYEREMKTYIPPK---GE-------'),
            ('b', '----M---QDRVKRPMNAFIVWSRDQRRKMALENPRMRN--SEISKQLGYQWKMLTEAEKWPFFQEAQKLQAMHREKYPNYKYRP---RRKAKMLPK'),
            ('c', 'MKKLK-KHPDFPKKPLTPYFRFFMEKRAKYAKLHPEMSN--LDLTKILSKKYKELPEKKKMKYIQDFQREKQEFERNLARFREDH---PDLIQNAKK'),
            ('d', '-------MH--IKKPLNAFMLYMKEMRANVVAESTLKES--AAINQILGRRWHALSREEQAKYYELARKERQLHMQLYPGWSARDNYGKKKKRKREK')
        ])

        crossover = SPXMSA(probability=1.0, remove_gap_columns=False)

        # run
        children = crossover.cross_parents(10, [msa_1, msa_2], [10, 10, 10, 10], [10, 10, 8, 8])

        # check
        self.assertTrue(children[0].is_valid_msa())
        self.assertTrue(children[1].is_valid_msa())

    def test_should_fill_sequences_with_gaps_to_reach_the_max_sequence_length(self):
        # setup
        problem = MSA(score_list=[])
        problem.identifiers = ['a', 'b']
        problem.number_of_variables = 2
        msa_1 = MSASolution(problem, msa=[('a', '-----GE'), ('b', 'KWPFFQEAQK')])
        msa_2 = MSASolution(problem, msa=[('a', '-----GE'), ('b', 'KWPFFQEAQK')])
        msa_3 = MSASolution(problem, msa=[('a', '-'), ('b', 'ABC')])

        crossover = SPXMSA(probability=1.0, remove_gap_columns=False)

        # run
        crossover.fill_sequences_with_gaps_to_reach_the_max_sequence_length(msa_1, 10, [-1, -1])
        crossover.fill_sequences_with_gaps_to_reach_the_max_sequence_length(msa_2, 10, [-1, 5])
        crossover.fill_sequences_with_gaps_to_reach_the_max_sequence_length(msa_3, 5, [-1, 1])

        # check
        self.assertEqual(["-----G---E", "KWPFFQEAQK"], msa_1.decode_alignment_as_list_of_sequences())
        self.assertEqual(["-----G---E", "KWPFFQEAQK"], msa_2.decode_alignment_as_list_of_sequences())
        self.assertEqual(["-----", "AB--C"], msa_3.decode_alignment_as_list_of_sequences())

    def test_should_find_max_sequence_length(self):
        # setup
        problem = MSA(score_list=[])
        problem.identifiers = ['a', 'b', 'c']
        problem.number_of_variables = 3
        msa = MSASolution(problem, msa=[('a', 'AAC'), ('b', 'AAAAAAAC'), ('c', 'C')])

        crossover = SPXMSA(probability=1.0, remove_gap_columns=False)

        # run
        max = crossover.find_length_of_the_largest_sequence(msa)

        # check
        self.assertEqual(8, max)


if __name__ == "__main__":
    unittest.main()
