import unittest
import Collapsinator as coll

##First Barcode short by one
class TestFuzzyBarcodes(unittest.TestCase):

 #   perfect_positions = [22,28,36,42]
    expected_barcode = "TTATACGCGCG"

    def get_positions(self,seq):
        inputargs = {}
        inputargs['m13oligo'] = True
        inputargs['positionalbarcodes'] = False
        counts = {}
        counts['getbarcode_pass_exactmatch'] = 0
        counts['getbarcode_fail_N'] = 0
        counts['getbarcode_pass_regexmatch'] = 0
        counts['getbarcode_fail_not2spacersfound'] = 0
        counts['getbarcode_fail_n1tooshort'] = 0
        counts['getbarcode_fail_n1toolong'] = 0
        counts['getbarcode_fail_n2pastend'] = 0
        counts['getbarcode_pass_fuzzymatch_rightlen'] = 0
        counts['getbarcode_pass_fuzzymatch_short'] = 0
        counts['getbarcode_pass_fuzzymatch_long'] = 0
        counts['getbarcode_pass_other'] = 0
        return coll.get_barcode(seq,inputargs,counts)

    def run_assertions(self,seq,positions,expected_positions):
        barcode = ( seq[expected_positions[0] : expected_positions[1]] + 
                    seq[expected_positions[2] : expected_positions[3]] )
        self.assertEqual(barcode, self.expected_barcode)
        self.assertEqual(positions, expected_positions)

    def test_M13_perfect(self):
        #M13 perfect spacers
        seq = 'GTCGTGACTGGGAAAACCCTGGTTATAGTCGTGATCGCGCGATAATGC'
        self.run_assertions(seq, self.get_positions(seq), [22,27,35,41])

    def test_M13_spcr1_1sub(self):
        # M13 first spacer one substitution
        seq = 'GTCGTAACTGGGAAAACCCTGGTTATAGTCGTGATCGCGCGAGAATGC'
        self.run_assertions(seq, self.get_positions(seq), [22,27,35,41])

    def test_M13_spcr1_2sub(self):
        #M13 first spacer two substitution
        seq = 'TTCGTAACTGGGAAAACCCTGGTTATAGTCGTGATCGCGCGAGAATGC'
        self.run_assertions(seq, self.get_positions(seq), [22,27,35,41])

    def test_M13_spcr1_1del(self):
        #M13 first spacer one deletion
        seq = 'GTCGTGATGGGAAAACCCTGGTTATAGTCGTGATCGCGCGAGAATGC'
        self.run_assertions(seq, self.get_positions(seq), [21,26,34,40])

    def test_M13_spcr1_1ins(self):
        #M13 first spacer one insertion
        seq = 'GTTCGTGACTGGGAAAACCCTGGTTATAGTCGTGATCGCGCGAGAATGC'
        self.run_assertions(seq, self.get_positions(seq), [23,28,36,42])

    def test_M13_spcr2_1sub(self):
        #M13 second spacer one substitution
        seq = 'GTCGTGACTGGGAAAACCCTGGTTATAGTAGTGATCGCGCGAGAATGC'
        self.run_assertions(seq, self.get_positions(seq), [22,27,35,41])

    def test_M13_spcr2_2sub(self):
        #M13 second spacer two substitutions
        seq = 'GTCGTGACTGGGAAAACCCTGGTTATAGAAGTGATCGCGCGAAATGC'
        self.run_assertions(seq, self.get_positions(seq), [22,27,35,41])

    # def test_M13_spcr2_1del(self):
    #     #M13 second spacer one deletion
    #     seq = 'GTCGTGACTGGGAAAACCCTGGTTATAGTCGTGTCGCGCGAGAATGC'
    #     self.run_assertions(seq, self.get_positions(seq), [22,27,34,40])

    # def test_M13_spcr2_1ins(self):
    #     #M13 second spacer one insertion
    #     seq = 'GTCGTGACTGGGAAAACCCTGGTTATAGTCGTGTATCGCGCGAGAATGC'
    #     self.run_assertions(seq, self.get_positions(seq), [22,27,36,42])

    def test_M13_spcr1_1sub_spcr2_1sub(self):
        #M13 first spacer one substitution, second spacer one substitution
        seq = 'GTCGTGACTGGGAAACCCCTGGTTATAGTCGAGATCGCGCGATAATGC'
        self.run_assertions(seq, self.get_positions(seq), [22,27,35,41])

    def test_M13_spcr1_2sub_spcr2_2sub(self):
        #M13 first spacer two substitutions, second spacer two substitutions
        seq = 'GTCGTAACTGGGAAACCCCTGGTTATACTCGAGATCGCGCGATAATGC'
        self.run_assertions(seq, self.get_positions(seq), [22,27,35,41])

    def test_M13_spcr1_1sub_spcr2_2sub(self):
        #M13 first spacer one substitution, second spacer two substitutions
        seq = 'GTCGTGACTGGGAAACCCCTGGTTATACTCGAGATCGCGCGATAATGC'
        self.run_assertions(seq, self.get_positions(seq), [22,27,35,41])

    # def test_M13_spcr1_1del_spcr2_1ins(self):
    #     #M13 first spacer one deletion, second spacer one insertion
    #     seq = 'GTCGTGACTGGAAAACCCTGGTTATAGTCGTGATTCGCGCGATAATGC'
    #     self.run_assertions(seq, self.get_positions(seq), [21,26,35,41])

    # def test_M13_spcr1_2sub_spcr2_1del(self):
    #     #M13 first spacer two substitutions, second spacer one deletion
    #     seq = 'GTCATGACTGGGATAACCCTGGTTATATCGTGATCGCGCGATAATGC'
    #     self.run_assertions(seq, self.get_positions(seq), [22,27,34,40])

if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(TestFuzzyBarcodes)
    unittest.TextTestRunner(verbosity=2).run(suite)