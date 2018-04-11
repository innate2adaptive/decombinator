import unittest
import Collapsinator as coll

##First Barcode short by one
class TestFuzzyBarcodes(unittest.TestCase):

 #   perfect_positions = [22,28,36,42]
    expected_barcode = "TTTACGCGCG"

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
        seq = 'GTCGTGACTGGGAAAACCCTGGTTTAGTCGTGATCGCGCGATAATGC'
        self.run_assertions(seq, self.get_positions(seq), [22,26,34,40])

    def test_M13_spcr1_1sub(self):
        # M13 first spacer one substitution
        seq = 'GTCGTAACTGGGAAAACCCTGGTTTAGTCGTGATCGCGCGAGAATGC'
        self.run_assertions(seq, self.get_positions(seq), [22,26,34,40])

    def test_M13_spcr1_2sub(self):
        #M13 first spacer two substitution
        seq = 'TTCGTAACTGGGAAAACCCTGGTTTAGTCGTGATCGCGCGAGAATGC'
        self.run_assertions(seq, self.get_positions(seq), [22,26,34,40])

    def test_M13_spcr1_1del(self):
        #M13 first spacer one deletion
        seq = 'GTCGTGATGGGAAAACCCTGGTTTAGTCGTGATCGCGCGAGAATGC'
        self.run_assertions(seq, self.get_positions(seq), [21,25,33,39])

    def test_M13_spcr1_1ins(self):
        #M13 first spacer one insertion
        seq = 'GTTCGTGACTGGGAAAACCCTGGTTTAGTCGTGATCGCGCGAGAATGC'
        self.run_assertions(seq, self.get_positions(seq), [23,27,35,41])

    def test_M13_spcr2_1sub(self):
        #M13 second spacer one substitution
        seq = 'GTCGTGACTGGGAAAACCCTGGTTTAGTAGTGATCGCGCGAGAATGC'
        self.run_assertions(seq, self.get_positions(seq), [22,26,34,40])

    def test_M13_spcr2_2sub(self):
        #M13 second spacer two substitutions
        seq = 'GTCGTGACTGGGAAAACCCTGGTTTAGAAGTGATCGCGCGAAATGC'
        self.run_assertions(seq, self.get_positions(seq), [22,26,34,40])

    # def test_M13_spcr2_1del(self):
    #     #M13 second spacer one deletion
    #     seq = 'GTCGTGACTGGGAAAACCCTGGTTTAGTCGTGTCGCGCGAGAATGC'
    #     self.run_assertions(seq, self.get_positions(seq), [22,26,33,39])

    # def test_M13_spcr2_1ins(self):
    #     #M13 second spacer one insertion
    #     seq = 'GTCGTGACTGGGAAAACCCTGGTTTAGTCGTGTATCGCGCGAGAATGC'
    #     self.run_assertions(seq, self.get_positions(seq), [22,26,35,41])

    def test_M13_spcr1_1sub_spcr2_1sub(self):
        #M13 first spacer one substitution, second spacer one substitution
        seq = 'GTCGTGACTGGGAAACCCCTGGTTTAGTCGAGATCGCGCGATAATGC'
        self.run_assertions(seq, self.get_positions(seq), [22,26,34,40])

    def test_M13_spcr1_2sub_spcr2_2sub(self):
        #M13 first spacer two substitutions, second spacer two substitutions
        seq = 'GTCGTAACTGGGAAACCCCTGGTTTACTCGAGATCGCGCGATAATGC'
        self.run_assertions(seq, self.get_positions(seq), [22,26,34,40])

    def test_M13_spcr1_1sub_spcr2_2sub(self):
        #M13 first spacer one substitution, second spacer two substitutions
        seq = 'GTCGTGACTGGGAAACCCCTGGTTTACTCGAGATCGCGCGATAATGC'
        self.run_assertions(seq, self.get_positions(seq), [22,26,34,40])

    # def test_M13_spcr1_1del_spcr2_1ins(self):
    #     #M13 first spacer one deletion, second spacer one insertion
    #     seq = 'GTCGTGACTGGAAAACCCTGGTTTAGTCGTGATTCGCGCGATAATGC'
    #     self.run_assertions(seq, self.get_positions(seq), [21,25,34,40])

    # def test_M13_spcr1_2sub_spcr2_1del(self):
    #     #M13 first spacer two substitutions, second spacer one deletion
    #     seq = 'GTCATGACTGGGATAACCCTGGTTTATCGTGATCGCGCGATAATGC'
    #     self.run_assertions(seq, self.get_positions(seq), [22,26,33,39])

if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(TestFuzzyBarcodes)
    unittest.TextTestRunner(verbosity=2).run(suite)