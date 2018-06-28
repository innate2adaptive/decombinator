import unittest
import Collapsinator as coll

##Correct Length Barcodes
class TestFuzzyBarcodes(unittest.TestCase):

    #positions for perfect run are [22,28,36,42]
    expected_barcode = "TATATATACGCGCG"

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
        seq = 'GTCGTGACTGGGAAAACCCTGGTATATATAGTCGTGATCGCGCGATAATGC'
        self.run_assertions(seq, self.get_positions(seq), [22,30,38,44])

    def test_M13_spcr1_1sub(self):
        # M13 first spacer one substitution
        seq = 'GTCGTAACTGGGAAAACCCTGGTATATATAGTCGTGATCGCGCGAGAATGC'
        self.run_assertions(seq, self.get_positions(seq), [22,30,38,44])

    def test_M13_spcr1_2sub(self):
        #M13 first spacer two substitution
        seq = 'TTCGTAACTGGGAAAACCCTGGTATATATAGTCGTGATCGCGCGAGAATGC'
        self.run_assertions(seq, self.get_positions(seq), [22,30,38,44])

    def test_M13_spcr1_1del(self):
        #M13 first spacer one deletion
        seq = 'GTCGTGATGGGAAAACCCTGGTATATATAGTCGTGATCGCGCGAGAATGC'
        self.run_assertions(seq, self.get_positions(seq), [21,29,37,43])

    def test_M13_spcr1_1ins(self): 
        #M13 first spacer one insertion
        seq = 'GTTCGTGACTGGGAAAACCCTGGTATATATAGTCGTGATCGCGCGAGAATGC'
        self.run_assertions(seq, self.get_positions(seq), [23,31,39,45])

    def test_M13_spcr2_1sub(self):
        #M13 second spacer one substitution
        seq = 'GTCGTGACTGGGAAAACCCTGGTATATATAGTAGTGATCGCGCGAGAATGC'
        self.run_assertions(seq, self.get_positions(seq), [22,30,38,44])

    def test_M13_spcr2_2sub(self):
        #M13 second spacer two substitutions
        seq = 'GTCGTGACTGGGAAAACCCTGGTATATATAGAAGTGATCGCGCGAAATGC'
        self.run_assertions(seq, self.get_positions(seq), [22,30,38,44])

    # def test_M13_spcr2_1del(self):
    #     #M13 second spacer one deletion
    #     seq = 'GTCGTGACTGGGAAAACCCTGGTATATATAGTCGTGTCGCGCGAGAATGC'
    #     self.run_assertions(seq, self.get_positions(seq), [22,30,37,43])

    # def test_M13_spcr2_1ins(self):
    #     #M13 second spacer one insertion
    #     seq = 'GTCGTGACTGGGAAAACCCTGGTATATATAGTCGTGTATCGCGCGAGAATGC'
    #     self.run_assertions(seq, self.get_positions(seq), [22,30,39,45])

    def test_M13_spcr1_1sub_spcr2_1sub(self):
        #M13 first spacer one substitution, second spacer one substitution
        seq = 'GTCGAGACTGGGAAAACCCTGGTATATATAGTCGTTATCGCGCGATAATGC'
        self.run_assertions(seq, self.get_positions(seq), [22,30,38,44])     

    def test_M13_spcr1_2sub_spcr2_2sub(self):
        #M13 first spacer two substitutions, second spacer two substitutions
        seq = 'GTCGTAACTGGGAAACCCCTGGTATATATACTCGAGATCGCGCGATAATGC'
        self.run_assertions(seq, self.get_positions(seq), [22,30,38,44])

    def test_M13_spcr1_1sub_spcr2_2sub(self):
        #M13 first spacer one substitution, second spacer two substitutions
        seq = 'GTCGTGACTGGGAAACCCCTGGTATATATACTCGAGATCGCGCGATAATGC'
        self.run_assertions(seq, self.get_positions(seq), [22,30,38,44])

    # def test_M13_spcr1_1del_spcr2_1ins(self):
    #     #M13 first spacer one deletion, second spacer one insertion
    #     seq = 'GTCGTGACTGGAAAACCCTGGTATATATAGTCGTGATTCGCGCGATAATGC'
    #     self.run_assertions(seq, self.get_positions(seq), [21,29,38,44])

    # def test_M13_spcr1_2sub_spcr2_1del(self):
    #     #M13 first spacer two substitutions, second spacer one deletion
    #     seq = 'GTCGTGACTGGGCATACCCTGGTATATATAGCGTGATCGCGCGATAATGC'
    #     self.run_assertions(seq, self.get_positions(seq), [22,30,37,43])


if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(TestFuzzyBarcodes)
    unittest.TextTestRunner(verbosity=2).run(suite)