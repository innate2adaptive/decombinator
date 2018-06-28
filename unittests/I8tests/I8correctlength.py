import unittest
import Collapsinator as coll

##Correct Length Barcodes
class TestFuzzyBarcodes(unittest.TestCase):

    #positions for perfect run are [8,14,22,28]
    expected_barcode = "TATATACGCGCG"

    def get_positions(self,seq):
        inputargs = {}
        inputargs['m13oligo'] = False
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

    def test_I8_perfect(self):
        #I8 perfect spacers
        seq = 'GTCGTGATTATATAGTCGTGATCGCGCGATAATGC'
        self.run_assertions(seq, self.get_positions(seq), [8,14,22,28])

    def test_I8_spcr1_1sub(self):
        # I8 first spacer one substitution
        seq = 'GTCGTGGTTATATAGTCGTGATCGCGCGAGAATGC'
        self.run_assertions(seq, self.get_positions(seq), [8,14,22,28])

    def test_I8_spcr1_2sub(self):
        #I8 first spacer two substitution
        seq = 'GTCGTAGTTATATAGTCGTGATCGCGCGAGAATGC'
        self.run_assertions(seq, self.get_positions(seq), [8,14,22,28])

    def test_I8_spcr1_1del(self):
        #I8 first spacer one deletion
        seq = 'GTGTGATTATATAGTCGTGATCGCGCGATAATGC'
        self.run_assertions(seq, self.get_positions(seq), [7,13,21,27])

    def test_I8_spcr1_1ins(self): 
        #I8 first spacer one insertion
        seq = 'GTCAGTGATTATATAGTCGTGATCGCGCGATAATGC'
        self.run_assertions(seq, self.get_positions(seq), [9,15,23,29])

    def test_I8_spcr2_1sub(self):
        #I8 second spacer one substitution
        seq = 'GTCGTGATTATATAGTAGTGATCGCGCGAGAATGC'
        self.run_assertions(seq, self.get_positions(seq), [8,14,22,28])

    def test_I8_spcr2_2sub(self):
        #I8 second spacer two substitutions
        seq = 'GTCGTGATTATATAGAAGTGATCGCGCGAAATGC'
        self.run_assertions(seq, self.get_positions(seq), [8,14,22,28])

    # def test_I8_spcr2_1del(self):
    #     #I8 second spacer one deletion
    #     seq = 'GTCGTGATTATATAGTCGTGTCGCGCGAGAATGC'
    #     self.run_assertions(seq, self.get_positions(seq), [8,14,21,27])

    # def test_I8_spcr2_1ins(self):
    #     #I8 second spacer one insertion
    #     seq = 'GTCGTGATTATATAGTCGTGTATCGCGCGAGAATGC'
    #     self.run_assertions(seq, self.get_positions(seq), [8,14,23,29])

    def test_I8_spcr1_1sub_spcr2_1sub(self):
        #I8 first spacer one substitution, second spacer one substitution
        seq = 'GTCGTGGTTATATAGTCGTTATCGCGCGATAATGC'
        self.run_assertions(seq, self.get_positions(seq), [8,14,22,28])     

    def test_I8_spcr1_2sub_spcr2_2sub(self):
        #I8 first spacer two substitutions, second spacer two substitutions
        seq = 'GACGTGAATATATACTCGAGATCGCGCGATAATGC'
        self.run_assertions(seq, self.get_positions(seq), [8,14,22,28])

    def test_I8_spcr1_1sub_spcr2_2sub(self):
        #I8 first spacer one substitution, second spacer two substitutions
        seq = 'GTCGTAATTATATACTCGAGATCGCGCGATAATGC'
        self.run_assertions(seq, self.get_positions(seq), [8,14,22,28])

    # def test_I8_spcr1_1del_spcr2_1ins(self):
    #     #I8 first spacer one deletion, second spacer one insertion
    #     seq = 'GTCTGATTATATAGTCGATGATCGCGCGATAATGC'
    #     self.run_assertions(seq, self.get_positions(seq), [7,13,22,28])

    # def test_I8_spcr1_2sub_spcr2_1del(self):
    #     #I8 first spacer two substitutions, second spacer one deletion
    #     seq = 'GTGGAGATTATATAGCGTGATCGCGCGATAATGC'
    #     self.run_assertions(seq, self.get_positions(seq), [8,14,21,27])


if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(TestFuzzyBarcodes)
    unittest.TextTestRunner(verbosity=2).run(suite)