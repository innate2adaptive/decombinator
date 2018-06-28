import unittest
import Collapsinator as coll

##Correct Length Barcodes
class TestFuzzyBarcodes(unittest.TestCase):

    #positions for perfect run are [8,14,22,28]
    expected_barcode = "TTATACGCGCG"

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
        seq = 'GTCGTGATTTATAGTCGTGATCGCGCGATAATGC'
        self.run_assertions(seq, self.get_positions(seq), [8,13,21,27])

    def test_I8_spcr1_1sub(self):
        # I8 first spacer one substitution
        seq = 'GTCGTGGTTTATAGTCGTGATCGCGCGAGAATGC'
        self.run_assertions(seq, self.get_positions(seq), [8,13,21,27])

    def test_I8_spcr1_2sub(self):
        #I8 first spacer two substitution
        seq = 'GTCGTAGTTTATAGTCGTGATCGCGCGAGAATGC'
        self.run_assertions(seq, self.get_positions(seq), [8,13,21,27])

    def test_I8_spcr1_1del(self):
        #I8 first spacer one deletion
        seq = 'GCGTGATTTATAGTCGTGATCGCGCGATAATGC'
        self.run_assertions(seq, self.get_positions(seq), [7,12,20,26])

    def test_I8_spcr1_1ins(self): 
        #I8 first spacer one insertion
        seq = 'GTACGTGATTTATAGTCGTGATCGCGCGATAATGC'
        self.run_assertions(seq, self.get_positions(seq), [9,14,22,28])

    def test_I8_spcr2_1sub(self):
        #I8 second spacer one substitution
        seq = 'GTCGTGATTTATAGTAGTGATCGCGCGAGAATGC'
        self.run_assertions(seq, self.get_positions(seq), [8,13,21,27])

    def test_I8_spcr2_2sub(self):
        #I8 second spacer two substitutions
        seq = 'GTCGTGATTTATAGAAGTGATCGCGCGAAATGC'
        self.run_assertions(seq, self.get_positions(seq), [8,13,21,27])

    # def test_I8_spcr2_1del(self):
    #     #I8 second spacer one deletion
    #     seq = 'GTCGTGATTTATAGTCGTGTCGCGCGAGAATGC'
    #     self.run_assertions(seq, self.get_positions(seq), [8,13,20,26])

    # def test_I8_spcr2_1ins(self):
    #     #I8 second spacer one insertion
    #     seq = 'GTCGTGATTTATAGTCGTGTATCGCGCGAGAATGC'
    #     self.run_assertions(seq, self.get_positions(seq), [8,13,22,28])

    def test_I8_spcr1_1sub_spcr2_1sub(self):
        #I8 first spacer one substitution, second spacer one substitution
        seq = 'GTCGTGGTTTATAGTCGTTATCGCGCGATAATGC'
        self.run_assertions(seq, self.get_positions(seq), [8,13,21,27])     

    def test_I8_spcr1_2sub_spcr2_2sub(self):
        #I8 first spacer two substitutions, second spacer two substitutions
        seq = 'GACGTGAATTATACTCGAGATCGCGCGATAATGC'
        self.run_assertions(seq, self.get_positions(seq), [8,13,21,27])

    def test_I8_spcr1_1sub_spcr2_2sub(self):
        #I8 first spacer one substitution, second spacer two substitutions
        seq = 'GTCGTAATTTATACTCGAGATCGCGCGATAATGC'
        self.run_assertions(seq, self.get_positions(seq), [8,13,21,27])

    # def test_I8_spcr1_1del_spcr2_1ins(self):
    #     #I8 first spacer one deletion, second spacer one insertion
    #     seq = 'GTCGGATTTATAGTCGTGATTCGCGCGATAATGC'
    #     self.run_assertions(seq, self.get_positions(seq), [7,12,21,27])

    # def test_I8_spcr1_2sub_spcr2_1del(self):
    #     #I8 first spacer two substitutions, second spacer one deletion
    #     seq = 'GTGGAGATTTATAGCGTGATCGCGCGATAATGC'
    #     self.run_assertions(seq, self.get_positions(seq), [8,13,20,26])


if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(TestFuzzyBarcodes)
    unittest.TextTestRunner(verbosity=2).run(suite)