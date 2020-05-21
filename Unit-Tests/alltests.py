import unittest
import Collapsinator as coll

from M13tests import M13correctlength as M13c
from M13tests import M13shortby1 as M13s1
from M13tests import M13shortby2 as M13s2
from M13tests import M13longby1 as M13l1
from M13tests import M13longby2 as M13l2

from I8tests import I8correctlength as I8c
from I8tests import I8shortby1 as I8s1
from I8tests import I8shortby2 as I8s2
from I8tests import I8longby1 as I8l1
from I8tests import I8longby2 as I8l2


if __name__ == '__main__':
    
    suite1 = unittest.TestLoader().loadTestsFromTestCase(M13c.TestFuzzyBarcodes)
    suite2 = unittest.TestLoader().loadTestsFromTestCase(M13s1.TestFuzzyBarcodes)
    suite3 = unittest.TestLoader().loadTestsFromTestCase(M13s2.TestFuzzyBarcodes)
    suite4 = unittest.TestLoader().loadTestsFromTestCase(M13l1.TestFuzzyBarcodes)
    suite5 = unittest.TestLoader().loadTestsFromTestCase(M13l2.TestFuzzyBarcodes)

    suite6 = unittest.TestLoader().loadTestsFromTestCase(I8c.TestFuzzyBarcodes)
    suite7 = unittest.TestLoader().loadTestsFromTestCase(I8s1.TestFuzzyBarcodes)
    suite8 = unittest.TestLoader().loadTestsFromTestCase(I8s2.TestFuzzyBarcodes)
    suite9 = unittest.TestLoader().loadTestsFromTestCase(I8l1.TestFuzzyBarcodes)
    suite10 = unittest.TestLoader().loadTestsFromTestCase(I8l2.TestFuzzyBarcodes)

    allsuites = unittest.TestSuite([ suite1, suite2, suite3, suite4, suite5,
    								 suite6, suite7, suite8, suite9, suite10 ])

    unittest.TextTestRunner(verbosity=2).run(allsuites)