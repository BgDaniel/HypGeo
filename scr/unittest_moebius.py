import unittest
from moebius import Moeb
from helpers.group import *

NUMBER_TESTS = 1000
IDENTITY = Moeb(1.0, .0, .0, 1.0)

class TestComplexNumber(unittest.TestCase):
    def test_associative_law(self):
        A = Moeb.rnd(-50, +50, NUMBER_TESTS)
        B = Moeb.rnd(-50, +50, NUMBER_TESTS)
        C = Moeb.rnd(-50, +50, NUMBER_TESTS)

        for i in range(0, NUMBER_TESTS):
            a, b, c = A[i], B[i], C[i]    
            self.assertEqual(associative_property(a, b, c), IDENTITY)

    def test_neutral_element(self):
        M = Moeb.rnd(-50, +50, NUMBER_TESTS)

        for m in M:
            self.assertEqual(neutral_element_r(m, IDENTITY), IDENTITY)
            self.assertEqual(neutral_element_l(m, IDENTITY), IDENTITY)

    def test_inv_element(self):
        M = Moeb.rnd(-50, +50, NUMBER_TESTS)

        for m in M:
            self.assertEqual(inv_element_r(m, m.inv()), IDENTITY)
            self.assertEqual(inv_element_r(m, m.inv()), IDENTITY)


if __name__ == '__main__':
    unittest.main()

