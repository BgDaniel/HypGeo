import unittest
from complex_plane import _i, ComplexNumber, RootOfUnity
from transfomations import *

NUMBER_TESTS = 1000
ZERO = ComplexNumber(.0, .0)

class TestTransfomations(unittest.TestCase):
    def test_moeb_to(self):
        U = ComplexNumber.rnd(.0, +50, NUMBER_TESTS)
        V = ComplexNumber.rnd(.0, +50, NUMBER_TESTS)

        for i in range(0, NUMBER_TESTS):
            u, v = U[i], V[i]
            self.assertEqual(moeb_to(u, v)(u) - v, ZERO)

    def test_refl_0(self):
        R = np.random.uniform(.0, 20.0, NUMBER_TESTS)

        for r in np.random.uniform(.0, 20.0, NUMBER_TESTS):
            self.assertEqual(refl_0(r) * refl_0(r), moeb_id)
            self.assertEqual(refl_0(r) * refl_0(r), moeb_id)





       
   
if __name__ == '__main__':
    unittest.main()