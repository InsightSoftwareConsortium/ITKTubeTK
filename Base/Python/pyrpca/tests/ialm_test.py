"""ialm_test.py
"""


import os
import sys
import unittest
import numpy as np


sys.path.insert(0, os.path.join(os.path.dirname(os.path.realpath(__file__)),'..'))
import core.ialm as ialm


class IALMTesting(unittest.TestCase):
    """Testcase for the IALM method.
    """

    def setUp(self):
        """Load data.
        """
        self._data = np.genfromtxt("im_outliers.dat");


    def test_recover(self):
        """Test recovery from outliers.
        """
        # run recovery
        lr, sp, _ = ialm.recover(self._data, None)
        # load baseline (no outliers)
        baseline = np.genfromtxt("im_baseline.dat")
        # Frobenius norm between recovered mat. and baseline
        d = np.linalg.norm(np.round(lr)-baseline, ord='fro')
        self.assertTrue(np.allclose(d,0.0))


if __name__ == '__main__':
    unittest.main()
