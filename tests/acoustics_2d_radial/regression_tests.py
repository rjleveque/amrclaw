"""
Regression tests for 2D acoustics.
"""

import sys
import unittest

import clawpack.amrclaw.test as test


class Acoustics2DTest(test.AMRClawRegressionTest):
    r"""Basic test for a 2D acoustics test case"""


    def runTest(self, save=False):

        # Write out data files
        self.load_rundata()
        self.write_rundata_objects()

        # Run code
        self.run_code()

        # Perform tests
        self.check_gauges(save=save, 
                          regression_data_path='regression_data_test2.txt')

        self.success = True



if __name__=="__main__":
    if len(sys.argv) > 1:
        if bool(sys.argv[1]):
            # Fake the setup and save out output
            test = Acoustics2DTest()
            try:
                test.setUp()
                test.runTest(save=True)
            finally:
                test.tearDown()
            sys.exit(0)
    unittest.main()
