"""
Regression tests for swril 3D advection.
"""

import sys
import os
import unittest

import clawpack.amrclaw.test as test


class Advection3DSwirlTest(test.AMRClawRegressionTest):
    r"""Basic test for a 3D advection test case"""


    def runTest(self, save=False):

        # Write out data files
        self.load_rundata()
        self.write_rundata_objects()

        # Run code
        self.run_code()

        # Perform tests
        data_path = os.path.join(self.test_path, "regression_data", 
                                                 "regression_data_test2.txt")
        self.check_frame(save=save, frame_num=1, regression_data_path=data_path)
        data_path = os.path.join(self.test_path, "regression_data", 
                                                 "regression_data_test3.txt")
        self.check_frame(save=save, frame_num=2, regression_data_path=data_path)

        self.check_gauges(save=save, gauge_id=1)
        self.check_gauges(save=save, gauge_id=2)

        self.success = True



if __name__=="__main__":
    if len(sys.argv) > 1:
        if bool(sys.argv[1]):
            # Fake the setup and save out output
            test = Advection3DSwirlTest()
            try:
                test.setUp()
                test.runTest(save=True)
            finally:
                test.tearDown()
            sys.exit(0)
    unittest.main()