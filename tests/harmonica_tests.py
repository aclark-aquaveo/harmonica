"""For testing."""

# 1. Standard python libraries
import os
import unittest

# 2. Third party libraries

# 3. Local libraries

__copyright__ = "(C) Copyright Aquaveo 2019"
__license__ = "All rights reserved"


WINDOWS_CI_TEST_DATA_DIR = 'C:/temp/harmonica_ci_test_data'


def ensure_ci_test_data_exists():
    from harmonica import config
    config['pre_existing_data_dir'] = WINDOWS_CI_TEST_DATA_DIR
    if not os.path.exists(WINDOWS_CI_TEST_DATA_DIR):
        # Download the resources to this test machine.
        pass


# TODO: Implement these tests once we have the Github repo mirrored to GitLab to use local runners.
class HarmonicaTests(unittest.TestCase):
    """
    Tests the tidal data class.
    """

    @classmethod
    def setUpClass(cls):
        # Change working directory to test location
        os.chdir(os.path.normpath(os.path.join(os.path.abspath(__file__), os.pardir)))
        # Make sure data files are present on the test machine
        ensure_ci_test_data_exists()

    def test_nodal_factor(self):
        """Test extracting constituent properties for all constituents"""
        pass

    def test_adcirc(self):
        """Test tidal extraction for the ADCIRC 2015 model"""
        pass

    def test_leprovost(self):
        """Test tidal extraction for the legacy LeProvost model"""
        pass

    def test_fes2014(self):
        """Test tidal extraction for the FES2014 model"""
        pass

    def test_tpxo8(self):
        """Test tidal extraction for the TPXO8 model"""
        pass

    def test_tpxo9(self):
        """Test tidal extraction for the TPXO9 model"""
        pass
