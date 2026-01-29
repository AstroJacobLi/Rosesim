import unittest
from unittest.mock import patch
import sys
import os
import shutil
from rosesim.scripts.sim_sky import main
# We need to import DATA_PATH to know where to check for output
from rosesim import DATA_PATH

class TestEmptySky(unittest.TestCase):
    def setUp(self):
        # Define the prefix used in the test
        self.prefix = 'empty_sky_test'
        self.output_dir = os.path.join(DATA_PATH, self.prefix)
        
        # Clean up before test if output directory exists from a previous run
        if os.path.exists(self.output_dir):
            shutil.rmtree(self.output_dir)

    def tearDown(self):
        # Clean up after test
        if os.path.exists(self.output_dir):
            shutil.rmtree(self.output_dir)

    def test_run_command(self):
        # The command arguments as specified by the user
        test_args = [
            "rosesim_sky",
            "--obs_ra=150.1049",
            "--obs_dec=2.2741",
            "--size=501",
            f"--prefix={self.prefix}",
            "--exptime=642",
            "--filters=['F106']",
            "--seed=42",
            "--include_bkg=False",
            "--include_star=False",
            "--psf_fov_arcsec=10"
        ]
        
        # Mock sys.argv to simulate running from command line
        with patch.object(sys, 'argv', test_args):
            # Run the main function
            main()

        # Check if the output directory was created
        self.assertTrue(os.path.isdir(self.output_dir), f"Output directory {self.output_dir} was not created")
        
        # Check if the catalog file was created
        cat_file = os.path.join(self.output_dir, 'temp', f'sky_table_{self.prefix}.ecsv')
        self.assertTrue(os.path.exists(cat_file), f"Catalog file {cat_file} was not created")

        # Check if the simulated images were created
        # The code generates files named like {band}_{exptime}s.asdf
        expected_filters = ['F106']
        exptime = 642
        for filt in expected_filters:
            img_file = os.path.join(self.output_dir, f"{filt}_{exptime}s.asdf")
            # Note: The creation of these files depends on external 'romanisim-make-l3' script execution.
            # If that script is missing or fails, these files won't exist.
            # Since the user specifically asked to test if they can run the command, 
            # verifying these files exist is the ultimate proof.
            self.assertTrue(os.path.exists(img_file), f"Image file {img_file} was not created")

if __name__ == "__main__":
    unittest.main()
