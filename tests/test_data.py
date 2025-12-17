import unittest
from unittest.mock import patch, MagicMock
import os
import shutil
import tempfile
from pathlib import Path

# Import the module to test
# We need to make sure we can import rosesim.data even if dependencies are missing (mocking them might be needed if imports happen at top level, but they are top level now)
import io
import sys
from unittest.mock import MagicMock

# Mock dependencies that might be missing
sys.modules["astropy"] = MagicMock()
# ... (omitted unchanged lines)
sys.modules["rosesim.utils"] = MagicMock()

# Improve tqdm mock to behave like an iterator
mock_tqdm_module = MagicMock()
# When tqdm(iterable, ...) is called, return the iterable
def tqdm_side_effect(iterable, *args, **kwargs):
    return iterable
mock_tqdm_module.tqdm.side_effect = tqdm_side_effect
sys.modules["tqdm"] = mock_tqdm_module

sys.modules["requests"] = MagicMock()
sys.modules["bs4"] = MagicMock()

# Ensure environment variable is set before import to avoid __init__ error
if "ROSESIM_DATA_PATH" not in os.environ:
    os.environ["ROSESIM_DATA_PATH"] = "/tmp/rosesim_test_data_path"

# Now import valid modules
from rosesim.data import fetch_data

class TestFetchData(unittest.TestCase):
    def setUp(self):
        # Create a temporary directory for data
        self.test_dir = tempfile.mkdtemp()
        self.env_key = "ROSESIM_DATA_PATH"
        
        # Save original env var
        self.original_env_val = os.environ.get(self.env_key)

    def tearDown(self):
        # Restore env var
        if self.original_env_val:
            os.environ[self.env_key] = self.original_env_val
        elif self.env_key in os.environ:
            del os.environ[self.env_key]
            
        # Remove temp dir
        shutil.rmtree(self.test_dir)

    def test_missing_env_var(self):
        """Test that fetch_data aborts if env var is missing."""
        if self.env_key in os.environ:
            del os.environ[self.env_key]
            
        # Capture stdout to verify the error message
        with patch('sys.stdout', new=io.StringIO()) as fake_out:
            fetch_data()
            output = fake_out.getvalue()
            self.assertIn("ROSESIM_DATA_PATH environment variable is not set", output)

    @patch('rosesim.data.BeautifulSoup')
    @patch('rosesim.data.requests.get')
    @patch('rosesim.data.urlretrieve')
    def test_download_files(self, mock_urlretrieve, mock_requests_get, mock_bs):
        """Test successful download logic."""
        os.environ[self.env_key] = self.test_dir
        
        # Mock response
        mock_response = MagicMock()
        mock_response.status_code = 200
        mock_response.text = "<html></html>"
        mock_requests_get.return_value = mock_response

        # Mock BeautifulSoup
        # soup.find_all("a") -> list of mocks
        mock_soup = mock_bs.return_value
        
        def create_link_mock(href):
            m = MagicMock()
            m.get.return_value = href
            return m
            
        links = [
            create_link_mock("?C=N;O=D"),
            create_link_mock("/parent/"),
            create_link_mock("config.dat"),
            create_link_mock("image.fits"),
            create_link_mock("model.asdf"),
            create_link_mock("ignore.txt"),
            create_link_mock("subdir/"),
            create_link_mock(None),
        ]
        mock_soup.find_all.return_value = links

        # Call function
        fetch_data()

        # Check downloads
        # Expected: config.dat, image.fits, model.asdf
        # Skipped: ?..., /..., ignore.txt (not valid ext), subdir/ (ends with /), None
        
        # Helper to check if a file was downloaded
        def was_downloaded(filename):
            expected_call_found = False
            for call in mock_urlretrieve.call_args_list:
                args, _ = call
                # args[0] is url (BASE_URL + dirname + / + filename)
                # args[1] is target path
                if str(args[1]).endswith(filename):
                    expected_call_found = True
                    break
            return expected_call_found

        self.assertTrue(was_downloaded("config.dat"), "config.dat should be downloaded")
        self.assertTrue(was_downloaded("image.fits"), "image.fits should be downloaded")
        self.assertTrue(was_downloaded("model.asdf"), "model.asdf should be downloaded")
        
        self.assertFalse(was_downloaded("ignore.txt"), "ignore.txt should be skipped")
        
    @patch('rosesim.data.BeautifulSoup')
    @patch('rosesim.data.requests.get')
    @patch('rosesim.data.urlretrieve')
    @patch('rosesim.data.DIRS', ["PARSEC"])
    def test_skip_existing_files(self, mock_urlretrieve, mock_requests_get, mock_bs):
        """Test that existing files are skipped."""
        os.environ[self.env_key] = self.test_dir
        
        # data.py downloads to ROSESIM_DATA_PATH / dirname / filename
        # First dir is PARSEC
        parsec_dir = Path(self.test_dir) / "PARSEC"
        parsec_dir.mkdir(parents=True, exist_ok=True)
        existing_file = parsec_dir / "existing.fits"
        existing_file.touch()
        
        mock_response = MagicMock()
        mock_response.status_code = 200
        mock_requests_get.return_value = mock_response

        mock_soup = mock_bs.return_value
        links = [
            MagicMock(get=MagicMock(return_value="existing.fits")),
            MagicMock(get=MagicMock(return_value="new.fits")),
        ]
        mock_soup.find_all.return_value = links

        fetch_data()
        
        # Verify call parameters
        found_new = False
        found_existing = False
        
        for call in mock_urlretrieve.call_args_list:
            args, _ = call
            target_path = str(args[1])
            if target_path.endswith("/new.fits"):
                found_new = True
            if target_path.endswith("/existing.fits"):
                found_existing = True
                
        self.assertTrue(found_new, "New file should be downloaded")
        self.assertFalse(found_existing, "Existing file should be skipped")

if __name__ == "__main__":
    unittest.main()
