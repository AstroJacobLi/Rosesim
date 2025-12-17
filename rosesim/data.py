from __future__ import annotations

import os
from pathlib import Path
from urllib.request import urlretrieve
from tqdm import tqdm
from bs4 import BeautifulSoup
import requests

BASE_URL = "https://tigress-web.princeton.edu/~jiaxuanl/Rosesim/"

DIRS = [
    "PARSEC",
    "TRILEGAL",
    "sky_jaguar_trilegal",
]

VALID_EXTENSIONS = (".dat", ".fits", ".asdf")

def fetch_data() -> Path:
    """
    Download required Rosesim data files if they are not already present.
    """
    data_dir = os.environ.get("ROSESIM_DATA_PATH")
    if data_dir is None:
        print("ROSESIM_DATA_PATH environment variable is not set.")
        print("It is recommended to set this variable to a persistent directory path.")
        print("Download aborted.")
        return

    data_dir = Path(data_dir).expanduser().resolve()
    data_dir.mkdir(parents=True, exist_ok=True)

    print(f"Downloading Rosesim data to {data_dir}")
    # Download directory recursively (simple index-based approach)
    for dirname in DIRS:
        _fetch_directory(dirname, data_dir)

    return data_dir


def _fetch_directory(dirname: str, data_dir: Path):
    """
    Download all files under a remote directory.
    Assumes directory listing is enabled on the server.
    """
    url = f"{BASE_URL}{dirname}/"

    r = requests.get(url)
    r.raise_for_status()

    soup = BeautifulSoup(r.text, "html.parser")

    files = []
    for link in soup.find_all("a"):
        href = link.get("href")
        if href is None:
            continue
        if href.startswith("?"):          # sorting links
            continue
        if ";" in href:                   # Apache query artifacts
            continue
        if not href.endswith(VALID_EXTENSIONS):
            continue
        files.append(href)

    desc = f"Downloading {dirname}"
    for fname in tqdm(files, desc=desc, unit="file"):
        relpath = f"{dirname}/{fname}"
        target = data_dir / relpath

        if target.exists():
            print(f"File {target} already exists, skipping.")
            continue

        target.parent.mkdir(parents=True, exist_ok=True)
        urlretrieve(f"{BASE_URL}{relpath}", target)