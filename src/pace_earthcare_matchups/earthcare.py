from dataclasses import dataclass
from datetime import datetime
from dateutil import parser
import os
from pathlib import Path
import requests

from pystac.item import Item
from tqdm import tqdm


def get_short_term_token(long_term_token: str) -> str:
    """TODO"""
    response = requests.post(
        "https://iam.maap.eo.esa.int/realms/esa-maap/protocol/openid-connect/token",
        data={
            "client_id": "offline-token",
            "client_secret": "p1eL7uonXs6MDxtGbgKdPVRAmnGxHpVE",
            "grant_type": "refresh_token",
            "refresh_token": long_term_token,
            "scope": "offline_access openid"
        }
    )
    response.raise_for_status()
    access_token = response.json().get("access_token")
    if not access_token:
        raise RuntimeError("Failed to retrieve access token from IAM response!")
    return access_token
        

def download_earthcare_item(item: Item, long_term_token: str, datadir: Path) -> Path:
    """TODO"""
    url_h5 = item.assets["enclosure_h5"].href
    # This is extremely unfortunate, but AWS bucket name length limits require
    #   the truncation of ONE character off of the end of the standard
    #   EarthCARE filenames. Therefore, truncate the unnecessary "ECA_" off of
    #   EarthCARE files.
    filename = item.assets["enclosure_h5"].title.removeprefix("ECA_")
    path_outfile = (datadir / filename)
    if path_outfile.exists():
        return path_outfile
    response = requests.get(
        url_h5,
        headers={"Authorization": f"Bearer {get_short_term_token(long_term_token)}"},
        stream=True,
    )
    response.raise_for_status()
    content_length = int(response.headers.get('content-length', 0))
    os.makedirs(datadir, exist_ok=True)
    with open(path_outfile, "wb") as outfile:
        pbar = tqdm(
            desc=f"{filename}",
            total=content_length // (2 ** 20),
            unit="MB",
        )
        # 10 MB chunks
        for chunk in response.iter_content(chunk_size=10 * (2 ** 20)):
            pbar.update(outfile.write(chunk) // (2 ** 20))
    return path_outfile


@dataclass
class EarthcareNameData:
    agency: str
    latency: str
    baseline: str
    category: str
    product: str
    level: str
    val_start: datetime
    val_end: datetime
    orbit_no: int
    frame_id: str

    def get_file_type(self) -> str:
        """TODO"""
        def _pad(s: str, plen: int = 4) -> str:
            return s + "_" * max(0, plen - len(s))
        return  _pad(self.category) + _pad(self.product) + self.level


def parse_earthcare_filename(filename: str | Path) -> EarthcareNameData:
    """TODO
    https://earthcarehandbook.earth.esa.int/article/product
    """
    stem = filename if isinstance(filename, str) else filename.stem
    stem = stem.removeprefix("ECA_")  # remove mission identifier if there
    stem_list = [s for s in stem.split("_") if s != ""]
    assert len(stem_list) == 7
    assert len(stem_list[0]) == 4
    agency, latency = stem_list[0][:2]
    baseline = stem_list[0][2:]
    category, product, level = stem_list[1:4]
    assert len(category) in [2, 3]
    assert len(product) in [2, 3]
    assert len(level) in [1, 2]
    val_start = parser.parse(stem_list[4])
    val_end = parser.parse(stem_list[5])
    orbit_no, frame_id = int(stem_list[6][:-1]), stem_list[6][-1]
    return EarthcareNameData(agency, latency, baseline, category, product,
                             level, val_start, val_end, orbit_no, frame_id)
