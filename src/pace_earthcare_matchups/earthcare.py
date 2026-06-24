"""This module handles EarthCARE tokens, downloading EarthCARE data, and parsing
filenames.
"""

from dataclasses import dataclass
from datetime import datetime, timedelta, UTC
from dateutil import parser
import os
from pathlib import Path
import requests

import h5py
import numpy as np
import numpy.typing as npt
from pystac.client import Client
from pystac.item import Item
from tqdm.notebook import tqdm

from pace_earthcare_matchups.path_utils import get_path, PATH_TOKEN


def get_short_term_token() -> str:
    """Get a short-term ESA MAAP token using your long-term token.

    :returns: A short-term ESA MAAP token.
    """
    # check for expiration time of the long-term token
    modtime = datetime.fromtimestamp(os.path.getmtime(PATH_TOKEN), UTC)
    time_until_expiration = timedelta(days=90) - (datetime.now(UTC) - modtime)
    if time_until_expiration < timedelta():
        raise RuntimeError(
            f"Long-term ESA MAAP token at {PATH_TOKEN} has expired! Please get a new "
            "long-term token from "
            "https://portal.maap.eo.esa.int/ini/services/auth/token/"
        )
    elif time_until_expiration < timedelta(days=7):
        print(
            f"Warning: Long-term ESA MAAP token at {PATH_TOKEN} will expire in "
            f"{time_until_expiration.days} days. You can get a new long-term token "
            "from https://portal.maap.eo.esa.int/ini/services/auth/token/"
        )
    long_term_token = open(PATH_TOKEN).read().rstrip("\n")
    response = requests.post(
        "https://iam.maap.eo.esa.int/realms/esa-maap/protocol/openid-connect/token",
        data={
            "client_id": "offline-token",
            "client_secret": "p1eL7uonXs6MDxtGbgKdPVRAmnGxHpVE",
            "grant_type": "refresh_token",
            "refresh_token": long_term_token,
            "scope": "offline_access openid",
        },
    )
    response.raise_for_status()
    access_token = response.json().get("access_token")
    if not access_token:
        raise RuntimeError("Failed to retrieve access token from IAM response!")
    return access_token


def download_earthcare_item(item: Item, datadir: Path) -> Path:
    """Download the EarthCARE file described by the provided STAC item.

    :param item: STAC entry corresponding to an EarthCARE file.
    :param datadir: Directory into which the EarthCARE file will be downloaded.
    :returns: Path to the downloaded EarthCARE file.
    """
    url_h5 = item.assets["enclosure_h5"].href
    # This is extremely unfortunate, but AWS bucket name length limits require
    #   the truncation of ONE character off of the end of the standard
    #   EarthCARE filenames. Therefore, truncate the unnecessary "ECA_" off of
    #   EarthCARE files.
    title = item.assets["enclosure_h5"].title
    assert title
    filename = title.removeprefix("ECA_")
    path_outfile = datadir / filename
    if path_outfile.exists():
        return path_outfile
    response = requests.get(
        url_h5,
        headers={"Authorization": f"Bearer {get_short_term_token()}"},
        stream=True,
    )
    response.raise_for_status()
    content_length = int(response.headers.get("content-length", 0))
    os.makedirs(datadir, exist_ok=True)
    with open(path_outfile, "wb") as outfile:
        pbar = tqdm(
            desc=f"{filename}",
            total=content_length // (2**20),
            unit="MB",
        )
        # 10 MB chunks
        for chunk in response.iter_content(chunk_size=10 * (2**20)):
            pbar.update(outfile.write(chunk) // (2**20))
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
        """Get the 10-character file type code of this file.

        :returns: 10-character file type code.
        """

        def _pad(s: str, plen: int = 4) -> str:
            """Pad a string with trailing underscores to a target length.

            :param s: String to pad.
            :param plen: Target length (default: 4).
            :returns: Padded string.
            """
            return s + "_" * max(0, plen - len(s))

        return _pad(self.category) + _pad(self.product) + self.level


def parse_earthcare_filename(filename: str | Path) -> EarthcareNameData:
    """Parse an EarthCARE filename or filepath.

    See https://earthcarehandbook.earth.esa.int/article/product for more details on the
    EarthCARE file naming conventions.

    :param filename: Name of or path to an EarthCARE file.
    :returns: Description of the EarthCARE file name.
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
    return EarthcareNameData(
        agency,
        latency,
        baseline,
        category,
        product,
        level,
        val_start,
        val_end,
        orbit_no,
        frame_id,
    )


def download_missing_earthcare_data(
    filepath: Path,
    client_esa: Client,
) -> None:
    """Download an EarthCARE file if it is not present at the expected local path.

    Searches for the file in the ESA STAC catalog using the filename metadata and
    downloads it to the appropriate local directory.

    :param filepath: Expected local path of the EarthCARE file.
    :param client_esa: pySTAC client to access the ESA MAAP catalog.
    """
    ec_namedata = parse_earthcare_filename(filepath)

    time_window = (
        ec_namedata.val_start + timedelta(seconds=1),
        ec_namedata.val_start + timedelta(seconds=2),
    )
    shortname = ec_namedata.category + "_" * (4 - len(ec_namedata.category))
    shortname += ec_namedata.product + "_" * (4 - len(ec_namedata.product))
    shortname += ec_namedata.level
    results = client_esa.search(
        collections=["EarthCAREL1Validated_MAAP", "EarthCAREL2Validated_MAAP"],
        datetime=time_window,
        method="GET",
        filter=f"productType = '{shortname}'",
    )
    items = list(results.items())
    assert len(items) == 1
    assert filepath.stem == items[0].id.removeprefix("ECA_")
    download_earthcare_item(
        item=items[0],
        datadir=get_path(items[0]).parent,
    )


def get_earthcare_latlon(
    filepath: Path,
) -> tuple[
    npt.NDArray[np.float32 | np.float64],
    npt.NDArray[np.float32 | np.float64],
]:
    """Get the latitude and longitude arrays from an EarthCARE file.

    :param filepath: Path to the EarthCARE file.
    :returns: Tuple of (latitude array, longitude array).
    """
    data_earthcare = h5py.File(filepath)
    science_data = data_earthcare["ScienceData"]
    assert isinstance(science_data, h5py.Group)
    if "latitude" in science_data:
        lat_earthcare = science_data["latitude"]
        lon_earthcare = science_data["longitude"]
    elif "sample_latitude" in science_data:
        lat_earthcare = science_data["sample_latitude"]
        lon_earthcare = science_data["sample_longitude"]
    elif "barycentre_latitude" in science_data:
        lat_earthcare = science_data["barycentre_latitude"]
        lon_earthcare = science_data["barycentre_longitude"]
    elif "Geo" in science_data:
        lat_earthcare = science_data["Geo/latitude"]
        lon_earthcare = science_data["Geo/longitude"]
    else:
        product_name = parse_earthcare_filename(filepath).get_file_type()
        raise NotImplementedError(
            f"Don't know how to parse EarthCARE product {product_name}!"
        )
    assert isinstance(lat_earthcare, h5py.Dataset)
    assert isinstance(lon_earthcare, h5py.Dataset)
    lat = lat_earthcare[()]
    lon = lon_earthcare[()]
    if hasattr(lat_earthcare, "fillvalue"):
        assert hasattr(lon_earthcare, "fillvalue")
        fillvalue = lat_earthcare.fillvalue
        lat[lat == fillvalue] = np.nan
        lon[lon == fillvalue] = np.nan
    return lat, lon
