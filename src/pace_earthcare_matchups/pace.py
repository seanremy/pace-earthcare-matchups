"""This module handles querying the NASA CMR for PACE data, downloading granules, and
parsing filenames.
"""

from dataclasses import dataclass
from datetime import datetime, timedelta
from dateutil import parser
import os
from pathlib import Path
import warnings
from zoneinfo import ZoneInfo

import earthaccess
from earthaccess.results import DataGranule
from maap.Result import Granule as MAAPGranule
from maap.maap import MAAP
import netCDF4
import numpy as np
from shapely import MultiPolygon, Polygon

from pace_earthcare_matchups.geospatial_utils import correct_polygon
from pace_earthcare_matchups.path_utils import PATH_DATA


CMR_HOST = "cmr.earthdata.nasa.gov"
SHORT_NAME_REPLACEMENTS = {
    "CLD": "CLOUD",
    "CLDMASK": "CLOUD_MASK",
}


def _earthaccess_login():
    """Log in to earthaccess, using a MAAP token if running in a MAAP environment."""
    if "MAAP_PGT" in os.environ:
        maap = MAAP(maap_host="api.maap-project.org")
        acc_info = maap.profile.account_info()
        assert isinstance(acc_info, dict)
        os.environ["EARTHDATA_TOKEN"] = acc_info["urs_token"]
    earthaccess.login(persist=True)


class Granule:
    """Generic representation of the relevant metadata from either a MAAP or an
    earthaccess result describing a PACE granule.

    Both tools query the CMR, so the underlying metadata is similar. After construction,
    instances expose: ``short_name``, ``beginning_datetime``, ``ending_datetime``,
    ``filepath``, and ``geospatial_bounds``.
    """

    def __init__(self, result: MAAPGranule | DataGranule):
        """Initialize a Granule from a MAAP or earthaccess result.

        :param result: A granule result from either ``maap.Result.Granule`` or
            ``earthaccess.results.DataGranule``.
        :raises TypeError: If ``result`` is neither a ``MAAPGranule`` nor a
            ``DataGranule``.
        """
        if isinstance(result, MAAPGranule):
            # short name
            self.short_name = result["Granule"]["Collection"]["ShortName"]
            # beginning and ending datetime
            range_datetime = result["Granule"]["Temporal"]["RangeDateTime"]
            self.beginning_datetime = datetime.strptime(
                range_datetime["BeginningDateTime"],
                "%Y-%m-%dT%H:%M:%S.000Z",
            ).replace(tzinfo=ZoneInfo("UTC"))
            self.ending_datetime = datetime.strptime(
                range_datetime["EndingDateTime"],
                "%Y-%m-%dT%H:%M:%S.000Z",
            ).replace(tzinfo=ZoneInfo("UTC"))
            # file location
            gran = result["Granule"]
            plat = gran["Platforms"]["Platform"]
            instrument = plat["Instruments"]["Instrument"]["ShortName"]
            level = gran["Collection"]["ShortName"].split("_")[2]
            filename = gran["DataGranule"]["ProducerGranuleId"]
            self.filepath = PATH_DATA / "PACE" / instrument / level / filename
            # geospatial bounds
            hsd = result["Granule"]["Spatial"]["HorizontalSpatialDomain"]
            points = hsd["Geometry"]["GPolygon"]["Boundary"]["Point"]
            coords = [[p["PointLongitude"], p["PointLatitude"]] for p in points]
            coords = np.array(coords).astype(float)
            poly = Polygon(coords)
            try:
                self.geospatial_bounds = correct_polygon(poly)
            except ValueError:
                warnings.warn(f"Broken geospatial bounds in {filename}")
                self.geospatial_bounds = poly
            # download function
            self._download = lambda: result.getData(str(self.filepath.parent))
        elif isinstance(result, DataGranule):
            _earthaccess_login()
            # short name
            self.short_name = result["umm"]["CollectionReference"]["ShortName"]
            # beginning and ending datetime
            range_datetime = result["umm"]["TemporalExtent"]["RangeDateTime"]
            self.beginning_datetime = datetime.strptime(
                range_datetime["BeginningDateTime"],
                "%Y-%m-%dT%H:%M:%SZ",
            ).replace(tzinfo=ZoneInfo("UTC"))
            self.ending_datetime = datetime.strptime(
                range_datetime["EndingDateTime"],
                "%Y-%m-%dT%H:%M:%SZ",
            ).replace(tzinfo=ZoneInfo("UTC"))
            # file location
            instrument = result["umm"]["Platforms"][0]["Instruments"][0]["ShortName"]
            level = result["umm"]["CollectionReference"]["ShortName"].split("_")[2]
            filename = result["umm"]["DataGranule"]["Identifiers"][0]["Identifier"]
            self.filepath = PATH_DATA / "PACE" / instrument / level / filename
            # geospatial bounds
            hsd = result["umm"]["SpatialExtent"]["HorizontalSpatialDomain"]
            gpolys = hsd["Geometry"]["GPolygons"]
            polys = []
            for gpoly in gpolys:
                points = gpoly["Boundary"]["Points"]
                coords = [[p["Longitude"], p["Latitude"]] for p in points]
                coords = np.array(coords).astype(float)
                polys.append(Polygon(coords))
            if len(polys) == 1:
                try:
                    self.geospatial_bounds = correct_polygon(polys[0])
                except ValueError:
                    warnings.warn(f"Broken geospatial bounds in {filename}")
                    self.geospatial_bounds = polys[0]
            else:
                self.geospatial_bounds = MultiPolygon(polys)
            # download function
            self._download = lambda: (
                _earthaccess_login(),
                earthaccess.download(
                    [result],
                    self.filepath.parent,
                    show_progress=False,
                ),
            )
        else:
            raise TypeError(
                f"Result of type {type(result)} is neither `maap.Result.Granule` nor `earthaccess.results.DataGranule`!"
            )

    def download(self) -> None:
        """Download this granule's data file to its expected local path.

        Does nothing if the file already exists.
        """
        if not self.filepath.exists():
            os.makedirs(self.filepath.parent, exist_ok=True)
            self._download()


def _query_cmr(
    short_name: str,
    temporal: tuple[datetime, datetime],
    limit: int,
    bbox: tuple[float, float, float, float] = (-180, -90, 180, 90),
) -> list[Granule]:
    """Query the NASA CMR for PACE granules matching the given filters.

    Uses MAAP by default, or earthaccess if the environment variable
    ``PACE_EARTHCARE_MATCHUPS_USE_EARTHACCESS`` is set to ``"1"``.

    :param short_name: PACE collection short name.
    :param temporal: Time range as a (start, end) tuple of datetimes.
    :param limit: Maximum number of granules to return.
    :param bbox: Lat/lon bounding box in W, S, E, N order.
    :returns: List of Granule objects matching the search criteria.
    """
    use_earthaccess = bool(
        int(os.getenv("PACE_EARTHCARE_MATCHUPS_USE_EARTHACCESS", "0"))
    )
    if not use_earthaccess:
        bbox_str = ",".join([str(n) for n in bbox])
        results_pace = MAAP().searchGranule(
            cmr_host=CMR_HOST,
            short_name=short_name,
            temporal=",".join(
                [
                    temporal[0].strftime("%Y-%m-%dT%H:%M:%SZ"),
                    temporal[1].strftime("%Y-%m-%dT%H:%M:%SZ"),
                ]
            ),
            bounding_box=bbox_str,
            limit=limit,
        )
    else:
        results_pace = earthaccess.search_data(
            short_name=short_name,
            temporal=temporal,
            bounding_box=bbox,
            count=limit,
        )
    return [Granule(r) for r in results_pace]


def get_simultaneous_pace_product(granule: Granule, shortname_pace: str) -> Granule:
    """Using one PACE granule, get a different product with the same timestamp.

    :param granule: A PACE granule's metadata.
    :param shortname_pace: The short name of a different PACE product to retrieve.
    :returns: Metadata of a PACE granule co-occurring with the provided granule.
    """
    result = _query_cmr(
        short_name=shortname_pace,
        temporal=(granule.beginning_datetime, granule.ending_datetime),
        limit=1,
    )[0]
    time_diff = (
        abs(result.beginning_datetime - granule.beginning_datetime),
        abs(result.ending_datetime - granule.ending_datetime),
    )
    assert time_diff[0].total_seconds() < 10 and time_diff[1].total_seconds() < 10
    return result


def get_nadir_idx_harp2_l1b(data_pace: netCDF4.Dataset) -> int:
    """Get the index of the smallest absolute viewing angle in HARP2 L1B data.

    :param data_pace: An open HARP2 L1B netCDF4 dataset.
    :returns: Index into the sensor view angle dimension corresponding to the nadir
        (near-zero) viewing angle.
    """
    view_angle = data_pace["sensor_views_bands/sensor_view_angle"]
    idx_nadir = np.argmin(np.abs(view_angle))
    return int(idx_nadir.item())


def get_pace_shortname(instrument: str, level: str, filestem: str) -> str:
    """Construct a PACE collection short name from instrument, level, and file stem.

    :param instrument: Instrument identifier (e.g., ``"OCI"``, ``"HARP2"``).
    :param level: Processing level string (e.g., ``"L1B"``, ``"L2"``).
    :param filestem: File stem (without extension) of the PACE granule.
    :returns: PACE collection short name (e.g., ``"PACE_OCI_L2_AOP"``).
    :raises ValueError: If the level is not ``"L1"`` or ``"L2"``.
    """
    shortname_pace = f"PACE_{instrument}_{level}"
    if level[1] == "1":
        shortname_pace += "_SCI"
    elif level[1] == "2":
        prod = filestem.split(".")[3]
        if prod in SHORT_NAME_REPLACEMENTS:
            prod = SHORT_NAME_REPLACEMENTS[prod]
        shortname_pace += "_" + prod
    else:
        raise ValueError(f"Level must be 'L1' or 'L2', but got: '{level}'")
    if filestem.split(".")[-1] == "NRT":
        shortname_pace += "_NRT"
    return shortname_pace


@dataclass
class PaceNameData:
    instrument: str
    start_time: datetime
    level: str
    product: str | None
    version: str | None


def parse_pace_filename(filename: str | Path) -> PaceNameData:
    """Parse a PACE filename or filepath.

    :param filename: Name of or path to a PACE file.
    :returns: Parsed components of the PACE file name.
    :raises ValueError: If conflicting product or version fields are found in the stem.
    """
    stem = filename if isinstance(filename, str) else filename.stem
    stem_list = [s for s in stem.split(".") if s != ""]
    instrument = stem_list[0].split("_")[1]
    start_time = parser.parse(f"{stem_list[1]}Z")  # Z to enforce UTC
    level = stem_list[2]
    assert level in ["L1B", "L1C", "L2"]
    product, version = None, None
    for s in stem_list[3:]:
        if s.startswith("V") and s[1:].replace("_", "").isnumeric():
            if isinstance(version, str):
                raise ValueError(
                    f"Got conflicting values {version} and {s} for version!"
                )
            version = s
        else:
            if isinstance(product, str):
                raise ValueError(
                    f"Got conflicting values {product} and {s} for product!"
                )
            product = s
    return PaceNameData(
        instrument,
        start_time,
        level,
        product,
        version,
    )


def download_missing_pace_data(filepath: Path) -> None:
    """Download a PACE file if it is not present at the expected local path.

    Searches for the file in the NASA CMR using the filename metadata and downloads
    it to the appropriate local directory.

    :param filepath: Expected local path of the PACE file.
    """
    pace_namedata = parse_pace_filename(filepath)
    shortname = get_pace_shortname(
        pace_namedata.instrument, pace_namedata.level, filepath.stem
    )
    time_window = (
        pace_namedata.start_time + timedelta(seconds=1),
        pace_namedata.start_time + timedelta(seconds=2),
    )

    results = _query_cmr(
        short_name=shortname,
        temporal=time_window,
        limit=1,
    )
    assert len(results) == 1
    assert results[0].filepath.name == filepath.name
    results[0].download()


def get_pace_latlon(
    filepath: Path,
) -> tuple[
    npt.NDArray[np.float32 | np.float64],
    npt.NDArray[np.float32 | np.float64],
]:
    """Get the latitude and longitude arrays from a PACE file.

    Note: In the case of HARP2 L1B, retrieves only the nadir lat/lon.

    :param filepath: Path to the PACE file.
    :returns: Tuple of (latitude array, longitude array).
    """
    data_pace = netCDF4.Dataset(filepath)
    subset = ()
    if data_pace.instrument == "HARP2" and data_pace.processing_level == "L1B":
        subset = get_nadir_idx_harp2_l1b(data_pace)
    for group in ["geolocation_data", "navigation_data"]:
        try:
            lat = data_pace[group + "/latitude"][subset].filled(fill_value=np.nan)
            lon = data_pace[group + "/longitude"][subset].filled(fill_value=np.nan)
            return lat, lon
        except KeyError:
            continue
    raise ValueError(
        "Provided file contained neither `geolocation_data` nor `navigation_data`"
    )
