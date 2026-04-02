"""This module handles querying the NASA CMR for PACE data, downloading granules, and
parsing filenames.
"""

from dataclasses import dataclass
from datetime import datetime, timedelta
from dateutil import parser
from pathlib import Path
import warnings
from zoneinfo import ZoneInfo

import earthaccess
from earthaccess.results import DataGranule
from maap.Result import Granule as MAAPGranule
from maap.maap import MAAP
import netCDF4
import numpy as np
import os
from shapely import MultiPolygon, Polygon

from pace_earthcare_matchups.geospatial_utils import correct_polygon
from pace_earthcare_matchups.path_utils import PATH_DATA


CMR_HOST = "cmr.earthdata.nasa.gov"


class Granule:
    """Granule serves as a generic representation of the relevant metadata from either a
    MAAP or an earthaccess result describing a PACE granule. Both tools query the CMR,
    so the underlying metadata is similar.
    """

    def __init__(self, result: MAAPGranule | DataGranule):
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
            earthaccess.login(persist=True)
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
                earthaccess.login(persist=True),
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

    def download(self):
        if not self.filepath.exists():
            os.makedirs(self.filepath.parent, exist_ok=True)
            self._download()


def _query_cmr(
    short_name: str,
    temporal: tuple[datetime, datetime],
    limit: int,
    bbox: tuple[float, float, float, float] = (-180, -90, 180, 90),
) -> list[Granule]:
    """TODO remove underscore, document"""
    use_earthaccess = bool(int(os.getenv("PACE_EARTHCARE_MATCHUPS_USE_EARTHACCESS", "0")))
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

    Args:
        granule: A PACE granule's MAAP metadata.
        shortname_pace: The shortname of a different PACE product to retrieve.

    Returns:
        result: MAAP metadata of a PACE granule co-occurring with the provided granule.
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

    Args:
        data_pace: A HARP2 L1B data file.

    Returns:
        idx_nadir: Index of the smallest absolute viewing angle.
    """
    view_angle = data_pace["sensor_views_bands/sensor_view_angle"]
    idx_nadir = np.argmin(np.abs(view_angle))
    return int(idx_nadir.item())


def get_pace_shortname(instrument: str, level: str) -> str:
    """TODO"""
    shortname_pace = f"PACE_{instrument}_{level}"
    if level[1] == "1":
        shortname_pace += "_SCI"
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

    Args:
        filename: Name of or path to a PACE file.

    Returns:
        pace_namedata: Description of the PACE file name.
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
                raise ValueError(f"Got conflicting values {version} and {s} for version!")
            version = s
        else:
            if isinstance(product, str):
                raise ValueError(f"Got conflicting values {product} and {s} for product!")
            product = s
    return PaceNameData(
        instrument,
        start_time,
        level,
        product,
        version,
    )


def download_missing_pace_data(filepath: Path) -> None:
    pace_namedata = parse_pace_filename(filepath)
    shortname = get_pace_shortname(pace_namedata.instrument, pace_namedata.level)
    time_window=(
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
