from datetime import datetime
import warnings
from zoneinfo import ZoneInfo

import earthaccess
from earthaccess.results import DataGranule
from maap.Result import Granule as MAAPGranule
from maap.maap import MAAP
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
            self._download()


def _query_cmr(
    short_name: str,
    temporal: tuple[datetime, datetime],
    limit: int,
    bbox: tuple[float, float, float, float] = (-180, -90, 180, 90),
):
    """TODO remove underscore, document"""
    use_earthaccess = bool(os.getenv("PACE_EARTHCARE_MATCHUPS_USE_EARTHACCESS"))
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
