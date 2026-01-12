"""Utilities for working with PACE and EarthCARE metadata."""

from datetime import datetime, timedelta
from zoneinfo import ZoneInfo

from maap.Result import Granule
import numpy as np
from pystac.item import Item
from shapely import (
    LineString,
    MultiLineString,
    MultiPolygon,
    Polygon,
)
from shapely.geometry import box

from pace_earthcare_matchups.geospatial_utils import correct_linestring, correct_polygon


MAAP_DT_FMT = "%Y-%m-%dT%H:%M:%SZ"
UTC = ZoneInfo("UTC")


def polygon_from_granule(granule: Granule) -> Polygon | MultiPolygon:
    """Get a bounding polygon from a PACE granule (metadata).

    Args:
        granule_pace: A PACE granule's MAAP metadata.

    Returns:
        poly: A polygon or multipolygon of the granule's geospatial bounds.
    """
    hsd = granule["Granule"]["Spatial"]["HorizontalSpatialDomain"]
    points = hsd["Geometry"]["GPolygon"]["Boundary"]["Point"]
    coords = [[point["PointLongitude"], point["PointLatitude"]] for point in points]
    coords = np.array(coords).astype(float)
    return correct_polygon(Polygon(coords))


def geometry_from_item(
    item: Item,
) -> LineString | MultiLineString | Polygon | MultiPolygon:
    """Get bounding geometry from an EarthCARE STAC item (metadata). EarthCARE bounds
    may be either a linestring, polygon, or multipolygon.

    Args:
        item: STAC item describing an EarthCARE file.

    Returns:
        geom: The geospatial bounds of an EarthCARE file.
    """
    #
    assert item.geometry
    coords = np.array(item.geometry["coordinates"])
    if item.geometry["type"] == "LineString":
        return correct_linestring(LineString(coords))
    elif item.geometry["type"] == "Polygon":
        assert coords.shape[0] == 1
        return correct_polygon(Polygon(coords[0]))
    else:
        raise ValueError(f"Unrecognized geometry type: {item.geometry['type']}")


def get_intersection_bbox(granule_pace: Granule, item_earthcare: Item) -> Polygon:
    """Get latlon bbox of PACE/EarthCARE metadata intersection.

    Args:
        granule_pace: A PACE granule's MAAP metadata.
        item_earthcare: STAC item describing an EarthCARE file.

    Returns:
        bbox: The bounding box of the PACE / EarthCARE metadata intersection.
    """
    poly_pace = polygon_from_granule(granule_pace)
    geom_earthcare = geometry_from_item(item_earthcare)
    inter = poly_pace.intersection(geom_earthcare)
    if isinstance(inter, MultiLineString):
        lon, lat = np.concatenate([np.array(g.coords.xy) for g in inter.geoms], axis=-1)
    elif isinstance(inter, MultiPolygon):
        lon, lat = np.concatenate(
            [np.array(g.exterior.coords.xy) for g in inter.geoms], axis=-1
        )
    elif isinstance(inter, LineString):
        lon, lat = np.array(inter.coords.xy)
    elif isinstance(inter, Polygon):
        lon, lat = np.array(inter.exterior.coords.xy)
    else:
        raise TypeError(f"Intersection cannot be of type {type(inter)}")
    return box(lon.min(), lat.min(), lon.max(), lat.max())


def get_datetime_range_from_granule(
    granule: Granule,
    padding: timedelta = timedelta(),
) -> tuple[datetime, datetime]:
    """Get the datetime range of a PACE granule.

    Args:
        granule: A PACE granule's MAAP metadata.
        padding: Padding to apply to the start and end of the datetime range. Positive
            padding expands the range, negative padding shrinks the range.

    Returns:
        dt_start: Starting datetime of the range.
        dt_end: Ending datetime of the range.
    """
    range_datetime = granule["Granule"]["Temporal"]["RangeDateTime"]
    dt_start = datetime.strptime(
        range_datetime["BeginningDateTime"],
        "%Y-%m-%dT%H:%M:%S.000Z",
    ).replace(tzinfo=UTC)
    dt_end = datetime.strptime(
        range_datetime["EndingDateTime"],
        "%Y-%m-%dT%H:%M:%S.000Z",
    ).replace(tzinfo=UTC)
    if dt_start > dt_end:
        raise ValueError("Start of datetime range must precede end of datetime range!")
    if padding < timedelta() and (dt_end - dt_start) > (-2 * padding):
        raise ValueError(
            "Negative padding must be smaller than half of the datetime range!"
        )
    return dt_start - padding, dt_end + padding


def filter_granules(
    granules: list, max_duration: timedelta = timedelta(hours=1)
) -> list:
    """Throw away granules with broken duration metadata.

    Args:
        granules: List of PACE granules' MAAP metadata.
        max_duration: The maximum duration that granules are allowed to have.

    Returns:
        filtered: The list containing only the valid granules.
    """
    filtered = []
    for granule in granules:
        dt_range = get_datetime_range_from_granule(granule)
        if dt_range[1] - dt_range[0] <= max_duration:
            filtered.append(granule)
    return filtered


def get_datetime_range_maap(start: datetime, end: datetime) -> str:
    """Get a datetime range as a string in the NASA MAAP format.

    Args:
        start: Starting datetime of the range.
        end: Ending datetime of the range.

    Returns:
        dt_str: Datetime range as a string.
    """
    return ",".join([dt.astimezone(UTC).strftime(MAAP_DT_FMT) for dt in [start, end]])
