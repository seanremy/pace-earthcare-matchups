"""TODO"""

from datetime import datetime, timedelta
from zoneinfo import ZoneInfo

from maap.Result import Granule
import numpy as np
from pystac.item import Item
from shapely import (
    Geometry,
    GeometryCollection,
    LineString,
    MultiLineString,
    MultiPolygon,
    Polygon,
)


MAAP_DT_FMT = "%Y-%m-%dT%H:%M:%SZ"
UTC = ZoneInfo("UTC")


# # TODO: convert shapely geometry to geojson somewhere else
# def geojson_from_granule(granule: Granule) -> dict:
#     """TODO"""
#     domain = granule['Granule']['Spatial']['HorizontalSpatialDomain']
#     points = domain['Geometry']['GPolygon']['Boundary']['Point']
#     coords = [[point['PointLongitude'], point['PointLatitude']] for point in points]
#     return {"type": "Feature", "geometry": {"type": "Polygon", "coordinates": [coords]}}


def polygon_from_granule(granule: Granule) -> Polygon:
    """TODO"""
    hsd = granule['Granule']['Spatial']['HorizontalSpatialDomain']
    points = hsd['Geometry']['GPolygon']['Boundary']['Point']
    coords = [[point['PointLongitude'], point['PointLatitude']] for point in points]
    coords = np.array(coords).astype(float)
    return Polygon(coords)


def geometry_from_item(item: Item) -> Geometry:
    """TODO"""
    # currently assume linestring
    return LineString(item.geometry['coordinates'])


def get_datetime_range_from_granule(
    granule: Granule,
    max_offset: timedelta = timedelta(),
) -> (datetime, datetime):
    """TODO"""
    range_datetime = granule['Granule']['Temporal']['RangeDateTime']
    dt_start = datetime.strptime(
        range_datetime['BeginningDateTime'],
        '%Y-%m-%dT%H:%M:%S.000Z',
    ).replace(tzinfo=UTC)
    dt_end = datetime.strptime(
        range_datetime['EndingDateTime'],
        '%Y-%m-%dT%H:%M:%S.000Z',
    ).replace(tzinfo=UTC)
    return dt_start - max_offset, dt_end + max_offset


def get_datetime_range_maap(start: datetime, end: datetime) -> str:
    """TODO"""
    return ",".join([dt.astimezone(UTC).strftime(MAAP_DT_FMT) for dt in [start, end]])
