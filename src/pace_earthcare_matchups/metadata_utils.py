"""TODO"""

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

from pace_earthcare_matchups.geospatial_utils import correct_polygon


MAAP_DT_FMT = "%Y-%m-%dT%H:%M:%SZ"
UTC = ZoneInfo("UTC")


# # TODO: convert shapely geometry to geojson somewhere else
# def geojson_from_granule(granule: Granule) -> dict:
#     """TODO"""
#     domain = granule['Granule']['Spatial']['HorizontalSpatialDomain']
#     points = domain['Geometry']['GPolygon']['Boundary']['Point']
#     coords = [[point['PointLongitude'], point['PointLatitude']] for point in points]
#     return {"type": "Feature", "geometry": {"type": "Polygon", "coordinates": [coords]}}


def polygon_from_granule(granule: Granule) -> Polygon | MultiPolygon:
    """TODO"""
    hsd = granule["Granule"]["Spatial"]["HorizontalSpatialDomain"]
    points = hsd["Geometry"]["GPolygon"]["Boundary"]["Point"]
    coords = [[point["PointLongitude"], point["PointLatitude"]] for point in points]
    coords = np.array(coords).astype(float)
    return correct_polygon(Polygon(coords))


def geometry_from_item(item: Item) -> LineString:
    """TODO"""
    # currently assume linestring
    assert item.geometry
    return LineString(item.geometry["coordinates"])


def get_intersection_bbox(granule_pace: Granule, item_earthcare: Item) -> Polygon:
    """Get latlon bbox of PACE/EarthCARE metadata intersection.
    TODO
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
    else:
        lon, lat = np.array(inter.coords.xy)
    return box(lon.min(), lat.min(), lon.max(), lat.max())


def get_datetime_range_from_granule(
    granule: Granule,
    max_offset: timedelta = timedelta(),
) -> tuple[datetime, datetime]:
    """TODO"""
    range_datetime = granule["Granule"]["Temporal"]["RangeDateTime"]
    dt_start = datetime.strptime(
        range_datetime["BeginningDateTime"],
        "%Y-%m-%dT%H:%M:%S.000Z",
    ).replace(tzinfo=UTC)
    dt_end = datetime.strptime(
        range_datetime["EndingDateTime"],
        "%Y-%m-%dT%H:%M:%S.000Z",
    ).replace(tzinfo=UTC)
    return dt_start - max_offset, dt_end + max_offset


def filter_granules(
    granules: list, max_duration: timedelta = timedelta(hours=1)
) -> list:
    """Throw away granules with broken duration metadata.
    TODO
    """
    filtered = []
    for granule in granules:
        dt_range = get_datetime_range_from_granule(granule)
        if dt_range[1] - dt_range[0] <= max_duration:
            filtered.append(granule)
    return filtered


def get_datetime_range_maap(start: datetime, end: datetime) -> str:
    """TODO"""
    return ",".join([dt.astimezone(UTC).strftime(MAAP_DT_FMT) for dt in [start, end]])
