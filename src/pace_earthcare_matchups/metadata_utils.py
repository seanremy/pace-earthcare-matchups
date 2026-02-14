"""Utilities for working with PACE and EarthCARE metadata."""

from datetime import datetime, timedelta

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
from pace_earthcare_matchups.pace import Granule


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
        granule_pace: A PACE granule's metadata.
        item_earthcare: STAC item describing an EarthCARE file.

    Returns:
        bbox: The bounding box of the PACE / EarthCARE metadata intersection.
    """
    geom_earthcare = geometry_from_item(item_earthcare)
    inter = granule_pace.geospatial_bounds.intersection(geom_earthcare)
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
        granule: A PACE granule's metadata.
        padding: Padding to apply to the start and end of the datetime range. Positive
            padding expands the range, negative padding shrinks the range.

    Returns:
        dt_start: Starting datetime of the range.
        dt_end: Ending datetime of the range.
    """
    if granule.beginning_datetime > granule.ending_datetime:
        raise ValueError("Start of datetime range must precede end of datetime range!")
    duration = granule.ending_datetime - granule.beginning_datetime
    if padding < timedelta() and (duration) > (-2 * padding):
        raise ValueError(
            "Negative padding must be smaller than half of the datetime range!"
        )
    return granule.beginning_datetime - padding, granule.ending_datetime + padding
