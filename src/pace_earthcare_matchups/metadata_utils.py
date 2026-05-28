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
    """Get bounding geometry from an EarthCARE STAC item (metadata).

    EarthCARE bounds may be either a line string, polygon, or multipolygon.

    :param item: STAC item describing an EarthCARE file.
    :returns: The geospatial bounds of the EarthCARE file.
    :raises ValueError: If the geometry type is not ``LineString`` or ``Polygon``.
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
    """Get the lat/lon bounding box of the PACE/EarthCARE metadata intersection.

    :param granule_pace: A PACE granule's metadata.
    :param item_earthcare: STAC item describing an EarthCARE file.
    :returns: Axis-aligned bounding box of the intersection between the PACE granule
        and EarthCARE file geometries.
    :raises TypeError: If the intersection geometry is not a supported type.
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

    :param granule: A PACE granule's metadata.
    :param padding: Padding to apply to the start and end of the datetime range.
        Positive padding expands the range, negative padding shrinks it.
    :returns: Tuple of (dt_start, dt_end), the padded start and end datetimes of the
        granule's temporal range.
    :raises ValueError: If the granule's start time is after its end time, or if
        negative padding would invert the range.
    """
    if granule.beginning_datetime > granule.ending_datetime:
        raise ValueError("Start of datetime range must precede end of datetime range!")
    duration = granule.ending_datetime - granule.beginning_datetime
    if padding < timedelta() and (duration) > (-2 * padding):
        raise ValueError(
            "Negative padding must be smaller than half of the datetime range!"
        )
    return granule.beginning_datetime - padding, granule.ending_datetime + padding
