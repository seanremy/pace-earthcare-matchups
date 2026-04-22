"""Utilities for plotting matchups."""

from pathlib import Path

import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.geoaxes import GeoAxes
import matplotlib.lines
import matplotlib.patches
import matplotlib.pyplot as plt
import numpy as np
import numpy.typing as npt
import seaborn as sns
from shapely import (
    LineString,
    MultiLineString,
    MultiPolygon,
    Polygon,
)

from pace_earthcare_matchups.geospatial_utils import geom_to_coords
from pace_earthcare_matchups.matchup import Matchup


def get_best_longitude_shift(matchups: list[Matchup]) -> float:
    """For a list of matchups, compute the longitude shift which 
    """
    # flatten all geometries into a single array of longitudes
    lon = []
    for matchup in matchups:
        bounds_pace = matchup.get_pace_bounds()
        lon.append(geom_to_coords(bounds_pace)[:, 0].ravel())
        for match in matchup.matches_earthcare:
            bounds_earthcare = match.get_earthcare_bounds()
            lon.append(geom_to_coords(bounds_earthcare)[:, 0].ravel())
    lon = np.concatenate(lon)
    assert lon.min() >= -180 and lon.max() <= 180
    # sort the longitudes
    lon = np.sort(lon)
    # get the differences in longitudes, accounting for the wrap
    diffs = np.concatenate([np.diff(lon), [lon[0] + 360 - lon[-1]]])
    # find the largest jump
    idx_jump = np.argmax(diffs)
    # determine the antemeridian
    if idx_jump == diffs.shape[0] - 1:
        # "edge case" where 180 already separates the data
        antemeridian = (lon[-1] + lon[0]) / 2 % 360 - 180
    else:
        antemeridian = (lon[idx_jump] + lon[idx_jump + 1]) / 2
    lon_shift = antemeridian % 360 - 180
    return lon_shift.item()


def plot_matchups(
    matchups: list[Matchup],
    figsize: tuple[int, int] | None = None,
    fig_filepath: Path | str | None = None,
) -> None:
    """Plot a list of matchups and optionally save the figure to file.

    Args:
        matchups: List of matchups to be plotted.
        figsize: Size of the matplotlib figure as an (x, y) tuple.
        fig_filepath: If provided, saves the figure to this path. Should end in ".png".
        central_longitude: Central longitude of the created plot. Default: 0.
    """
    if isinstance(fig_filepath, str):
        fig_filepath = Path(fig_filepath)
    lon_shift = get_best_longitude_shift(matchups)
    ax = plt.figure(figsize=figsize).add_subplot(
        projection=ccrs.PlateCarree(central_longitude=lon_shift)
    )
    assert isinstance(ax, GeoAxes)
    ax.stock_img(alpha=0.8)
    ax.add_feature(cfeature.COASTLINE)

    plot_elements = {}
    # first, assign a color palette to each instrument
    labels_pace = set()
    labels_earthcare = set()
    for matchup in matchups:
        labels_pace.add(matchup.shortname_pace)
        for match in matchup.matches_earthcare:
            labels_earthcare.add(
                f"EarthCARE_{match.filepath_earthcare.parts[-2].split('_')[0]}"
            )
    palette = dict(
        zip(
            list(labels_pace) + list(labels_earthcare),
            sns.color_palette("colorblind", len(labels_pace) + len(labels_earthcare)),
        )
    )
    # keep track of min/max lat/lon for cropping the plot later
    minlat, maxlat = np.inf, -np.inf
    minlon, maxlon = np.inf, -np.inf

    def _update_extent(coords: npt.NDArray):
        nonlocal minlat, maxlat, minlon, maxlon
        minlat = min(minlat, np.min(coords[..., 1]))
        maxlat = max(maxlat, np.max(coords[..., 1]))
        minlon = min(minlon, np.min(coords[..., 0]))
        maxlon = max(maxlon, np.max(coords[..., 0]))

    # helper function for plotting the bounds of either PACE or EarthCARE data
    def _plot_bounds(
        bounds: LineString | MultiLineString | Polygon | MultiPolygon, label: str
    ) -> matplotlib.lines.Line2D | matplotlib.patches.Polygon:
        if isinstance(bounds, MultiLineString | MultiPolygon):
            geoms = bounds.geoms
        else:
            geoms = [bounds]
        plot_element = None
        for geom in geoms:
            if isinstance(geom, LineString):
                coords_geom = np.array(geom.coords)
                coords_geom[..., 0] = (coords_geom[..., 0] - lon_shift + 180) % 360 - 180
                _update_extent(coords_geom)
                plot_element = ax.plot(
                    coords_geom[..., 0],
                    coords_geom[..., 1],
                    color=palette[label],
                    linewidth=2,
                    label=label,
                )[0]
            elif isinstance(geom, Polygon):
                coords_geom = np.array(geom.exterior.coords)
                coords_geom[..., 0] = (coords_geom[..., 0] - lon_shift + 180) % 360 - 180
                _update_extent(coords_geom)
                plot_element = ax.fill(
                    coords_geom[..., 0],
                    coords_geom[..., 1],
                    edgecolor=palette[label],
                    facecolor=palette[label] + (0.2,),
                    linewidth=1,
                    label=label,
                )[0]
            else:
                raise TypeError(
                    f"Did not know how to handle geometry type: {type(geoms[0])}"
                )
        assert plot_element
        return plot_element

    # loop, plotting matchups / matches, and track lines/polygons for the legend
    for matchup in matchups:
        bounds_pace = matchup.get_pace_bounds()
        plot_elements[matchup.shortname_pace] = _plot_bounds(
            bounds_pace, matchup.shortname_pace
        )
        for match in matchup.matches_earthcare:
            bounds_earthcare = match.get_earthcare_bounds()
            label_earthcare = (
                f"EarthCARE_{match.filepath_earthcare.parts[-2].split('_')[0]}"
            )
            plot_elements[label_earthcare] = _plot_bounds(
                bounds_earthcare, label_earthcare
            )
    ax.set_title("PACE / EarthCARE Matchups")
    ax.legend(handles=plot_elements.values())
    ax.set_xlim(max(-180, minlon - 5), min(180, maxlon + 5))
    
    ax.set_ylim(max(-90, minlat - 5), min(90, maxlat + 5))

    if fig_filepath:
        plt.savefig(fig_filepath, dpi=200, bbox_inches="tight")
