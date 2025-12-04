"""TODO"""

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

from pace_earthcare_matchups.matchup import Matchup


def plot_matchups(
    matchups: list[Matchup],
    figsize: tuple[int, int] | None = None,
    fig_filepath: Path | str | None = None,
) -> None:
    """TODO"""
    if isinstance(fig_filepath, str):
        fig_filepath = Path(fig_filepath)
    ax = plt.figure(figsize=figsize).add_subplot(projection=ccrs.PlateCarree())
    assert isinstance(ax, GeoAxes)
    ax.stock_img(alpha=0.8)
    ax.add_feature(cfeature.COASTLINE)

    plot_elements = {}
    # first, assign a color palette to each instrument
    labels_pace = set()
    labels_earthcare = set()
    for matchup in matchups:
        labels_pace.add("_".join(matchup.filepath_pace.parts[-4:-1]))
        for match in matchup.matches_earthcare:
            labels_earthcare.add(f"EarthCARE_{match.filepath_earthcare.parts[-2].split('_')[0]}")
    palette = dict(zip(
        list(labels_pace) + list(labels_earthcare),
        sns.color_palette("colorblind", len(labels_pace) + len(labels_earthcare)),
    ))
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
        bounds: LineString | MultiLineString | Polygon | MultiPolygon,
        label: str
    ) -> matplotlib.lines.Line2D | matplotlib.patches.Polygon:
        if isinstance(bounds, MultiLineString | MultiPolygon):
            geoms = bounds.geoms
        else:
            geoms = [bounds]
        plot_element = None
        for geom in geoms:
            if isinstance(geom, LineString):
                coords_geom = np.array(geom.coords)
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
                raise TypeError(f"Did not know how to handle geometry type: {type(geoms[0])}")
        assert plot_element
        return plot_element
    # loop, plotting matchups / matches, and track lines/polygons for the legend
    for matchup in matchups:
        bounds_pace = matchup.get_pace_bounds()
        label_pace = "_".join(matchup.filepath_pace.parts[-4:-1])
        plot_elements[label_pace] = _plot_bounds(bounds_pace, label_pace)
        for match in matchup.matches_earthcare:
            bounds_earthcare = match.get_earthcare_bounds()
            label_earthcare = f"EarthCARE_{match.filepath_earthcare.parts[-2].split('_')[0]}"
            plot_elements[label_earthcare] = _plot_bounds(bounds_earthcare, label_earthcare)
    ax.set_title(f"PACE / EarthCARE Matchups")
    ax.legend(handles=plot_elements.values())
    ax.set_xlim(max(-180, minlon - 5), min(180, maxlon + 5))
    ax.set_ylim(max(-90, minlat - 5), min(90, maxlat + 5))
    
    if fig_filepath:
        plt.savefig(fig_filepath, dpi=100, bbox_inches="tight")
