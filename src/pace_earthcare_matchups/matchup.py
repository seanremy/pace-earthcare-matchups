"""This module handles metadata matchups and matchups between PACE and EarthCARE. Note
the difference between a "metadata matchup" and a "matchup". Metadata matchups use only
the metadata descriptors stored in the metadata repositories (NASA CMR and ESA STAC).
Matchups use the data arrays in the files themselves, which is more accurate, but
requires the files to be downloaded locally. A two-stage approach of performing metadata
comparisons first vastly reduces the number of files which must be locally downloaded.
"""

from collections import defaultdict
from dataclasses import dataclass
from datetime import datetime, timedelta
import os
from pathlib import Path
import re
import shutil
from typing import Callable
from zoneinfo import ZoneInfo

import h5py
import netCDF4
import numpy as np
import numpy.typing as npt
from pystac_client import Client
from pystac.item import Item
from shapely import (
    LineString,
    MultiLineString,
    MultiPolygon,
    Point,
    Polygon,
    to_geojson,
)
from tqdm.notebook import tqdm

from pace_earthcare_matchups.earthcare import (
    download_earthcare_item,
    download_missing_earthcare_data,
    parse_earthcare_filename,
)
from pace_earthcare_matchups.geospatial_utils import (
    correct_linestring,
    correct_polygon,
    get_centering_function,
    get_outer_ring,
)
from pace_earthcare_matchups.metadata_utils import (
    get_datetime_range_from_granule,
    get_intersection_bbox,
)
from pace_earthcare_matchups.pace import (
    Granule,
    _query_cmr,
    download_missing_pace_data,
    get_nadir_idx_harp2_l1b,
    get_pace_shortname,
)
from pace_earthcare_matchups.path_utils import PATH_DATA, PATH_TOKEN, get_path
from pace_earthcare_matchups.supported_products import (
    EARTHCARE_SHORTNAMES,
    PACE_SHORTNAMES,
)


@dataclass
class MetaMatchEarthcare:
    """Represents the metadata of an EarthCARE file matched to a PACE file's metadata.

    Args:
        item: STAC item describing an EarthCARE file.
        bbox: Bounding box where the EarthCARE bounding geometry intersects the parent
            PACE granule's bounding geometry.
    """

    item: Item
    bbox: Polygon


@dataclass
class MetaMatchup:
    """Represents a metadata match between a PACE and EarthCARE file.

    Args:
        granule_pace: A PACE granule's MAAP metadata.
        matches_earthcare: List of EarthCARE metadata matched to the PACE granule.
    """

    granule_pace: Granule
    matches_earthcare: list[MetaMatchEarthcare]


@dataclass
class MatchEarthcare:
    """Represents an EarthCARE file to which a PACE file was matched.

    Args:
        filepath_earthcare: Path of downloaded EarthCARE data.
        mask: Mask into the data where it overlaps geospatially with PACE data.
    """

    filepath_earthcare: Path
    mask: npt.NDArray

    def get_earthcare_bounds(
        self, pts_per_side: int = 20
    ) -> LineString | MultiLineString | Polygon | MultiPolygon:
        """Get the linestring or polygon bounds for this MatchEarthcare's data.

        Args:
            pts_per_side: Number of points per side in a polygon approximation of the
                latitude/longitude bounds of 2D EarthCARE data.
        """
        data_earthcare = h5py.File(self.filepath_earthcare)
        lat_data = data_earthcare["ScienceData/latitude"]
        lon_data = data_earthcare["ScienceData/longitude"]
        assert isinstance(lat_data, h5py.Dataset)
        assert isinstance(lon_data, h5py.Dataset)
        lat, lon = lat_data[()], lon_data[()]
        if len(lat.shape) == 1:
            return correct_linestring(LineString(np.stack([lon, lat], axis=-1)))
        elif len(lat.shape) == 2:
            # crop missing columns
            missing = ((lat > 1e35) + (lon > 1e35)).any(axis=0)
            lat, lon = lat[:, ~missing], lon[:, ~missing]
            idx_sides = [
                np.linspace(0, idx - 1, pts_per_side).astype(int) for idx in lat.shape
            ]
            idx0 = np.concatenate(
                [
                    idx_sides[0][:-1],
                    np.ones(pts_per_side - 1, dtype=int) * (lat.shape[0] - 1),
                    idx_sides[0][-1:0:-1],
                    np.zeros(pts_per_side - 1, dtype=int),
                ]
            )
            idx1 = np.concatenate(
                [
                    np.zeros(pts_per_side - 1, dtype=int),
                    idx_sides[1][:-1],
                    np.ones(pts_per_side - 1, dtype=int) * (lat.shape[1] - 1),
                    idx_sides[1][-1:0:-1],
                ]
            )
            return correct_polygon(
                Polygon(np.stack([lon[idx0, idx1], lat[idx0, idx1]], axis=-1))
            )
        else:
            raise ValueError(
                f"Expected lat/lon arrays to have 1 or 2 dimensions, got {len(lat.shape[0])}"
            )


@dataclass
class Matchup:
    """Represents a PACE file and a set of EarthCARE files which overlap with it.

    Args:
        filepath_pace: Path of downloaded PACE data.
        matches_earthcare: A list of MatchEarthcare overlapping the PACE file.
    """

    filepath_pace: Path
    shortname_pace: str
    matches_earthcare: list[MatchEarthcare]

    def save(self) -> None:
        """Save this Matchup to file."""
        data_pace = netCDF4.Dataset(self.filepath_pace)
        pace_type = f"{data_pace.instrument}_{data_pace.processing_level}"
        stem_pace = self.filepath_pace.stem
        path_pace = PATH_DATA / "matchups" / pace_type / stem_pace
        for match in self.matches_earthcare:
            namedata = parse_earthcare_filename(match.filepath_earthcare)
            mdir = path_pace / namedata.get_file_type()
            mpath = mdir / match.filepath_earthcare.stem
            os.makedirs(mpath.parent, exist_ok=True)
            np.save(mpath, match.mask, allow_pickle=False)

    def get_pace_bounds(self) -> Polygon | MultiPolygon:
        """Get the polygon bounds for this Matchup's PACE data."""
        data_pace = netCDF4.Dataset(self.filepath_pace)
        if hasattr(data_pace, "geospatial_bounds"):
            poly_str = (
                data_pace.geospatial_bounds.removeprefix("POLYGON")
                .lstrip(" ")
                .removeprefix("((")
                .removesuffix("))")
            )
            poly_arr = np.array(
                [pt.rstrip(" ").split(" ")[::-1] for pt in re.split(", *", poly_str)]
            ).astype(float)
            reorder_latlon_l2 = (
                data_pace.processing_version <= "3.1"
                and self.shortname_pace.rstrip("_NRT")
                in [
                    "PACE_OCI_L2_AOP",
                    "PACE_OCI_L2_BGC",
                    "PACE_OCI_L2_CLOUD",
                    "PACE_OCI_L2_IOP",
                    "PACE_OCI_L2_LANDVI",
                    "PACE_OCI_L2_PAR",
                    "PACE_OCI_L2_SFREFL",
                ]
            )
            if reorder_latlon_l2:
                poly_arr = poly_arr[..., ::-1]
            return correct_polygon(Polygon(poly_arr))
        elif self.shortname_pace == "PACE_HARP2_L1B_SCI":
            idx_nadir = get_nadir_idx_harp2_l1b(data_pace)
            lat = data_pace["geolocation_data/latitude"][idx_nadir].filled(
                fill_value=np.nan
            )
            lon = data_pace["geolocation_data/longitude"][idx_nadir].filled(
                fill_value=np.nan
            )
            latlon_ring = get_outer_ring(np.stack([lon, lat], axis=-1))
            return correct_polygon(Polygon(latlon_ring))
        else:
            raise NotImplementedError


def get_meta_matchup_from_granule(
    client_esa: Client,
    granule_pace: Granule,
    shortnames_earthcare: list[str],
    time_offset: timedelta = timedelta(seconds=5),
) -> MetaMatchup | None:
    """Get a metadata matchup from the provided PACE granule.

    Args:
        client_esa: pySTAC client to access ESA data.
        granule_pace: A PACE granule's metadata.
        shortnames_earthcare: List of EarthCARE collection short names.
        time_offset: This offset will be subtracted from the start time and added to the
            end time of the time range.

    Returns:
        meta_matchup: A metadata matchup with EarthCARE metadata overlapping the
            provided PACE granule.
    """
    assert isinstance(shortnames_earthcare, list)

    dt_range_offset = get_datetime_range_from_granule(granule_pace, time_offset)
    # check for granules with broken timestamps
    if dt_range_offset[1] - dt_range_offset[0] > timedelta(days=1):
        return None
    items_ec = []
    for shortname_ec in shortnames_earthcare:
        results_esa = client_esa.search(
            collections=["EarthCAREL1Validated_MAAP", "EarthCAREL2Validated_MAAP"],
            datetime=dt_range_offset,
            intersects=to_geojson(granule_pace.geospatial_bounds),
            method="GET",
            filter=f"productType = '{shortname_ec}'",
        )
        items_ec += list(results_esa.items())
    if len(items_ec) > 0:
        # arrange by file type and start time, keep the newest baseline
        items_dict = defaultdict(list)
        for item in items_ec:
            items_dict[item.id[9:36]].append(item)
        items_ec = [
            sorted(items, key=lambda i: i.id[6:8])[-1] for items in items_dict.values()
        ]
        meta_matches = []
        for item in items_ec:
            inter_bbox = get_intersection_bbox(granule_pace, item)
            meta_matches.append(
                MetaMatchEarthcare(
                    item=item,
                    bbox=inter_bbox,
                )
            )
        return MetaMatchup(granule_pace=granule_pace, matches_earthcare=meta_matches)


def get_matchup_mask(
    lat_pace: npt.NDArray[np.float32],
    lon_pace: npt.NDArray[np.float32],
    lat_ec: npt.NDArray[np.float32],
    lon_ec: npt.NDArray[np.float32],
) -> npt.NDArray:
    """Get the mask of the overlap between PACE and EarthCARE lat/lon arrays. The mask
    is in the shape of the provided PACE data. PACE data is expected to be a 2D array.

    TODO: account for multi-angle L1B lat/lon arrays.
    TODO: break into smaller functions.

    Args:
        lat_pace: PACE latitude array.
        lon_pace: PACE longitude array.
        lat_ec: EarthCARE latitude array.
        lon_ec: EarthCARE longitude array.

    Returns:
        ec_in_granule: Mask into the PACE data where it overlaps with the provided
            EarthCARE data.
    """
    num_pts_per_edge = 10
    # TODO: account for dimensions of different products
    vspace = (lat_pace.shape[0] - 1) / (num_pts_per_edge - 1)
    hspace = (lat_pace.shape[1] - 1) / (num_pts_per_edge - 1)
    pts_left = []
    pts_bot = []
    pts_right = []
    pts_top = []
    for i in range(1, num_pts_per_edge - 1):
        pts_left.append((int(i * vspace), 0))
        pts_bot.append((lat_pace.shape[0] - 1, int(i * hspace)))
        pts_right.append(
            (lat_pace.shape[0] - 1 - int(i * vspace), lat_pace.shape[1] - 1)
        )
        pts_top.append((0, lat_pace.shape[1] - 1 - int(i * hspace)))
    poly_idx = [(0, 0)]
    poly_idx += pts_left
    poly_idx += [(lat_pace.shape[0] - 1, 0)]
    poly_idx += pts_bot
    poly_idx += [(lat_pace.shape[0] - 1, lat_pace.shape[1] - 1)]
    poly_idx += pts_right
    poly_idx += [(0, lat_pace.shape[1] - 1)]
    poly_idx += pts_top
    poly_idx += [(0, 0)]

    centering_fn = get_centering_function(lat_pace, lon_pace)
    lat_rot_ec, lon_rot_ec = centering_fn(lat_ec, lon_ec)

    # get and rotate the lat/lon at these polygon locations, then take the mean for each angle
    lat_poly = lat_pace[[p[0] for p in poly_idx], [p[1] for p in poly_idx]]
    lon_poly = lon_pace[[p[0] for p in poly_idx], [p[1] for p in poly_idx]]
    lat_rot_poly, lon_rot_poly = centering_fn(lat_poly, lon_poly)

    # convert to shapely polygon
    pace_poly = Polygon(zip(lon_rot_poly, lat_rot_poly))

    # filter out with a bbox first for speed
    lat_min_poly = np.nanmin(lat_rot_poly)
    lon_min_poly = np.nanmin(lon_rot_poly)
    lat_max_poly = np.nanmax(lat_rot_poly)
    lon_max_poly = np.nanmax(lon_rot_poly)
    latlon_rot_bbox_mask = (
        (lat_rot_ec >= lat_min_poly)
        * (lat_rot_ec <= lat_max_poly)
        * (lon_rot_ec >= lon_min_poly)
        * (lon_rot_ec <= lon_max_poly)
    )

    # get mask of where the EarthCARE points are in the PACE granule
    ec_in_granule = np.zeros_like(lat_ec, dtype=bool)
    if latlon_rot_bbox_mask.any():
        pace_contains = np.vectorize(
            lambda p: pace_poly.contains(Point(p)), signature="(n)->()"
        )
        ec_in_granule[latlon_rot_bbox_mask] = pace_contains(
            np.stack(
                [
                    lon_rot_ec[latlon_rot_bbox_mask],
                    lat_rot_ec[latlon_rot_bbox_mask],
                ],
                axis=-1,
            )
        )
    return ec_in_granule


def get_matchup(
    meta_matchup: MetaMatchup,
) -> tuple[Matchup, list[Path]]:
    """Get a matchup from a metadata matchup. Downloads the associated EarthCARE and
    PACE data, then get the mask describing their overlap.

    Args:
        meta_matchup: A metadata matchup.

    Returns:
        matchup: The matchup derived from the provided metadata matchup.
        TODO
    """
    paths_added = []
    paths_earthcare = []
    for meta_match in meta_matchup.matches_earthcare:
        path_earthcare = get_path(meta_match.item)
        if not path_earthcare.exists():
            download_earthcare_item(
                item=meta_match.item,
                datadir=path_earthcare.parent,
            )
            paths_added.append(path_earthcare)
        paths_earthcare.append(path_earthcare)
    path_pace = get_path(meta_matchup.granule_pace)
    if not path_pace.exists():
        os.makedirs(path_pace.parent, exist_ok=True)
        meta_matchup.granule_pace.download()
        paths_added.append(path_pace)
    data_pace = netCDF4.Dataset(path_pace)
    if "geolocation_data" in data_pace.groups:
        lat_pace = data_pace["geolocation_data/latitude"][:].filled(fill_value=np.nan)
        lon_pace = data_pace["geolocation_data/longitude"][:].filled(fill_value=np.nan)
    elif "navigation_data" in data_pace.groups:
        lat_pace = data_pace["navigation_data/latitude"][:].filled(fill_value=np.nan)
        lon_pace = data_pace["navigation_data/longitude"][:].filled(fill_value=np.nan)
    else:
        print(data_pace)
        raise NotImplementedError("Don't know how to parse this PACE product!")

    matches = []
    for path_earthcare in paths_earthcare:
        data_earthcare = h5py.File(path_earthcare)
        lat_earthcare = data_earthcare["ScienceData/latitude"]
        lon_earthcare = data_earthcare["ScienceData/longitude"]
        assert isinstance(lat_earthcare, h5py.Dataset)
        assert isinstance(lon_earthcare, h5py.Dataset)
        if meta_matchup.granule_pace.short_name == "PACE_HARP2_L1B_SCI":
            # in the case of HARP2 L1B, use the nadir view angle's geolocation as bounds
            idx_nadir = get_nadir_idx_harp2_l1b(data_pace)
            match_mask = get_matchup_mask(
                lat_pace[idx_nadir],
                lon_pace[idx_nadir],
                lat_earthcare[()],
                lon_earthcare[()],
            )
        else:
            match_mask = get_matchup_mask(
                lat_pace,
                lon_pace,
                lat_earthcare[()],
                lon_earthcare[()],
            )
        matches.append(
            MatchEarthcare(
                filepath_earthcare=path_earthcare,
                mask=match_mask,
            )
        )
    return Matchup(
        filepath_pace=path_pace,
        shortname_pace=meta_matchup.granule_pace.short_name,
        matches_earthcare=matches,
    ), paths_added


def get_matchups(
    shortname_pace: str,
    shortnames_earthcare: list[str],
    temporal: tuple[datetime, datetime],
    time_offset: timedelta = timedelta(),
    bbox: tuple[float, float, float, float] = (-180, -90, 180, 90),
    limit: int = 10,
    search_batch_size: int = 20,
    verbose: bool = True,
    save: bool = True,
    filter_fn: Callable | None = None,
) -> list[Matchup]:
    """Get a list of matchups using provided search arguments.

    This function searches for PACE data matching the provided filters in batches,
    incrementing the start time of the search window so as to search through the entire
    range in batches. For each item in a batch of PACE results, this function finds
    and downloads matching EarthCARE data, computes the mask of overlaps, and optionally
    saves the matchup data to disk.

    Args:
        shortname_pace: PACE collection short name.
        shortnames_earthcare: EarthCARE collection short names.
        temporal: The time range in which to retrieve data. Times are assumed to be UTC.
        time_offset: This offset will be subtracted from the start time and added to the
            end time of the time range.
        bbox: Lat/lon bounding box in W, S, E, N order by which to limit the search.
        limit: Limit on how many matchups to download.
        search_batch_size: How many PACE files to get per batch.
        verbose: If true, print update messages.
        save: Whether to save matchups to disk.
        filter_fn: A callable function which operates on a single matchup, returning
            True if the matchup is to be kept, False otherwise.

    Returns:
        matchups: List of retrieved matchups.
    """
    assert shortname_pace in PACE_SHORTNAMES
    assert isinstance(shortnames_earthcare, list)
    for shortname in shortnames_earthcare:
        assert shortname in EARTHCARE_SHORTNAMES

    temporal = (
        temporal[0].replace(tzinfo=ZoneInfo("UTC")),
        temporal[1].replace(tzinfo=ZoneInfo("UTC")),
    )
    results_pace = _query_cmr(
        short_name=shortname_pace,
        temporal=temporal,
        bbox=bbox,
        limit=search_batch_size,
    )
    client_esa = Client.open("https://catalog.maap.eo.esa.int/catalogue/")

    time_start = temporal[0]
    matches = []
    while len(results_pace) > 0:
        if verbose:
            t_start_str = datetime.strftime(results_pace[0].beginning_datetime, "%Y-%m-%d|%H:%M:%S")
            t_end_str = datetime.strftime(results_pace[-1].ending_datetime, "%Y-%m-%d|%H:%M:%S")
            print(f"{len(matches)}/{limit} matches so far.")
            print(f"Batch of {len(results_pace)} new {shortname_pace} results found "
                  f"from {t_start_str} -> {t_end_str}.")
            pbar = tqdm(results_pace, desc=f"Searching batch")
        else:
            pbar = results_pace
        for result_pace in pbar:
            dt_range = get_datetime_range_from_granule(result_pace)
            assert dt_range[1] > dt_range[0]
            meta_match = get_meta_matchup_from_granule(
                client_esa=client_esa,
                granule_pace=result_pace,
                shortnames_earthcare=shortnames_earthcare,
                time_offset=time_offset,
            )
            if meta_match:
                match, paths_added = get_matchup(meta_match)
                if filter_fn is None or filter_fn(match):
                    if save:
                        match.save()
                    matches.append(match)
                    if len(matches) >= limit:
                        break
                else:
                    for path in paths_added:
                        os.remove(path)
            time_start = max(time_start, dt_range[1] + timedelta(seconds=1))
        if len(matches) >= limit:
            break
        if time_start >= temporal[1]:
            break
        results_pace = _query_cmr(
            short_name=shortname_pace,
            temporal=(time_start, temporal[1]),
            bbox=bbox,
            limit=search_batch_size,
        )
    if verbose:
        print(f"Found {len(matches)} total matches!")
    return matches


def load_matchup(
    filepath: Path,
    download_missing: bool = False,
    client_esa: Client | None = None,
    shortnames_earthcare: list[str] | None = None,
) -> Matchup:
    """TODO
    :param filepath: Standardized path to the PACE name element of a matchup's mask filepath.
    :type filepath: Path
    :param download_missing:
    """
    assert filepath.stem.startswith("PACE_")

    parts = filepath.parts
    instrument_pace, level_pace = parts[-3], parts[-2]
    filestem_pace = parts[-1]
    shortname_pace = get_pace_shortname(instrument_pace, level_pace)
    filepath_pace = PATH_DATA / "PACE" / instrument_pace / level_pace / f"{filestem_pace}.nc"
    if download_missing and not filepath_pace.exists():
        download_missing_pace_data(filepath_pace)

    products_ec = sorted([fp.name for fp in filepath.glob("*/")])
    if shortnames_earthcare:
        products_ec = [p for p in products_ec if p in shortnames_earthcare]

    long_term_token = open(PATH_TOKEN).read()
    matches_earthcare = []
    for prod in products_ec:
        filepaths_mask = sorted((filepath / prod).glob("*"))
        for fp in filepaths_mask:
            filepath_earthcare = PATH_DATA / "EarthCARE" / prod / f"{fp.stem}.h5"
            if download_missing and not filepath_earthcare.exists():
                if not isinstance(long_term_token, str):
                    raise ValueError("You must provide a long_term_token to download missing EarthCARE data!")
                if not isinstance(client_esa, Client):
                    raise ValueError("You must provide a pystac Client to download missing EarthCARE data!")
                download_missing_earthcare_data(filepath_earthcare, long_term_token, client_esa)
            matches_earthcare.append(MatchEarthcare(
                filepath_earthcare=filepath_earthcare,
                mask=np.load(fp, allow_pickle=False),
            ))
    return Matchup(
        filepath_pace=filepath_pace,
        shortname_pace=shortname_pace,
        matches_earthcare=matches_earthcare,
    )


def delete_matchup(
    matchup: Matchup,
    delete_associated_files: bool = False,
) -> None:
    """TODO
    
    Warning: If you have any other data stored in the matchup path, it will be
    deleted too! For this reason it is highly suggested not to store data in
    the matchups data directory.

    TODO
    """
    matchup_path = PATH_DATA / "matchups" / matchup.filepath_pace.parts[-3] / matchup.filepath_pace.parts[-2] / matchup.filepath_pace.stem
    shutil.rmtree(matchup_path)
    if not delete_associated_files:
        return
    os.remove(matchup.filepath_pace)
    for match in matchup.matches_earthcare:
        os.remove(match.filepath_earthcare)


def get_all_matchup_paths() -> list[Path]:
    """TODO"""
    return sorted((PATH_DATA / "matchups").glob("*/*/"))
