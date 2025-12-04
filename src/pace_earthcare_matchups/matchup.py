"""TODO"""

from collections import defaultdict
from dataclasses import dataclass
from datetime import datetime, timedelta
import os
from pathlib import Path
import warnings

import h5py
from maap.maap import MAAP
from maap.Result import Granule
import netCDF4
import numpy as np
import numpy.typing as npt
from pystac_client import Client
from pystac.item import Item
from shapely import (
    Point,
    Polygon,
    to_geojson,
)
from tqdm import tqdm

from pace_earthcare_matchups.earthcare import (
    download_earthcare_item,
    parse_earthcare_filename,
)
from pace_earthcare_matchups.geospatial_utils import (
    get_centering_function,
)
from pace_earthcare_matchups.metadata_utils import (
    UTC,
    filter_granules,
    get_datetime_range_from_granule,
    get_datetime_range_maap,
    get_intersection_bbox,
    polygon_from_granule,
)
from pace_earthcare_matchups.path_utils import PATH_DATA, get_path
from pace_earthcare_matchups.supported_products import (
    EARTHCARE_SHORTNAMES,
    PACE_SHORTNAMES,
)


CMR_HOST = "cmr.earthdata.nasa.gov"


@dataclass
class MetaMatchEarthcare:
    """TODO"""

    item: Item
    bbox: Polygon


@dataclass
class MetaMatchup:
    """TODO"""

    granule_pace: Granule
    matches_earthcare: list[MetaMatchEarthcare]


@dataclass
class MatchEarthcare:
    filepath_earthcare: Path
    mask: npt.NDArray[np.bool]


@dataclass
class Matchup:
    filepath_pace: Path
    matches_earthcare: list[MatchEarthcare]

    def save(self) -> None:
        data_pace = netCDF4.Dataset(self.filepath_pace)
        pace_type = f"{data_pace.instrument}_{data_pace.processing_level}"
        stem_pace = self.filepath_pace.stem
        pace_path = PATH_DATA / "matchups" / pace_type / stem_pace
        for match in self.matches_earthcare:
            namedata = parse_earthcare_filename(match.filepath_earthcare)
            mdir = pace_path / namedata.get_file_type()
            mpath = mdir / match.filepath_earthcare.stem
            os.makedirs(mpath.parent, exist_ok=True)
            np.save(mpath, match.mask, allow_pickle=False)


def get_meta_matchup_from_granule(
    maap: MAAP,
    client_esa: Client,
    granule_pace: Granule,
    shortnames_earthcare: list[str],
    time_offset: timedelta = timedelta(seconds=5),
) -> MetaMatchup | None:
    """TODO"""
    assert isinstance(shortnames_earthcare, list)

    try:
        pace_poly = polygon_from_granule(granule_pace)
    except ValueError:
        warnings.warn(
            f"Broken geospatial bounds in {granule_pace['Granule']['DataGranule']['ProducerGranuleId']}"
        )
        return
    dt_range_offset = get_datetime_range_from_granule(granule_pace, time_offset)
    # check for granules with broken timestamps
    if dt_range_offset[1] - dt_range_offset[0] > timedelta(days=1):
        return None
    items_ec = []
    for shortname_ec in shortnames_earthcare:
        results_esa = client_esa.search(
            collections=["EarthCAREL1Validated_MAAP", "EarthCAREL2Validated_MAAP"],
            datetime=dt_range_offset,
            intersects=to_geojson(pace_poly),
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
) -> npt.NDArray[np.bool]:
    """TODO: document and break into smaller functions
    TODO: handle 2D EarthCARE products with a mask
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

    # get mask of where the EarthCARE track is in the PACE granule
    ec_in_granule = np.zeros_like(lat_ec, dtype=bool)
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
    # # get the longest contiguous
    # ec_in_granule_diff = np.diff(ec_in_granule.astype(int), prepend=np.zeros(1), append=np.zeros(1))
    # ec_in_granule_start_idxs = np.where(ec_in_granule_diff == 1)[0]
    # ec_in_granule_end_idxs = np.where(ec_in_granule_diff == -1)[0]
    # ec_contig_lens = ec_in_granule_end_idxs - ec_in_granule_start_idxs
    # longest_contig_idx = np.argmax(ec_contig_lens)
    # ec_in_granule_start_idx = ec_in_granule_start_idxs[longest_contig_idx]
    # ec_in_granule_end_idx = ec_in_granule_end_idxs[longest_contig_idx]
    # return ec_in_granule_start_idx, ec_in_granule_end_idx


def get_matchup(
    meta_matchup: MetaMatchup,
    long_term_token: str,
):
    """TODO"""
    paths_earthcare = []
    for meta_match in meta_matchup.matches_earthcare:
        path_earthcare = get_path(meta_match.item)
        if not path_earthcare.exists():
            download_earthcare_item(
                item=meta_match.item,
                long_term_token=long_term_token,
                datadir=get_path(meta_match.item).parent,
            )
        paths_earthcare.append(path_earthcare)
    path_pace = get_path(meta_matchup.granule_pace)
    if not path_pace.exists():
        os.makedirs(path_pace.parent, exist_ok=True)
        meta_matchup.granule_pace.getData(str(path_pace.parent))
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
        # start_idx, end_idx = get_matchup_start_end_idx(lat_pace, lon_pace, lat_earthcare, lon_earthcare)
        match_mask = get_matchup_mask(
            lat_pace, lon_pace, lat_earthcare[()], lon_earthcare[()]
        )
        matches.append(
            MatchEarthcare(
                filepath_earthcare=path_earthcare,
                mask=match_mask,
            )
        )
    return Matchup(filepath_pace=path_pace, matches_earthcare=matches)


def get_matchups(
    maap: MAAP,
    client_esa: Client,
    shortname_pace: str,
    shortnames_earthcare: list[str],
    temporal: tuple[datetime, datetime],
    long_term_token: str,
    time_offset: timedelta = timedelta(),
    limit: int = 20,
    verbose: bool = True,
    save: bool = True,
) -> list[Matchup]:
    """TODO"""
    assert shortname_pace in PACE_SHORTNAMES
    assert isinstance(shortnames_earthcare, list)
    if not temporal[0].tzinfo:
        temporal = (temporal[0].astimezone(UTC), temporal[1].astimezone(UTC))
    for shortname in shortnames_earthcare:
        assert shortname in EARTHCARE_SHORTNAMES
    results_pace = filter_granules(
        maap.searchGranule(
            cmr_host=CMR_HOST,
            short_name=shortname_pace,
            temporal=get_datetime_range_maap(*temporal),
            limit=20,
        )
    )
    time_start = temporal[0]
    matches = []
    while len(results_pace) > 0:
        if verbose:
            print(
                f"{len(matches)} matches so far, searching {len(results_pace)} {shortname_pace} results"
            )
        pbar = tqdm(results_pace) if verbose else results_pace
        for result_pace in pbar:
            dt_range = get_datetime_range_from_granule(result_pace)
            assert dt_range[1] > dt_range[0]
            meta_match = get_meta_matchup_from_granule(
                maap=maap,
                client_esa=client_esa,
                granule_pace=result_pace,
                shortnames_earthcare=shortnames_earthcare,
                time_offset=time_offset,
            )
            if meta_match:
                match = get_matchup(meta_match, long_term_token)
                if save:
                    match.save()
                matches.append(match)
                if len(matches) >= limit:
                    break
            time_start = max(time_start, dt_range[1] + timedelta(seconds=1))
        if len(matches) >= limit:
            break
        if time_start >= temporal[1]:
            break
        results_pace = filter_granules(
            maap.searchGranule(
                cmr_host=CMR_HOST,
                short_name=shortname_pace,
                temporal=get_datetime_range_maap(time_start, temporal[1]),
                limit=20,
            )
        )
    if verbose:
        print(f"Found {len(matches)} total matches!")
    return matches
