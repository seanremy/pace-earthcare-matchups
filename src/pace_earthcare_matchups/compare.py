"""This module handles intercomparison using already downloaded matchups.

Thanks to Andrzej Wasilewski and Snorre Stamnes for their help with the
intercomparison code!
"""

from collections import defaultdict

import netCDF4
import numpy as np
import numpy.typing as npt
from scipy.spatial import KDTree

from pace_earthcare_matchups.earthcare import get_earthcare_latlon
from pace_earthcare_matchups.geospatial_utils import geo2ecef
from pace_earthcare_matchups.matchup import Matchup
from pace_earthcare_matchups.pace import get_pace_latlon


def interp_inverse_distance_weighting(
    img_src: npt.NDArray,
    lat_src: npt.NDArray,
    lon_src: npt.NDArray,
    alt_src: npt.NDArray | None,
    lat_dest: npt.NDArray,
    lon_dest: npt.NDArray,
    alt_dest: npt.NDArray | None,
    k: int = 4,
    remove_oob: bool = True,
) -> npt.NDArray:
    """Interpolate a geolocated image to new coords with inverse distance weighting.

    To speed up the k-nearest-neighbor search, the source coordinates are indexed with a
    kd-tree. By default, this function sets out-of-bounds locations to NaN. A location
    is considered out-of-bounds if any of the k-nearest neighbors include the image
    boundary. Higher k values will lead to more of the image being "cropped".

    :param img_src: The source image, shape (U1, V1) or (U1, V1, N).
    :param lat_src: The source latitudes in degrees, shape (U1, V1).
    :param lon_src: The source longitudes in degrees, shape (U1, V1).
    :param alt_src: The source altitudes in meters above the reference ellipsoid
        (WGS84), shape (U1, V1), or None to assume sea level altitude.
    :param lat_dest: The destination latitudes in degrees, shape (U2, V2).
    :param lon_dest: The destination longitudes in degrees, shape (U2, V2).
    :param alt_dest The destination altitudes in meters above the reference ellipsoid
        (WGS84), shape (U2, V2), or None to assume sea level altitude.
    :param k: The number of nearest neighbors to sample for inverse distance weighting.
    :param remove_oob: A flag that controls whether out-of-bounds locations will be set
        to NaN, default=True.

    Returns:
        img_dest: The interpolated image, shape (U2, V2) or (U2, V2, N).
    """
    if alt_src is None or alt_dest is None:
        assert alt_src == alt_dest
        alt_src = np.zeros_like(lat_src)
        alt_dest = np.zeros_like(lat_dest)
    assert len(img_src.shape) in [2, 3]
    assert len(lat_src.shape) == 2
    assert len(lat_dest.shape) == 2
    assert img_src.shape[:2] == lat_src.shape == lon_src.shape == alt_src.shape
    assert lat_dest.shape == lon_dest.shape == alt_dest.shape
    shp_src = img_src.shape
    shp_dest = lat_dest.shape
    if len(img_src.shape) == 2:
        img_src = img_src[..., None]
    else:
        shp_dest += (shp_src[-1],)
    # project latitude and longitude to ECEF for better distances
    xyz_src = np.stack(
        geo2ecef(lat_src.ravel(), lon_src.ravel(), alt_src.ravel()), axis=-1
    )
    xyz_dest = np.stack(
        geo2ecef(lat_dest.ravel(), lon_dest.ravel(), alt_dest.ravel()), axis=-1
    )
    # kdtree for the source locations, query at the destination
    non_nan = ~np.isnan(xyz_src).any(axis=1)
    kdtree = KDTree(xyz_src[non_nan])
    dists, match_idx = kdtree.query(xyz_dest.reshape((-1, 3)), k=k)
    match_idx = np.where(non_nan)[0][match_idx]
    # inverse distance weighting
    wts_interp = np.zeros_like(dists)
    wts_interp[dists > 0] = 1 / dists[dists > 0]
    wts_interp_sum = wts_interp.sum(axis=1)
    np.divide(
        wts_interp,
        wts_interp_sum[:, None],
        where=wts_interp_sum[:, None] > 0,
        out=wts_interp,
    )
    vals_interp = img_src.reshape((-1,) + img_src.shape[2:])[match_idx]
    img_dest = (wts_interp[..., None] * vals_interp).sum(axis=1)
    # determine "out-of-bounds" (computing the full convex hull is too slow)
    # unfortunately, this includes any point whose k-nearest neighbors include an edge
    if remove_oob:
        midx_min = match_idx.min(axis=1)
        midx_max = match_idx.max(axis=1)
        oob_top = midx_min // shp_src[1] == 0
        oob_left = midx_min % shp_src[1] == 0
        oob_bot = midx_max // shp_src[1] == shp_src[0] - 1
        oob_right = midx_max % shp_src[1] == shp_src[1] - 1
        oob = oob_top + oob_left + oob_bot + oob_right
        img_dest[oob] = np.nan
    return img_dest.reshape(shp_dest)


def get_comparison_dict(
    matchup: Matchup,
    pace_variables: list[str],
) -> defaultdict[str, defaultdict[int, dict[str, npt.NDArray]]]:
    """Build a nested comparison dictionary by interpolating PACE variables onto
    EarthCARE geolocations.

    For each EarthCARE product in the matchup, and for each individual EarthCARE
    granule matched to the PACE file, PACE variables are interpolated to the EarthCARE
    lat/lon positions using inverse distance weighting. The interpolated output matches
    the dimensionality of the EarthCARE product; if matching with ATLID or CPR, the
    output is a curtain, with MSI, it is an image, and with BBR, it is a

    :param matchup: A :class:`~pace_earthcare_matchups.matchup.Matchup` containing the
        PACE filepath, EarthCARE filepaths, and per-granule overlap masks.
    :param pace_variables: List of NetCDF variable paths to read from the PACE file and
        interpolate onto EarthCARE coordinates (e.g. `["observation_data/i", ...]`).
    :returns: A nested dict with structure
        ``{earthcare_product: {match_index: {key: array}}}``, where each inner dict
            contains:

        - ``"latitude"`` -- EarthCARE latitudes at matched locations, shape ``(M,)``
            or ``(M, 1)``.
        - ``"longitude"`` -- EarthCARE longitudes at matched locations, shape ``(M,)``
            or ``(M, 1)``.
        - ``"start"`` / ``"end"`` -- slice indices into the full EarthCARE array
          (only present for 2-D masks).
        - one entry per name in *pace_variables* -- the PACE variable interpolated to
          EarthCARE locations, same leading shape as the latitude/longitude arrays.
    """

    matches = defaultdict(list)
    for match in matchup.matches_earthcare:
        ec_short_name = match.filepath_earthcare.parent.name
        matches[ec_short_name].append(match)

    pace_arrs = []
    with netCDF4.Dataset(matchup.filepath_pace) as data_pace:
        for pace_var in pace_variables:
            pace_arrs.append(data_pace[pace_var][()].filled(fill_value=np.nan))
        pace_arr = np.stack(pace_arrs, axis=-1)
        # TODO: handle multi-angle cases
        lat_pace, lon_pace = get_pace_latlon(matchup.filepath_pace)

    comp_dict = defaultdict(lambda: defaultdict(lambda: dict()))
    for ec_short_name, match_list in matches.items():
        for i, match in enumerate(match_list):
            # TODO: handle BBR's 3 angles
            lat_ec, lon_ec = get_earthcare_latlon(match.filepath_earthcare)
            if len(match.mask.shape) == 2:
                mask_idx = np.where(match.mask)
                start = mask_idx[0].min()
                end = mask_idx[0].max() + 1
                lat_ec = lat_ec[start:end]
                lon_ec = lon_ec[start:end]
                comp_dict[ec_short_name][i]["start"] = start
                comp_dict[ec_short_name][i]["end"] = end
            else:
                lat_ec = lat_ec[match.mask]
                lon_ec = lon_ec[match.mask]
                lat_ec, lon_ec = lat_ec[:, None], lon_ec[:, None]
            comp_dict[ec_short_name][i]["latitude"] = lat_ec
            comp_dict[ec_short_name][i]["longitude"] = lon_ec
            # TODO: support more than IDW
            pace_arr_interp = interp_inverse_distance_weighting(
                pace_arr,
                lat_pace,
                lon_pace,
                None,
                lat_ec,
                lon_ec,
                None,
            )
            for j, pace_var in enumerate(pace_variables):
                comp_dict[ec_short_name][i][pace_var] = pace_arr_interp[..., j]
    return comp_dict
