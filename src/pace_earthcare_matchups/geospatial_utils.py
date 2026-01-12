"""Set of utilities for geospatial operations. For a reference on coordinate frames, see
'PACE OCI Calibration and Geolocation Operational Algorithm Description,' Fred Patt
https://pace.oceansciences.org/docs/PACE_Vol_13_all_combined.pdf
"""

from collections.abc import Callable

import numpy as np
import numpy.typing as npt
from shapely import (
    LineString,
    MultiLineString,
    MultiPolygon,
    Polygon,
)


# WGS-84 constants
WGS_84_A = 6378137.0  # semimajor axis
WGS_84_B = 6356752.314245  # semiminor axis
WGS_84_E = (WGS_84_A**2 - WGS_84_B**2) / (WGS_84_A**2)  # first eccentricity
WGS_84_E2 = (WGS_84_A**2 - WGS_84_B**2) / (WGS_84_B**2)  # second eccentricity
WGS_84_F = (WGS_84_A - WGS_84_B) / WGS_84_A  # flattening


def vincenty_distance(
    latlon1: npt.NDArray[np.float32 | np.float64],
    latlon2: npt.NDArray[np.float32 | np.float64],
    tol: float = 1e-12,
    max_iters: int = 10,
) -> tuple[
    npt.NDArray[np.float32 | np.float64],
    npt.NDArray[np.float32 | np.float64],
    npt.NDArray[np.float32 | np.float64],
]:
    """Compute geodesic distance using Vincenty's formulae and the WGS-84 ellipsoid.

    See: https://en.wikipedia.org/wiki/Vincenty%27s_formulae#Inverse_problem
    Args:
        latlon1: Starting points' latitudes and longitudes as a numpy array of
            shape (2, N).
        latlon2: Destination points' latitudes and longitudes as a numpy array
            of shape (2, N).
        tol: Tolerance in meters. When the updates are less than tol, the
            iteration ends (default: 1e-12).
        max_iters: Maximum number of iterations to perform.

    Returns:
        s: Geodesic distance in meters between provided points.
        alpha1: Forward azimuths at starting points.
        alpha2: Forward azimuths at destination points.
    """
    assert isinstance(latlon1, np.ndarray) and latlon1.dtype == np.float64
    assert isinstance(latlon2, np.ndarray) and latlon2.dtype == np.float64
    assert isinstance(tol, float)
    assert isinstance(max_iters, int)
    assert latlon1.shape[0] == 2 and latlon2.shape[0] == 2
    assert latlon1.shape[1] == latlon2.shape[1]

    lat1, lat2 = latlon1[0] * np.pi / 180, latlon2[0] * np.pi / 180
    lon1, lon2 = latlon1[1] * np.pi / 180, latlon2[1] * np.pi / 180
    U1 = np.arctan((1 - WGS_84_F) * np.tan(lat1))
    U2 = np.arctan((1 - WGS_84_F) * np.tan(lat2))
    L = lon2 - lon1

    lambd = L
    lambd_diff = 1000
    num_iters = 0

    sin_sigma = np.zeros_like(U1)
    cos_sigma = np.zeros_like(U1)
    sigma = np.zeros_like(U1)
    cos2_alpha = np.zeros_like(U1)
    cos_2sigmam = np.zeros_like(U1)

    while (np.abs(lambd_diff) > tol).any():
        if num_iters > max_iters:
            raise Warning(
                f"Exceeded {max_iters} iterations without lambda changing by less than "
                f"{tol:.1e}"
            )

        sin_sigma = np.sqrt(
            (np.cos(U2) * np.sin(lambd)) ** 2
            + (np.cos(U1) * np.sin(U2) - np.sin(U1) * np.cos(U2) * np.cos(lambd)) ** 2
        )
        cos_sigma = np.sin(U1) * np.sin(U2) + np.cos(U1) * np.cos(U2) * np.cos(lambd)
        sigma = np.arctan2(sin_sigma, cos_sigma)

        sin_alpha = np.cos(U1) * np.cos(U2) * np.sin(lambd) / sin_sigma
        cos2_alpha = 1 - sin_alpha**2
        cos_2sigmam = cos_sigma - (2 * np.sin(U1) * np.sin(U2)) / cos2_alpha

        C = (WGS_84_F / 16) * cos2_alpha * (4 + WGS_84_F * (4 - 3 * cos2_alpha))
        lambd_i = L + (1 - C) * WGS_84_F * sin_alpha * (
            sigma
            + C * sin_sigma * (cos_2sigmam + C * cos_sigma * (-1 + 2 * cos_2sigmam**2))
        )
        lambd_diff = lambd_i - lambd
        lambd = lambd_i
        num_iters += 1

    u2 = cos2_alpha * (WGS_84_A**2 - WGS_84_B**2) / WGS_84_B**2
    A = 1 + (u2 / 16384) * (4096 + u2 * (-768 + u2 * (320 - 175 * u2)))
    B = (u2 / 1024) * (256 + u2 * (-128 + u2 * (74 - 47 * u2)))
    delta_sigma = (
        B
        * sin_sigma
        * (
            cos_2sigmam
            + (1 / 4)
            * B
            * (
                cos_sigma * (-1 + 2 * cos_2sigmam**2)
                - (1 / 6)
                * B
                * cos_2sigmam
                * (-3 + 4 * sin_sigma**2)
                * (-3 + 4 * cos_2sigmam**2)
            )
        )
    )

    s = WGS_84_B * A * (sigma - delta_sigma)
    alpha1 = np.arctan2(
        np.cos(U2) * np.sin(lambd),
        np.cos(U1) * np.sin(U2) - np.sin(U1) * np.cos(U2) * np.cos(lambd),
    )
    alpha2 = np.arctan2(
        np.cos(U1) * np.sin(lambd),
        -np.sin(U1) * np.cos(U2) + np.cos(U1) * np.sin(U2) * np.cos(lambd),
    )

    return s, alpha1 * 180 / np.pi, alpha2 * 180 / np.pi


def vincenty_point_along_geodesic(
    latlon1: npt.NDArray[np.float32 | np.float64],
    alpha1: npt.NDArray[np.float32 | np.float64],
    s: npt.NDArray[np.float32 | np.float64],
    tol: float = 1e-6,
    max_iters: int = 10,
) -> tuple[
    tuple[npt.NDArray[np.float32 | np.float64], npt.NDArray[np.float32 | np.float64]]
    | npt.NDArray[np.float32 | np.float64],
    npt.NDArray[np.float32 | np.float64],
]:
    """Compute destination locations along geodesics defined by starting locations,
    azimuth angles, and distances.

    See: https://en.wikipedia.org/wiki/Vincenty%27s_formulae#Direct_problem

    Args:
        latlon1: Starting points' latitudes and longitudes as a numpy array of shape
            (2, N).
        alpha1: Forward azimuths at initial points as a numpy array of shape (N,).
        s: Distances to travel along the geodesics as a numpy array of shape (N,).
        tol: Tolerance in meters. When the updates are less than tol, the
            iteration ends (default: 1e-12).
        max_iters: Maximum number of iterations to perform.

    Returns:
        latlon2: Estimated latitude and longitude as a numpy array of shape (2, N).
        alpha2: Forward azimuths at destination points as a numpy array of shape(N,).
    """
    assert isinstance(latlon1, np.ndarray) and latlon1.dtype == np.float64
    assert isinstance(alpha1, np.ndarray) and alpha1.dtype == np.float64
    assert isinstance(s, np.ndarray) and s.dtype == np.float64
    assert isinstance(tol, float)
    assert isinstance(max_iters, int) and max_iters > 0
    assert latlon1.shape[0] == 2
    assert latlon1.shape[1] == alpha1.shape[0] and latlon1.shape[1] == s.shape[0]

    lat1, lon1 = latlon1[0] * np.pi / 180, latlon1[1] * np.pi / 180
    alpha1 = alpha1 * np.pi / 180

    U1 = np.arctan((1 - WGS_84_F) * np.tan(lat1))
    sigma1 = np.arctan2(np.tan(U1), np.cos(alpha1))
    sin_alpha = np.cos(U1) * np.sin(alpha1)
    u2 = (1 - sin_alpha**2) * (WGS_84_A**2 - WGS_84_B**2) / WGS_84_B**2
    A = 1 + (u2 / 16384) * (4096 + u2 * (-768 + u2 * (320 - 175 * u2)))
    B = (u2 / 1024) * (256 + u2 * (-128 + u2 * (74 - 47 * u2)))

    sigma = s / (WGS_84_B * A)
    sigma_diff = 1000
    num_iters = 0

    cos_2sigmam = 0

    while (np.abs(sigma_diff) > tol).any():
        if num_iters > max_iters:
            raise Warning(
                f"Exceeded {max_iters} iterations without sigma changing by less than "
                f"{tol:.1e}"
            )

        cos_2sigmam = np.cos(2 * sigma1 + sigma)
        delta_sigma = (
            B
            * np.sin(sigma)
            * (
                cos_2sigmam
                + (1 / 4)
                * B
                * (
                    np.cos(sigma) * (-1 + 2 * cos_2sigmam**2)
                    - (1 / 6)
                    * B
                    * cos_2sigmam
                    * (-3 + 4 * np.sin(sigma) ** 2)
                    * (-3 + 4 * cos_2sigmam**2)
                )
            )
        )
        sigma_i = s / (WGS_84_B * A) + delta_sigma
        sigma_diff = sigma_i - sigma
        sigma = sigma_i
        num_iters += 1

    lat2 = np.arctan2(
        np.sin(U1) * np.cos(sigma) + np.cos(U1) * np.sin(sigma) * np.cos(alpha1),
        (1 - WGS_84_F)
        * np.sqrt(
            sin_alpha**2
            + (np.sin(U1) * np.sin(sigma) - np.cos(U1) * np.cos(sigma) * np.cos(alpha1))
            ** 2
        ),
    )
    lambd = np.arctan2(
        np.sin(sigma) * np.sin(alpha1),
        np.cos(U1) * np.cos(sigma) - np.sin(U1) * np.sin(sigma) * np.cos(alpha1),
    )
    C = (
        (WGS_84_F / 16)
        * (1 - sin_alpha**2)
        * (4 + WGS_84_F * (4 - 3 * (1 - sin_alpha**2)))
    )
    L = lambd - (1 - C) * WGS_84_F * sin_alpha * (
        sigma
        + C
        * np.sin(sigma)
        * (cos_2sigmam + C * np.cos(sigma) * (-1 + 2 * cos_2sigmam**2))
    )
    lon2 = L + lon1
    alpha2 = np.arctan2(
        sin_alpha,
        -np.sin(U1) * np.sin(sigma) + np.cos(U1) * np.cos(sigma) * np.cos(alpha1),
    )

    lat2, lon2 = lat2 * 180 / np.pi, lon2 * 180 / np.pi
    if isinstance(latlon1, tuple):
        latlon2 = (lat2, lon2)
    else:
        latlon2 = np.stack([lat2, lon2])
    return latlon2, alpha2


def get_antimeridian_intersection(
    latlon1: npt.NDArray[np.float32 | np.float64],
    latlon2: npt.NDArray[np.float32 | np.float64],
    tol: float = 1e-6,
    max_iters: int = 1000,
) -> npt.NDArray[np.float32 | np.float64]:
    """Get the intersection between the great circle defined by two endpoints and the
    antimeridian (dateline) using Vincenty's formulae.

    Args:
        latlon1: Starting points' latitudes and longitudes as a numpy array of
            shape (2, N).
        latlon2: Destination points' latitudes and longitudes as a numpy array
            of shape (2, N).
        tol: Tolerance in meters. When the updates are less than tol, the
            iteration ends (default: 1e-12).
        max_iters: Maximum number of iterations to perform.

    Returns:
        lat_am: Latitude of the great circle's intersection with the antimeridian.
    """
    dists, alpha1, _ = vincenty_distance(latlon1, latlon2)

    lon_sign = np.sign(latlon1[1])
    lon_dist = np.abs(180 - np.abs(latlon1[1])) + np.abs(180 - np.abs(latlon2[1]))

    # initial guess
    t = (180 - np.abs(latlon1[1])) / lon_dist
    latlon_am, _ = vincenty_point_along_geodesic(latlon1, alpha1, t * dists)
    am_lon_diffs = lon_sign * (180 - latlon_am[1] % 360)
    converged = np.abs(am_lon_diffs) < tol
    i = 0
    while i < max_iters and not converged.all():
        t = t * 0.5 + 0.5 * (t + am_lon_diffs / lon_dist)
        # assert (t >= 0).all() and (t <= 1).all()  # boundary condition
        latlon_am, _ = vincenty_point_along_geodesic(latlon1, alpha1, t * dists)
        am_lon_diffs = lon_sign * (180 - latlon_am[1] % 360)
        converged = np.abs(am_lon_diffs) < tol
        i += 1
    return latlon_am[0]


def correct_polygon(poly: Polygon) -> Polygon | MultiPolygon:
    """Correct a geospatial polygon. If it contains the North xor the South pole, add
    vertices to minimally affect the polygon's area while avoiding intersection with the
    antimeridian (dateline). Otherwise, if it crosses the antimeridian without including
    a pole, split it into two polygons, one on each side of the antimeridian.

    Args:
        poly: A polygon in (longitude, latitude) order.

    Returns:
        poly_correct: Corrected polygon / multipolygon.
    """
    lonlat = np.array(poly.exterior.coords.xy)
    lon_jumps = np.where(np.abs(np.diff(lonlat[0])) > 180)[0]
    # if we have 0 jumps, check this polygon is valid, then return it
    if lon_jumps.shape[0] == 0:
        if not poly.is_valid:
            raise ValueError(
                "This polygon appears to be malformed for reasons unrelated to the antimeridian."
            )
        return poly
    # if we have 1 jump, this granule includes a pole
    if lon_jumps.shape[0] == 1:
        latlon1 = lonlat[:, lon_jumps][::-1]
        latlon2 = lonlat[:, lon_jumps + 1][::-1]
        lat_am = get_antimeridian_intersection(latlon1, latlon2)
        # add a point at the antimeridian and a ring of points near the pole
        lon_ring = -np.arange(-180, 181, 60) * np.sign(lonlat[0, lon_jumps[0]])
        polar_lat = np.repeat(np.sign(lat_am) * 89.999, lon_ring.shape[0])
        lonlat_pole = np.concatenate(
            [
                lonlat[:, : lon_jumps[0] + 1],
                np.stack(
                    [
                        lon_ring,
                        polar_lat,
                    ],
                    axis=0,
                ),
                lonlat[:, lon_jumps[0] + 1 :],
            ],
            axis=1,
        )
        return Polygon(lonlat_pole.T)
    # if we have 2 jumps, this granule crosses the dateline but does not include a pole
    elif lon_jumps.shape[0] == 2:
        latlon1 = lonlat[:, lon_jumps][::-1]
        latlon2 = lonlat[:, lon_jumps + 1][::-1]
        lat_am = get_antimeridian_intersection(latlon1, latlon2)
        # split into two polygons at the antimeridian
        lonlat_split1 = np.concatenate(
            [
                lonlat[:, : lon_jumps[0] + 1],
                np.stack(
                    [
                        np.zeros_like(lat_am) + np.sign(lonlat[0, lon_jumps[0]]) * 180,
                        lat_am,
                    ],
                    axis=0,
                ),
                lonlat[:, lon_jumps[1] + 1 :],
            ],
            axis=1,
        )
        lonlat_split2 = np.concatenate(
            [
                np.stack(
                    [
                        np.zeros_like(lat_am) + np.sign(lonlat[0, lon_jumps[1]]) * 180,
                        lat_am,
                    ],
                    axis=0,
                )[:, ::-1],
                lonlat[:, lon_jumps[0] + 1 : lon_jumps[1] + 1],
            ],
            axis=1,
        )
        return MultiPolygon([Polygon(lonlat_split1.T), Polygon(lonlat_split2.T)])
    else:
        raise ValueError(
            "Was not expecting a polygon with more than 2 jumps "
            "in longitude exceeding 180 degrees."
        )


def correct_linestring(line: LineString) -> LineString | MultiLineString:
    """Correct a geospatial line string. If it crosses the antimeridian (dateline),
    split it into two line strings, one on either side.

    Args:
        line: A line string in (longitude, latitude) order.

    Returns:
        line_correct: Corrected line string / multi-line string.
    """
    lonlat = np.array(line.coords.xy)
    lon_jumps = np.where(np.abs(np.diff(lonlat[0])) > 180)[0]
    if lon_jumps.shape[0] == 0:
        return line
    # add points at the antimeridian
    latlon1 = lonlat[:, lon_jumps][::-1]
    latlon2 = lonlat[:, lon_jumps + 1][::-1]
    lat_am = get_antimeridian_intersection(latlon1, latlon2)

    lonlat_split1 = np.stack([np.sign(latlon1[1]) * 180, lat_am], axis=1)
    lonlat_split2 = np.stack([np.sign(latlon2[1]) * 180, lat_am], axis=1)

    coords = [c.T for c in np.split(lonlat, lon_jumps + 1, axis=1)]
    for idx in range(lon_jumps.shape[0]):
        coords[idx] = np.concatenate([coords[idx], lonlat_split1[idx, None]], axis=0)
        coords[idx + 1] = np.concatenate(
            [lonlat_split2[idx, None], coords[idx + 1]], axis=0
        )
    return MultiLineString([LineString(c) for c in coords])


def get_centered_latlon(
    lat: npt.NDArray[np.float32 | np.float64],
    lon: npt.NDArray[np.float32 | np.float64],
) -> tuple[float, float]:
    """Get latitude and longitude arrays centered so that the mean longitude is 0.
    To resolve ambiguity arising from the circularity of longitude, this function
    assumes the maximum range of longitudes when properly centered is less than 180.

    Args:
        lat: Latitude array.
        lon: Longitude array.

    Returns:
        central_lat: Centered latitudes.
        central_lon: Centered longitudes.
    """
    crosses_180 = (np.nanmax(lon) - np.nanmin(lon)) > 180
    if crosses_180:
        central_lon = (np.nanmean((lon + 180) / 360) - 180).item()
    else:
        central_lon = np.nanmean(lon).item()
    central_lat = np.nanmean(lat).item()
    return central_lat, central_lon


def central_latlon_to_rot_mtx(
    central_lat: float,
    central_lon: float,
) -> npt.NDArray[np.float64]:
    """Get a rotation matrix from Earth-centered rotating to center on the specified
    latitude and longitude.

    Args:
        central_lat: Centered latitude, in degrees.
        central_lon: Centered longitude, in degrees.

    Returns:
        rot_mtx: Rotation matrix.
    """
    theta = -central_lon * np.pi / 180
    rot_mtx_z = np.array(
        [
            [np.cos(theta), -np.sin(theta), 0],
            [np.sin(theta), np.cos(theta), 0],
            [0, 0, 1],
        ]
    )
    theta = central_lat * np.pi / 180
    rot_mtx_y = np.array(
        [
            [np.cos(theta), 0, np.sin(theta)],
            [0, 1, 0],
            [-np.sin(theta), 0, np.cos(theta)],
        ]
    )
    return rot_mtx_y @ rot_mtx_z


def get_centering_function(
    lat: npt.NDArray[np.float32 | np.float64],
    lon: npt.NDArray[np.float32 | np.float64],
) -> Callable:
    """Get a function that rotates subsequent lat/lon arrays to the central lat/lon of
    the provided lat/lon arrays.

    Args:
        lat: Latitude array, in degrees.
        lon: Longitude array, in degrees.

    Returns:
        centering_function: Function that rotates lat/lon arrays to center on the
            provided lat/lon array.
    """
    central_lat, central_lon = get_centered_latlon(lat, lon)
    rot_mtx = central_latlon_to_rot_mtx(central_lat, central_lon)[None]

    def centering_function(
        lat: npt.NDArray[np.float32 | np.float64],
        lon: npt.NDArray[np.float32 | np.float64],
    ) -> tuple[
        npt.NDArray[np.float32 | np.float64],
        npt.NDArray[np.float32 | np.float64],
    ]:
        cos_lat = np.cos(np.radians(lat))
        cos_lon = np.cos(np.radians(lon))
        sin_lat = np.sin(np.radians(lat))
        sin_lon = np.sin(np.radians(lon))
        xyz = np.stack([cos_lat * cos_lon, cos_lat * sin_lon, sin_lat], axis=-1)
        shp = xyz.shape
        xyz = np.reshape(xyz, (-1, 3))
        xyz_rot = np.zeros_like(xyz)
        num_batches = int(np.ceil(lat.size / 1e6))
        for i in range(num_batches):
            start = int(i * 1e6)
            end = int(min(lat.size, (i + 1) * 1e6))
            xyz_rot[start:end] = (rot_mtx @ xyz[start:end, :, None])[..., 0]
        xyz_rot = np.reshape(xyz_rot, shp)
        lat_rot = np.degrees(np.arcsin(xyz_rot[..., 2]))
        lon_rot = np.degrees(np.arctan2(xyz_rot[..., 1], xyz_rot[..., 0]))
        return lat_rot, lon_rot

    return centering_function
