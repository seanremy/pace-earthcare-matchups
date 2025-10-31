---
jupytext:
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.18.1
kernelspec:
  display_name: Python 3 (ipykernel)
  language: python
  name: python3
---

# MAAP search basics

+++

1. [X] Get all PACE L1B-L2 short names
2. [X] Get all EarthCARE L1B-L2 short names
3. [X] Format the PACE geospatial bounds for STAC
4. [X] Search EarthCARE data overlapping a PACE product
5. [X] Download EarthCARE files
6. [ ] Get overlaps for arbitrary combinations

```{code-cell} ipython3
%load_ext autoreload
%autoreload 2
```

```{code-cell} ipython3
!pip install cartopy netCDF4
```

```{code-cell} ipython3
from collections import defaultdict
from datetime import datetime, timedelta
import json
import os
from pathlib import Path
import requests
import sys
import warnings

import boto3
import cartopy.crs as ccrs
import fsspec
import h5py
from maap.maap import MAAP
from maap.Result import Granule
import matplotlib.pyplot as plt
import netCDF4
import numpy as np
from pystac_client import Client
from pystac.item import Item
from requests.exceptions import RequestException
from shapely import (
    GeometryCollection,
    LineString,
    MultiLineString,
    MultiPolygon,
    Polygon,
    Point,
    to_geojson,
)
from tqdm.notebook import tqdm
import xarray as xr

if os.getcwd().split("/")[-1] == "notebooks":
    os.chdir("..")
PATH_HERE = Path(os.getcwd())
PATH_SRC = PATH_HERE / "src" / "pace_earthcare_matchups"
if str(PATH_SRC) not in sys.path:
    sys.path.insert(0, str(PATH_SRC))
from geospatial_utils import correct_linestring, correct_polygon
from matchup import Matchup, MatchEarthCARE
from metadata_utils import (
    geometry_from_item,
    get_datetime_range_from_granule,
    get_datetime_range_maap,
    polygon_from_granule,
    UTC,
)
from supported_products import (
    EARTHCARE_SHORTNAMES,
    PACE_SHORTNAMES,
    PACEPAX_SHORTNAMES,
)


# constants
PATH_DATA = PATH_HERE / "data"
CMR_HOST = "cmr.earthdata.nasa.gov"
ESA_CATALOGUE = "https://catalog.maap.eo.esa.int/catalogue/"
_ESA_MAAP_TOKEN = open(PATH_HERE / "token.txt").read()
FS_ESA = fsspec.filesystem("https", headers={"Authorization": f"Bearer {_ESA_MAAP_TOKEN}"})


# arguments
SHORTNAME_EARTHCARE = "ATL_EBD_2A"
TIME_START = datetime(year=2025, month=6, day=11, hour=0, minute=0, second=0, tzinfo=UTC)
TIME_END = datetime(year=2025, month=9, day=19, hour=23, minute=59, second=59, tzinfo=UTC)
TIME_OFFSET = timedelta(minutes=5)
# TIME_START = datetime(year=2025, month=9, day=2, hour=23, minute=58, second=0, tzinfo=UTC)
# TIME_END = datetime(year=2025, month=9, day=2, hour=23, minute=59, second=59, tzinfo=UTC)
# TIME_OFFSET = timedelta(minutes=0)
```

### Get the clients to search both the NASA and ESA MAAPs

```{code-cell} ipython3
PACE_SHORTNAMES
```

```{code-cell} ipython3
maap = MAAP()
client_esa = Client.open(ESA_CATALOGUE)

# # get all L1 / L2 collections:
# pace_collections_all = {int(l): maap.searchCollection(cmr_host=CMR_HOST, platform='PACE', processing_level_id=l) for l in ['1', '2']}
# # search for just the supported shortnames (this takes a while)
# pace_collections = {sn: maap.searchCollection(cmr_host=CMR_HOST, short_name=sn)[0] for sn in PACE_SHORTNAMES}
```

### Search in a loop to make sure you get all of the results

```{code-cell} ipython3
matches = []

# use a tiny bbox in the center of the Andes to get granules with high elevations
bbox_andes = '-69,-20,-67,-18'  # [W, S, E, N]

results_pace = maap.searchGranule(
    cmr_host=CMR_HOST,
    short_name="PACE_HARP2_L1B_SCI",
    # concept_id="C3555841897-OB_CLOUD",
    temporal = get_datetime_range_maap(TIME_START, TIME_END),
    # bounding_box=bbox_andes,  # TODO: remove
    limit=20,
)
time_start = TIME_START
while len(results_pace) > 1:
    print(f"Looking through {len(results_pace)} PACE results")
    for result_pace in results_pace:
        granule_filename = result_pace['Granule']['DataGranule']['ProducerGranuleId']
        pace_poly = correct_polygon(polygon_from_granule(result_pace))
        if not pace_poly.is_valid:
            warnings.warn(f"Polygon for {granule_filename} is not valid, even after correction")
            continue
        dt_range = get_datetime_range_from_granule(result_pace)
        elapsed_seconds = (dt_range[1] - dt_range[0]).total_seconds()
        if elapsed_seconds > 299.0:
            warnings.warn(f"Granule's duration of {elapsed_seconds} seconds is too long.")
            continue
        dt_range_offset = get_datetime_range_from_granule(result_pace, TIME_OFFSET)
        results_esa = client_esa.search(
            collections=['EarthCAREL1Validated_MAAP', 'EarthCAREL2Validated_MAAP'],
            datetime=dt_range_offset,
            intersects=to_geojson(pace_poly),
            method="GET",
            filter=f"productType = '{SHORTNAME_EARTHCARE}'",
        )
        num_matched = results_esa.matched()
        print(f"Matched {num_matched}")
        if num_matched > 0:
            matches.append((result_pace, list(results_esa.items())))
        time_start = max(time_start, dt_range[1])

    time_start += timedelta(seconds=1)
    if time_start > TIME_END:
        break
    results_pace = maap.searchGranule(
        cmr_host=CMR_HOST,
        short_name="PACE_HARP2_L1B_SCI",
        # concept_id="C3555841897-OB_CLOUD",
        temporal=get_datetime_range_maap(time_start, TIME_END),
        # bounding_box=bbox_andes,  # TODO: remove
        limit=20,
    )
```

### Filter the matches into real and broken, display the broken ones

+++

### List out filenames and display all of the matches

```{code-cell} ipython3
earth_poly = Polygon(np.array([
    [-180, -90],
    [180, -90],
    [180, 90],
    [-180, 90],
    [-180, -90],
]))

for match_idx, (result_pace, items_esa) in enumerate(matches):
    pace_poly = correct_polygon(polygon_from_granule(result_pace))
    esa_lines = [correct_linestring(geometry_from_item(item)) for item in items_esa]

    matching_esa_idx = [i for i, l in enumerate(esa_lines) if pace_poly.intersects(l)]
    matching_esa_lines = [esa_lines[i] for i in matching_esa_idx]
    matching_esa_urls = [items_esa[i].assets["enclosure_"].href for i in matching_esa_idx]

    print(match_idx)
    print(result_pace['Granule']['DataGranule']['ProducerGranuleId'])
    print([Path(url).name for url in matching_esa_urls])

    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.coastlines()
    pace_polys = pace_poly.geoms if isinstance(pace_poly, MultiPolygon) else [pace_poly]
    for poly in pace_polys:
        ax.plot(*poly.exterior.xy)
    for esa_line in matching_esa_lines:
        lines = esa_line.geoms if isinstance(esa_line, MultiLineString) else [esa_line]
        for line in lines:    
            ax.plot(*line.coords.xy)
    plt.show()
    print('----------------------------')
```

### Download an EarthCARE item

```{code-cell} ipython3
result_pace, items_esa = real_matches[5]

def download_earthcare_item(item: Item) -> str:
    url_h5 = item.assets["enclosure_"].href
    filename = item.assets["enclosure_"].title
    response = requests.get(
        url_h5,
        headers={"Authorization": f"Bearer {_ESA_MAAP_TOKEN}"},
        stream=True,
    )
    response.raise_for_status()
    content_length = int(response.headers.get('content-length', 0))
    path_outdir = PATH_DATA / item.properties["platform"] / item.properties["product:type"]
    path_outfile = path_outdir / filename
    
    os.makedirs(path_outdir, exist_ok=True)
    with open(path_outfile, "wb") as outfile:
        pbar = tqdm(
            desc=f"{filename}",
            total=content_length // (2 ** 20),
            unit="MB",
        )
        # 10 MB chunks
        for chunk in response.iter_content(chunk_size=10 * (2 ** 20)):
            pbar.update(outfile.write(chunk) // (2 ** 20))

    return path_outfile

for item in items_esa:
    download_earthcare_item(item)
```

### HARP2_L1B overlap with ATL_EBD_2A using centered lat/lon

```{code-cell} ipython3
FILENAME_EBD = items_esa[0].assets['enclosure_'].title
PATH_EBD_FILE = PATH_DATA / "EarthCARE" / "ATL_EBD_2A" / FILENAME_EBD
ec_h5 = h5py.File(PATH_EBD_FILE)

FILENAME_PACE = result_pace['Granule']['DataGranule']['ProducerGranuleId']
result_pace.getData(PATH_DATA / "PACE" / "PACE_HARP2_L1B_SCI")
harp2_l1b_data = netCDF4.Dataset(PATH_DATA / "PACE" / "PACE_HARP2_L1B_SCI" / FILENAME_PACE)
```

```{code-cell} ipython3
# get HARP2 bounding polygon from netcdf data
lat_harp2 = harp2_l1b_data['geolocation_data/latitude'][:].filled(fill_value=np.nan)
lon_harp2 = harp2_l1b_data['geolocation_data/longitude'][:].filled(fill_value=np.nan)
seconds_harp2 = harp2_l1b_data['navigation_data/seconds_of_day'][:].filled(fill_value=np.nan)
# TO DO: handle seconds correctly if granule crosses the day line
```

```{code-cell} ipython3
lat_ec = ec_h5['ScienceData/latitude'][:]
lon_ec = ec_h5['ScienceData/longitude'][:]
time_ec = ec_h5['ScienceData/time'][:]
```

```{code-cell} ipython3
# detect if crossing 180 lon:
cross_180 = (np.nanmax(lon_harp2) - np.nanmin(lon_harp2)) > 180
if cross_180:
    central_lon = np.nanmean((lon_harp2 + 180) / 360) - 180
else:
    central_lon = np.nanmean(lon_harp2)
central_lat = np.nanmean(lat_harp2)

def central_latlon_to_rot_mtx(central_lat, central_lon):
    """Get a rotation matrix to apply to spherical cartesian points to center
    on a provided latitude and longitude.
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

def get_centering_function(central_lat, central_lon):
    rot_mtx = central_latlon_to_rot_mtx(central_lat, central_lon)[None]
    def centering_function(lat, lon):
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

centering_fn = get_centering_function(central_lat, central_lon)
# lat_rot_harp2, lon_rot_harp2 = centering_fn(lat_harp2, lon_harp2)
lat_rot_ec, lon_rot_ec = centering_fn(lat_ec, lon_ec)
```

```{code-cell} ipython3
# select some pixel indices along each edge of the HARP2 granule to form a makeshift polygon
num_pts_per_edge = 10
vspace = (lat_harp2.shape[1] - 1) / (num_pts_per_edge - 1)
hspace = (lat_harp2.shape[2] - 1) / (num_pts_per_edge - 1)
pts_left = []
pts_bot = []
pts_right = []
pts_top = []
for i in range(1, num_pts_per_edge - 1):
    pts_left.append((int(i * vspace), 0))
    pts_bot.append((lat_harp2.shape[1] - 1, int(i * hspace)))
    pts_right.append((lat_harp2.shape[1] - 1 - int(i * vspace), lat_harp2.shape[2] - 1))
    pts_top.append((0, lat_harp2.shape[2] - 1 - int(i * hspace)))
poly_idx = [(0, 0)]
poly_idx += pts_left
poly_idx += [(lat_harp2.shape[1] - 1, 0)]
poly_idx += pts_bot
poly_idx += [(lat_harp2.shape[1] - 1, lat_harp2.shape[2] - 1)]
poly_idx += pts_right
poly_idx += [(0, lat_harp2.shape[2] - 1)]
poly_idx += pts_top
poly_idx += [(0, 0)]

# get and rotate the lat/lon at these polygon locations, then take the mean for each angle
lat_poly = lat_harp2[:, [p[0] for p in poly_idx], [p[1] for p in poly_idx]]
lon_poly = lon_harp2[:, [p[0] for p in poly_idx], [p[1] for p in poly_idx]]
lat_rot_poly, lon_rot_poly = centering_fn(lat_poly, lon_poly)
lat_rot_poly, lon_rot_poly = np.nanmean(lat_rot_poly, axis=0), np.nanmean(lon_rot_poly, axis=0)

# convert to shapely polygon
harp2_poly = Polygon(zip(lon_rot_poly, lat_rot_poly))

# get mask of where the EarthCARE track is in the HARP2 granule
harp2_contains = np.vectorize(lambda p: harp2_poly.contains(Point(p)), signature='(n)->()')
ec_in_granule = harp2_contains(np.stack([lon_rot_ec, lat_rot_ec], axis=-1))
# get the longest contiguous
ec_in_granule_diff = np.diff(ec_in_granule.astype(int), prepend=np.zeros(1), append=np.zeros(1))
ec_in_granule_start_idxs = np.where(ec_in_granule_diff == 1)[0]
ec_in_granule_end_idxs = np.where(ec_in_granule_diff == -1)[0]
ec_contig_lens = ec_in_granule_end_idxs - ec_in_granule_start_idxs
longest_contig_idx = np.argmax(ec_contig_lens)
ec_in_granule_start_idx = ec_in_granule_start_idxs[longest_contig_idx]
ec_in_granule_end_idx = ec_in_granule_end_idxs[longest_contig_idx]

display(GeometryCollection([
    harp2_poly,
    LineString(np.stack([
        lon_rot_ec,
        lat_rot_ec
    ], axis=-1)[ec_in_granule_start_idx:ec_in_granule_end_idx]),
]))
```

### Write matchup to file

```{code-cell} ipython3
matchup = Matchup(
    filename_PACE = FILENAME_PACE,
    filename_EarthCARE = FILENAME_EBD,
    idx_start_EarthCARE = ec_in_granule_start_idx.item(),
    idx_end_EarthCARE = ec_in_granule_end_idx.item(),
)
matchup_filename = FILENAME_PACE.split('.')[1] + ".json"
matchup_filepath = PATH_DATA / "matchups" / "HARP2_L1B_ATL_EBD_2A" / matchup_filename
json.dump(vars(matchup), open(matchup_filepath, 'w'))
```

```{code-cell} ipython3

```

```{code-cell} ipython3

```

```{code-cell} ipython3

```

```{code-cell} ipython3
# # get all L1 / L2 collections:
# pace_collections_all = {int(l): maap.searchCollection(cmr_host=CMR_HOST, platform='PACE', processing_level_id=l) for l in ['1', '2']}
# # search for just the supported shortnames (this takes a while)
# pace_collections = {sn: maap.searchCollection(cmr_host=CMR_HOST, short_name=sn)[0] for sn in PACE_SHORTNAMES}
```

```{code-cell} ipython3

```

```{code-cell} ipython3

```

```{code-cell} ipython3

```

```{code-cell} ipython3

```
