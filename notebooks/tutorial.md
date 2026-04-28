---
jupytext:
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.19.1
kernelspec:
  display_name: Python 3 (ipykernel)
  language: python
  name: python3
---

```{code-cell} ipython3
%load_ext autoreload
%autoreload 2
```

```{code-cell} ipython3
from datetime import datetime, timedelta

from pace_earthcare_matchups.matchup import get_matchups
from pace_earthcare_matchups.plotting import plot_matchups
```

### Define some search parameters
The only mandatory arguments for `get_matchups` are the PACE and EarthCARE shortnames, as well as the `temporal` argument.

Let's look for overlaps between a PACE OCI level 2 product and two different EarthCARE products. The dates here are arbitrary.

```{code-cell} ipython3
TIME_START = datetime(year=2025, month=9, day=1, hour=0, minute=0, second=0)
TIME_END = datetime(year=2025, month=9, day=2, hour=0, minute=0, second=0)

matchups = get_matchups(
    shortname_pace="PACE_OCI_L2_CLOUD",
    shortnames_earthcare=["ATL_CTH_2A", "AM__CTH_2B"],
    temporal=(TIME_START, TIME_END),
)

plot_matchups(matchups, figsize=(12, 4))
```

### Additional search arguments
There are additional options when searching for matchups!

1. `time_offset`: Different scientific applications have different requirements for the simultaneity of data. Clouds, for example, can move quickly. For matchups of cloud products you might wish to only get matchups where the EarthCARE file's time window overlaps the PACE file's time window, which is the default behavior of `get_matchups`. Ocean matchups typically do not require such strict time windows. You can control the strictness of the time overlaps with `time_offset`, which defines a padding that is applied to the time window of every PACE file before searching for matching EarthCARE files. Here, for demonstration, let's pad the window by 5 minutes.
2. `bbox`: PACE and EarthCARE experience near-simultaneous overlaps in the southern mid-latitudes. Let's limit our search to only PACE granules that include some data between 50$^\circ$ south and the equator, using the `bbox` argument.
3. `limit`: This argument controls how many matchups the function will retrieve before halting. The default is 10, but here we'll set it to just 5.
4. `filter_fn`: If the above options are not powerful enough for you, `filter_fn` allows you to provide a function that screens matchups before they are accepted. It operates on a single `Matchup` and returns `True` if it should be kept, `False` otherwise. Here we'll only keep a matchup if the ATLID track contains at least 2000 total observations within the PACE file's geospatial extent, added across all matched `ATL_CTH_2A` files.

```{code-cell} ipython3
# Keep southern hemisphere, non-polar overlaps only
BBOX = (-180, -50, 180, 0)  # W, S, E, N order

# get a time offset to pad each PACE granule's time window
TIME_OFFSET = timedelta(minutes=5)

def matchup_has_good_overlap(matchup):
    atlid_mask_sum = 0
    for match in matchup.matches_earthcare:
        match_filename = match.filepath_earthcare.name
        product = match_filename[5:15]
        if product == "ATL_CTH_2A":
            atlid_mask_sum += match.mask.sum()
    return atlid_mask_sum >= 2000

matchups = get_matchups(
    shortname_pace="PACE_OCI_L2_CLOUD",
    shortnames_earthcare=["ATL_CTH_2A", "AM__CTH_2B"],
    temporal=(TIME_START, TIME_END),
    bbox=BBOX,
    time_offset=TIME_OFFSET,
    limit=5,
    filter_fn=matchup_has_good_overlap,
)

plot_matchups(matchups, figsize=(12, 4), fig_filepath="../assets/matchup_example.png")
```

### Loading matchups from disk
Once you've found matchups, they are automatically saved to disk. You can use PEM to see all the matchups you have on disk like so:

```{code-cell} ipython3
from pace_earthcare_matchups.matchup import get_all_matchup_paths, load_matchup

matchups_from_disk = [load_matchup(p) for p in get_all_matchup_paths()]
print(f"You have {len(matchups_from_disk)} matchups saved.")
```

### Working with matchups:
Let's do something a bit more advanced now, and start working with the actual data products. We'll need a couple libraries to read the data. First, let's filter the matchups we just loaded to only the ones using the PACE L2 cloud product. If you haven't changed anything, this should give you the same matchups from before.

```{code-cell} ipython3
matchups_cld = [m for m in matchups_from_disk if m.shortname_pace == "PACE_OCI_L2_CLOUD"]
```

Next, let's just make sure we can load the first PACE granule. We'll need the netCDF4 library, though you can also use xarray if you're familiar with that. Every matchup has a `filepath_pace` attribute pointing to its local path. You should see the top-level netCDF information displayed after running the next cell.

```{code-cell} ipython3
import netCDF4

matchup = matchups_cld[0]  # take just the first matchup as an example
data_pace = netCDF4.Dataset(matchup.filepath_pace)
data_pace
```

Now, let's check that we can load the EarthCARE data as well. We'll use h5py. Every matchup has a list of matching EarthCARE data called `matches_earthcare`, and every one of these matches has a `filepath_earthcare`.

```{code-cell} ipython3
import h5py

data_earthcare = h5py.File(matchup.matches_earthcare[0].filepath_earthcare)
data_earthcare
```

PEM conveniently computes a geospatial overlap mask between a PACE granule and each matching EarthCARE file. We retrieved two EarthCARE products: one from ATLID, and one which combines ATLID and MSI. The ATLID masks will be 1-dimensional, corresponding the ATLID curtain, and the ATLID/MSI combined masks will be 2-dimensional images.

```{code-cell} ipython3
for match in matchup.matches_earthcare:
    print(match.mask.shape)
```

Here's a simple plot with the PACE OCI cloud-top height on the left, and the ATLID / MSI combined cloud-top height on the right:

```{code-cell} ipython3
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import numpy as np


# helper function to fill hdf5 data with NaN at its fill value
def read_h5_dataset(dataset):
    arr = dataset[()]
    arr[arr == dataset.fillvalue] = np.nan
    return arr

# get the PACE cloud-top height (convert km -> m)
cth_pace = data_pace["geophysical_data/cth"][()].filled(fill_value=np.nan) * 1000
lat_pace = data_pace["navigation_data/latitude"][()].data
lon_pace = data_pace["navigation_data/longitude"][()].data
# rotate by 180 if needed
rotate_180 = np.ptp(lon_pace) > 350
if rotate_180:
    lon_pace = lon_pace % 360 - 180

# keep just the AM__CTH_2B combined CTH files
data_am = []
for match in matchup.matches_earthcare:
    data_earthcare = h5py.File(match.filepath_earthcare)
    if data_earthcare["HeaderData/FixedProductHeader/File_Type"][()].decode() == "AM__CTH_2B":
        bounds = match.get_earthcare_bounds()
        ext = np.array(bounds.exterior.coords)
        if rotate_180:
            ext[..., 0] = ext[..., 0] % 360 - 180
        data_am.append((ext, data_earthcare))

fig, axs = plt.subplots(1, 2, figsize=(10, 10), subplot_kw={"projection": ccrs.PlateCarree()})
axs[0].pcolormesh(lon_pace, lat_pace, cth_pace, clim=(0, 2e4), shading="gouraud")
for ext, data_earthcare in data_am:
    lat = data_earthcare["ScienceData/latitude"][()]
    lon = data_earthcare["ScienceData/longitude"][()]
    if rotate_180:
        lon = lon % 360 - 180
    cth_msi = read_h5_dataset(data_earthcare["ScienceData/cloud_top_height_MSI"])
    axs[1].pcolormesh(lon, lat, cth_msi, clim=(0, 2e4), shading="gouraud")
    axs[0].plot(
        ext[..., 0],
        ext[..., 1],
        linewidth=2,
    )[0]
axs[0].set_title(matchup.filepath_pace.stem)
axs[1].set_title("Matched MSI cloud-top height")
plt.show();
```
