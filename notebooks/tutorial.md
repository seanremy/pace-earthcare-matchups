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

### Comparisons
So far we've only worried about getting matchups on the file level. Next, we'll learn how to use PEM to compare our data on the level of individual observations.

Fortunately, PEM makes this easy. All we need is the names of the PACE and EarthCARE variables we want to compare.

We'll start with the `get_comparison_dict` function, which takes a single matchup and a list of PACE variable paths. This function samples the PACE variables of interest to each EarthCARE sensor's frame. For now, we'll just ask for the OCI cloud-top height "cth", which can be found in the "geophysical_data" group.

```{code-cell} ipython3
from pace_earthcare_matchups.compare import get_comparison_dict

matchups_cld = [m for m in matchups_from_disk if m.shortname_pace == "PACE_OCI_L2_CLOUD"]
matchup = matchups_cld[0]  # take just the first matchup
# the next line may take a few seconds to interpolate all the data
comp = get_comparison_dict(
    matchup,
    ["geophysical_data/cth"],
)
```

The retrieved dictionary has our CTH variable resampled for comparison with both of the EarthCARE instruments we specified earlier. For MSI, the array will be 2D; for ATLID, it will be 1D. It also contains the latitude and longitudes of the matched points. Also, it keeps track of the filepath to the matching EarthCARE file.

```{code-cell} ipython3
print("MSI-interpolated shapes:")
for c in comp["AM__CTH_2B"].values():
    print("\t", c["filepath_earthcare"].name)
    print("\t", c["geophysical_data/cth"].shape)
print("ATLID-interpolated shapes:")
for c in comp["ATL_CTH_2A"].values():
    print("\t", c["filepath_earthcare"].name)
    print("\t", c["geophysical_data/cth"].shape)
```

### Plotting
Now we have everything we need to generate some comparison plots. First, comparing OCI and MSI:

```{code-cell} ipython3
import cartopy.crs as ccrs
import h5py
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import numpy as np

# helper function to fill hdf5 data with NaN at its fill value
def read_h5_dataset(dataset):
    arr = dataset[()]
    arr[arr == dataset.fillvalue] = np.nan
    return arr

c = comp["AM__CTH_2B"][0]

data = h5py.File(c["filepath_earthcare"])
cth_msi = read_h5_dataset(data["ScienceData/cloud_top_height_MSI"])[c["start"]:c["end"]] / 1000  # m -> km
cth_pace = c["geophysical_data/cth"]

fig = plt.figure(figsize=(9, 8), constrained_layout=True)
gs = gridspec.GridSpec(62, 62, figure=fig)
ax_oci = fig.add_subplot(gs[:, :20], projection=ccrs.PlateCarree())
ax_msi = fig.add_subplot(gs[:, 20:40], projection=ccrs.PlateCarree())
ax_cbar1 = fig.add_subplot(gs[:, 40:41])
ax_diff = fig.add_subplot(gs[:, 41:61], projection=ccrs.PlateCarree())
ax_cbar2 = fig.add_subplot(gs[:, 61:])

# fig, axs = plt.subplots(1, 3, figsize=(10, 10), subplot_kw={"projection": ccrs.PlateCarree()})
cbar1 = ax_oci.pcolormesh(c["longitude"], c["latitude"], cth_pace, clim=(0, 15), shading="gouraud")
ax_oci.set_title(f"OCI cloud-top height (km)\n{matchup.filepath_pace.stem.split('.')[1]}")
ax_msi.pcolormesh(c["longitude"], c["latitude"], cth_msi, clim=(0, 15), shading="gouraud")
ax_msi.set_title(f"MSI cloud-top height (km)\n{c['filepath_earthcare'].name.split('_')[5]}")
cbar2 = ax_diff.pcolormesh(c["longitude"], c["latitude"], cth_msi - cth_pace, clim=(-2, 2), cmap="RdBu", shading="gouraud")
ax_diff.set_title("MSI $-$ OCI (km)")
for ax in [ax_oci, ax_msi, ax_diff]:
    ax.set_xlim(c["longitude"].min(), c["longitude"].max())
    ax.set_ylim(c["latitude"].min(), c["latitude"].max())
    ax.set_facecolor("black")
fig.colorbar(cbar1, cax=ax_cbar1)
fig.colorbar(cbar2, cax=ax_cbar2)
plt.show();
```

Finally, a line plot comparing OCI with the ATLID CTH:

```{code-cell} ipython3
c = comp["ATL_CTH_2A"][0]

data = h5py.File(c["filepath_earthcare"])
cth_atl = read_h5_dataset(data["ScienceData/ATLID_cloud_top_height"])[c["start"]:c["end"]] / 1000  # m -> km

plt.figure(figsize=(16,  4))
plt.plot(c["latitude"][:, 0], c["geophysical_data/cth"], label="OCI")
plt.plot(c["latitude"][:, 0], cth_atl, label="ATLID")
plt.legend()
plt.title("Cloud-top height (km)")
plt.show()
```
