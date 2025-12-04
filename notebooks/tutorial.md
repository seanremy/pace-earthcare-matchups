---
jupytext:
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.18.1
kernelspec:
  display_name: Python 3
  language: python
  name: python3
---

```{code-cell} ipython3
%load_ext autoreload
%autoreload 2
```

```{code-cell} ipython3
from datetime import datetime, timedelta

from maap.maap import MAAP
from pystac_client import Client

from pace_earthcare_matchups.matchup import get_matchups
from pace_earthcare_matchups.path_utils import PATH_TOKEN
from pace_earthcare_matchups.plotting import plot_matchups
```

```{code-cell} ipython3
TIME_START = datetime(year=2025, month=9, day=1, hour=0, minute=0, second=0)
TIME_END = datetime(year=2025, month=9, day=2, hour=0, minute=0, second=0)

# Keep southern hemisphere, non-polar overlaps only
BBOX = (-180, -50, 180, 0)  # W, S, E, N order

CMR_HOST = "cmr.earthdata.nasa.gov"
_ESA_MAAP_TOKEN = open(PATH_TOKEN).read()
ESA_CATALOGUE = "https://catalog.maap.eo.esa.int/catalogue/"

maap = MAAP()
client_esa = Client.open(ESA_CATALOGUE)
```

```{code-cell} ipython3
matchups = get_matchups(
    maap=maap,
    client_esa=client_esa,
    long_term_token=_ESA_MAAP_TOKEN,
    shortname_pace="PACE_OCI_L2_CLOUD",
    shortnames_earthcare=["ATL_CTH_2A", "AM__CTH_2B"],
    temporal=(TIME_START, TIME_END),
    bbox=BBOX,
    limit=5,
)
```

```{code-cell} ipython3
data_pace = netCDF4.Dataset(matchups[3].filepath_pace)
data_pace["navigation_data/latitude"][:].max()
```

```{code-cell} ipython3
for matchup in matchups:
    print(matchup.filepath_pace)
    for match_earthcare in matchup.matches_earthcare:
        print(f"\t{match_earthcare.filepath_earthcare}")
```

```{code-cell} ipython3
# plot_matchups(matchups)
plot_matchups(matchups, fig_filepath="../assets/matchup_example.png")
```

```{code-cell} ipython3

```
