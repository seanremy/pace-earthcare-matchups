import os
from pathlib import Path

from maap.Result import Granule
from pystac.item import Item

PATH_ROOT = (Path(__file__) / ".." / ".." / "..").resolve()
PATH_DATA = Path(os.getenv("PACE_EARTHCARE_DATA_PATH", PATH_ROOT / "data")).resolve()
PATH_TOKEN = (PATH_ROOT / "token.txt").resolve()


def get_path(obj: Granule | Item) -> Path:
    """Get the local path of a serializable object. Serializable objects have uniform
    specifications defined here for their paths, up to a root folder. The root folder
    defaults to "{repo_root}/data", or can be configured with the environment variable
    "PACE_EARTHCARE_DATA_PATH".

    Args:
        obj: An object of a type that has a specification for paths, currently only
            Granule and Item.

    Returns:
        path: Local path of the object.
    """
    if isinstance(obj, Granule):
        plat = obj["Granule"]["Platforms"]["Platform"]
        instrument = plat["Instruments"]["Instrument"]["ShortName"]
        level = obj["Granule"]["Collection"]["ShortName"].split("_")[2]
        filename = obj["Granule"]["DataGranule"]["ProducerGranuleId"]
        return PATH_DATA / "PACE" / instrument / level / filename
    elif isinstance(obj, Item):
        product_type = obj.properties["product:type"]
        title = obj.assets["enclosure_h5"].title
        assert title
        # This is extremely unfortunate, but AWS bucket name length limits
        #   require the truncation of ONE character off of the end of the
        #   standard EarthCARE filenames. Therefore, truncate the unnecessary
        #   "ECA_" prefix off of EarthCARE files.
        filename = title.removeprefix("ECA_")
        return PATH_DATA / "EarthCARE" / product_type / filename
    else:
        raise TypeError(f"No path specification for objects of type {type(obj)}!")
