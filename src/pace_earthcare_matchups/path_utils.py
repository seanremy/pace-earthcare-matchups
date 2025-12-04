import os
from pathlib import Path

from maap.Result import Granule
from pystac.item import Item

PATH_ROOT = Path(os.getenv('PACE_ROOT_PATH', (Path(__file__) / ".." / ".." / "..").resolve()))
PATH_DATA = (PATH_ROOT / "data").resolve()
PATH_TOKEN = (PATH_ROOT / "token.txt").resolve()


def get_path(obj: Granule | Item) -> Path:
    """TODO"""
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
