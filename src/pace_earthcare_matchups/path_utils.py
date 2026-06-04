import os
from pathlib import Path

from pystac.item import Item

PATH_ROOT = (Path(__file__) / ".." / ".." / "..").resolve()
PATH_DATA = Path(os.getenv("PACE_EARTHCARE_DATA_PATH", PATH_ROOT / "data")).resolve()
PATH_TOKEN = Path(os.getenv("ESA_MAAP_TOKEN_PATH", PATH_ROOT / "token.txt")).resolve()


def get_path(obj: object) -> Path:
    """Get the local path of a serializable object.

    Serializable objects have uniform path specifications defined here, relative to a
    configurable root folder. The root folder defaults to ``{repo_root}/data``, or can
    be set via the ``PACE_EARTHCARE_DATA_PATH`` environment variable.

    :param obj: An object whose local path should be resolved. Currently supports
        ``pystac.item.Item`` (EarthCARE) and any object with a ``filepath`` attribute
        e.g., ``Granule`` (PACE).
    :returns: Local path of the object.
    """
    if isinstance(obj, Item):
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
        assert hasattr(obj, "filepath")
        return getattr(obj, "filepath")
