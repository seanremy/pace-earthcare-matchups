from maap.maap import MAAP
from maap.Result import Granule

from pace_earthcare_matchups.matchup import CMR_HOST
from pace_earthcare_matchups.metadata_utils import (
    get_datetime_range_from_granule,
    get_datetime_range_maap,
)


def get_simultaneous_pace_product(granule: Granule, shortname_pace: str) -> Granule:
    """Using one PACE granule, get a different product with the same timestamp.

    Args:
        granule: A PACE granule's MAAP metadata.
        shortname_pace: The shortname of a different PACE product to retrieve.

    Returns:
        result: MAAP metadata of a PACE granule co-occurring with the provided granule.
    """
    maap = MAAP()
    dt_range = get_datetime_range_from_granule(granule)
    result = maap.searchGranule(
        cmr_host=CMR_HOST,
        short_name=shortname_pace,
        temporal=get_datetime_range_maap(*get_datetime_range_from_granule(granule)),
    )[0]
    dt_range_result = get_datetime_range_from_granule(result)
    assert abs((dt_range_result[0] - dt_range[0]).total_seconds()) < 10
    assert abs((dt_range_result[1] - dt_range[1]).total_seconds()) < 10
    return result
