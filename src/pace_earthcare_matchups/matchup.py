"""TODO"""

from dataclasses import dataclass


@dataclass
class MatchEarthCARE:
    """TODO"""
    filename: str
    idx_start: int
    idx_end: int


@dataclass
class Matchup:
    """TODO"""
    filename_PACE: str
    matches_EarthCARE: list[MatchEarthCARE]
