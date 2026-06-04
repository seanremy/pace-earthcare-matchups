.. PACE / EarthCARE Matchups documentation master file, created by
   sphinx-quickstart on Thu Jun  4 18:42:37 2026.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

PACE / EarthCARE Matchups (PEM) API documentation
=================================================

.. _PACE: https://pace.gsfc.nasa.gov/
.. _EarthCARE: https://earth.esa.int/eogateway/missions/earthcare

**PACE / EarthCARE Matchups (PEM)** is a Python library for simplifying the intercomparison of PACE_ and EarthCARE_ data. With a variety of utilities, PEM simplifies all stages of PACE / EarthCARE intercomparisons, deduplicating coding effort and improving reproducibility.


.. note::
   This project is under active development. Documentation is a work in progress.


Basic Usage
-----------

The most commonly used functions are presented here. Complete API documentation will be available soon.

.. autofunction:: pace_earthcare_matchups.matchup.get_matchups
.. autofunction:: pace_earthcare_matchups.matchup.get_all_matchup_paths
.. autofunction:: pace_earthcare_matchups.matchup.load_matchup
.. autofunction:: pace_earthcare_matchups.plotting.plot_matchup
