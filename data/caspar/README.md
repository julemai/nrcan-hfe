# CaSPAr Data

Common folder to potentially hold data downloaded from CaSPAr.

The shapefile ``RDRS_v21_hfe_outline.zip'' has been used to request
RDRS v2.1 from CaSPAr. Only the precipitation variable based on CaPA
``RDRS_v2.1_A_PR0_SFC'' has been requested for all days available (Jan
1, 1980 to Dec 31, 2018). Each file contains 24 time steps (1pm UTC
until noon UTC the day after).

Example files are provided:
* ``src/test-data/1980011312.nc'' containing hourly precipitation
  fields for Jan 13, 1980 1:00pm until Jan 14, 1980 12:00pm across Canada
* ``src/test-data/1980011412.nc'' containing hourly precipitation
  fields for Jan 14, 1980 1:00pm until Jan 15, 1980 12:00pm across Canada
