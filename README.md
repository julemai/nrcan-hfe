# nrcan-hfe
Codes to feed data for flood events into NRCan's HFE database


## Python environment

It is recommended to use `pyenv` to setup a Python environment. An installation guide can be found [here](https://realpython.com/lessons/installing-pyenv/). Additionally, the GDAL library needs to be installed before setting up the Python environment. Please use, for example, `brew install gdal` or your native package installer to get GDAL. Test if it is installed by checking the version number using `gdal-config --version`. Then proceed to setup the Python environment.

### A. Manual setup and installation (recommended since basemap is used)

```python
pyenv virtualenv 3.8.5 env-3.8.5-nrcan
pyenv activate env-3.8.5-nrcan

pip install GDAL==$(gdal-config --version) --global-option=build_ext --global-option="-I/usr/include/gdal"
pip install argparse
pip install numpy
pip install scipy
pip install geopandas
pip install netCDF4
pip install xarray
pip install pygrib
pip install matplotlib
python -m pip install basemap
pip install -U pytest

# optional
pip install ipython
```

### B. Automatic setup and installation (not tested)

```python
pyenv virtualenv 3.8.5 env-3.8.5-nrcan
pyenv activate env-3.8.5-nrcan

pip install -r /path/to/nrcan-hfe/src/requirements.txt
```


## Testing codes

All tests for ```pytest``` are setup under ```src/tests```. Please run the following command to check that all tests pass before pushing any changes:
```python
# check all tests (before push to this repo)
pytest

# check specific test (for development)
pytest src/tests/test_b1_read_geomet_grib2.py
```

All functions have build-in docstring tests that are checked by running the individual scripts. They will result in error messages if any test fails.

```python
python a1_request_geomet_grib2.py 
python a2_request_caspar_nc.py 
python b1_read_geomet_grib2.py 
python b2_read_caspar_nc.py 
python a3_request_hfe_json.py 
python b3_read_hfe_json.py 
python cx_plot_data.py 
python dx_interpolate_data.py 
python ex_determine_bbox.py 
```

