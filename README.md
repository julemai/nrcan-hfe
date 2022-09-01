# nrcan-hfe
Codes to feed data for flood events into NRCan's HFE database


## Python environment

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

# optional
pip install ipython
```

## Testing codes

All functions have build-in docstring tests that are checked by running the indiviudal scripts. They will result in error messages if any test fails.

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

