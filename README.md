# nrcan-hfe
Codes to feed data for flood events into NRCan's HFE database


## Python environment

```python
pyenv virtualenv 3.8.5 env-3.8.5-nrcan
pyenv activate env-3.8.5-nrcan
```

* `pip install GDAL==$(gdal-config --version) --global-option=build_ext --global-option="-I/usr/include/gdal"`
* `pip install argparse`
* `pip install numpy`
* `pip install geopandas`
* `pip install netCDF4`
* `pip install pygrib`
