# nrcan-hfe
Codes to provide precipitation data of flood events into NRCan's HFE database.


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

### C. Using conda

Helpful commands can be found [here](https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html). Syntax on how to install various packages can be found [here](https://anaconda.org/anaconda/xarray) by using search bar to find packages.

```bash
conda create -n env-conda-3.8 python=3.8
conda activate env-conda-3.8

conda install -n env-conda-3.8 -c conda-forge gdal 
conda install -n env-conda-3.8 -c conda-forge argparse
conda install -n env-conda-3.8 -c conda-forge numpy
conda install -n env-conda-3.8 -c conda-forge scipy
conda install -n env-conda-3.8 -c conda-forge geopandas
conda install -n env-conda-3.8 -c conda-forge netcdf4
conda install -n env-conda-3.8 -c conda-forge xarray 
conda install -n env-conda-3.8 -c conda-forge pygrib
conda install -n env-conda-3.8 -c conda-forge matplotlib
conda install -n env-conda-3.8 -c conda-forge basemap
conda install -n env-conda-3.8 -c conda-forge basemap-data-hires
conda install -n env-conda-3.8 -c conda-forge pytest

# optional
conda install -n env-conda-3.8 -c conda-forge ipython 
```

### D. Setting up a Pythin environment on Compute Canada systems (tested on Graham)

```bash
cd nrcan-hfe
mkdir env-3.8

module purge
module load StdEnv/2020 netcdf gcc/9.3.0 gdal/3.0.4
module load mpi4py proj
module load python/3.8

virtualenv --no-download env-3.8
source env-3.8/bin/activate

pip install --no-index --upgrade pip
pip install numpy --no-index
pip install GDAL --no-index
pip install argparse --no-index
pip install geopandas --no-index
pip install netCDF4 --no-index
pip install scipy --no-index
pip install xarray --no-index
pip install pygrib --no-index
pip install matplotlib --no-index
pip install ipython --no-index

python -m pip install basemap --no-index
pip install -U pytest --no-index
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
python fx_determine_dates.py 
python gx_identify_precipitation_event.py
```

