# Historical Flood Event Python Toolkit

This library is a collection of tools helpful to analyse flood events of the Historial Flood Event Database issued by Natural Resources Canada (NRCan). A documentation of the tools and their usage can be found [here](https://github.com/julemai/nrcan-hfe/blob/main/doc/HFE_database_2022Q03_documentation.pdf). Below are instructions on how to setup a Python environment with requirements to run tools of this collection. It also conatins instructions on how to test the various functions. Questions can be directed to <a href="mailto:juliane.mai@uwaterloo.ca">Julie Mai</a> or through the issues section [here](https://github.com/julemai/nrcan-hfe/issues).

<p align="center">
   <img alt="Umbrella icons created by Freepik - Flaticon. See here: https://www.flaticon.com/free-icons/umbrella" src="https://github.com/julemai/nrcan-hfe/blob/main/doc/logo/HFE_logo.png" width="65%" />
</p>

## Citation
<a href="https://doi.org/10.5281/zenodo.7121539"><img src="https://zenodo.org/badge/DOI/10.5281/zenodo.7121539.svg" alt="DOI"></a><br>
The code is published under [Zenodo](https://doi.org/10.5281/zenodo.7121538). Please cite the following if you are using it:

Mai, Juliane. (2022). <br>
Historical Flood Event Python Toolkit (v1.0). <br>
Zenodo. <br>
https://doi.org/10.5281/zenodo.7121539


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

## Getting started

Analyse two single-point occurrences listed in `data/hfe/historical_flood.json`:
```python
# Maybe load some modules on your system?! (see, for example, the ones listed for Graham under D.)

# Load your Python environment.
pyenv activate env-3.8.5-nrcan

# Analyse occurrences.
# Note: These two feaures have been selected because they require data from Geomet.
#       If you pick other ones, you might need to request RDRS v2.1 from CaSPAr first.
python src/analyse_occurrence.py --ifeatures "873, 1092" --bbox_buffer 0.5 --dates_buffer 5.0,0.0 --tmpdir "/tmp/"
```

Analyse two multi-point events listed in `data/hfe/historical_flood_event.json`:
```python
# Maybe load some modules on your system?! (see, for example, the ones listed for Graham under D.)

# Load your Python environment.
pyenv activate env-3.8.5-nrcan

# Analyse events.
# Note: These two feaures have been selected because they require data from Geomet.
#       If you pick other ones, you might need to request RDRS v2.1 from CaSPAr first.
python src/analyse_event.py --ifeatures "0, 178" --bbox_buffer 0.5 --dates_buffer 5.0,0.0 --tmpdir "/tmp/"
```

Please refer to the [documentation](https://github.com/julemai/nrcan-hfe/blob/main/doc/HFE_database_2022Q03_documentation.pdf), for an exaplantion of the arguments passed to these functions as well as a documentation of the entire toolkit.

