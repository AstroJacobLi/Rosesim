# RoSE-Sim: Roman semi-resolved galaxy simulator

A Python package for image simulations of semi-resolved dwarf galaxies for Roman

<!-- insert demo.png  -->
![Concept of Rosesim](demo.png)

## Installation

You can install locally with:
```bash
pip install -e .
```
You need to set `ROSESIM_DATA_PATH` in your environment variables, e.g., add `export ROSESIM_DATA_PATH=/scratch/gpfs/JENNYG/jiaxuanl/Data/SBF/Rosesim/` to `.bashrc`.

Then in Python, you can download the data needed for Rosesim with:
```python
import rosesim
rosesim.fetch_data()
```
The data will be downloaded to the directory specified by `ROSESIM_DATA_PATH`. The following sections describes the data files needed for Rosesim.

### Stellar Population Synthesis Models
The [PARSEC isochrones](https://stev.oapd.inaf.it/cgi-bin/cmd) are required for the stellar population synthesis. I have prepared them for you. The isochrones for Roman are in **Vega** magnitudes. You have to add the zeropoint offsets to convert them to **AB** magnitudes. The zeropoint offsets are stored in `rosesim.Roman_zp_AB_Vega_mist`. Unfortunately, `artpop` doesn't support the PARSEC isochrones to be interpolated (yet), so the available simple stellar populations are:

**Available filters:**  
['F062', 'F087', 'F106', 'F129', 'F158', 'F184', 'F146', 'F213']

**Available log(age / yr):**  
[8.0, 8.1, 8.2, 8.3, 8.4, 8.5, 8.6, 8.7, 8.8, 8.9,  
 9.0, 9.1, 9.2, 9.3, 9.4, 9.5, 9.6, 9.7, 9.8, 9.9, 10.0, 10.1]

**Available metallicities [Fe/H]:**  
[−2.00, −1.75, −1.50, −1.25, −1.00, −0.75, −0.50, −0.25, 0.00, +0.25]

In PARSEC, I set the resolution of the thermal pulse cycles in the COLIBRI section: `ninTPC=20` as detailed in Marigo et al. (2017). This is to make AGB stars more resolved, see [Lee et al. (2025)](https://ui.adsabs.harvard.edu/abs/2025ApJ...995..135L/abstract) for more details.

### Background Sky Model
You also need a realistic sky background (including MW stars and background galaxies). I have also prepared a sky model for you. It takes a simulated galaxy catalog called [JAGUAR](https://fenrir.as.arizona.edu/jaguar/download_jaguar_files.html) and the Milky Way star catalog from the [TRILEGAL](https://stev.oapd.inaf.it/trilegal) model. 

> [!CAUTION]
> JAGUAR is for JWST filters. I need to figure out the conversion from JWST to Roman filters.

The sky model is made following `notebook/Rosesim/01_simulate_JAGUAR_sky.ipynb`, or using the script `rosesim/scripts/sim_sky.py`. If you wanna use other catalogs for background galaxies or MW stars, take a look at the `RomanSky.load_jaguar_bkg` and `RomanSky.load_trilegal_star` functions in `rosesim/rose.py`!

In `ROSESIM_DATA_PATH/sky_jaguar_trilegal`, you will find the sky model I prepared for you. The catalog of galaxies and stars in that sky model is under `ROSESIM_DATA_PATH/sky_jaguar_trilegal/temp/`. You can read the ASDF file and write it to FITS (to view it in DS9) following the example:

```python
import rosesim
sky_dm = rosesim.read_L3_asdf('./F158_642s.asdf')
rosesim.asdf_to_fits(sky_dm, 'F158_642s.fits', subtract_bkg=True)
```


## Usage

Run as a script:

To construct the sky model with JAGUAR galaxies and TRILEGAL stars, use:
```bash
rosesim_sky --obs_ra=150.1049 --obs_dec=2.2741 --size=5001 --prefix='sky_jaguar_trilegal' --exptime=642 --filters="['F106', 'F129', 'F158']" --seed=42 --include_bkg=True --include_star=True --psf_fov_arcsec=10
```
If you only want to make an empty sky model, you can use the following command:
```bash
rosesim_sky --obs_ra=150.1049 --obs_dec=2.2741 --size=5001 --prefix='sky_jaguar_trilegal_new' --exptime=642 --filters="['F106', 'F129', 'F158']" --seed=42 --include_bkg=False --include_star=False
```

To simulate a single galaxy, use:
```bash
rosesim_gal --obs_ra=150.1049 --obs_dec=2.2741 --distance=5 --age=1.0 --log_m_star=4 --exptime=642
```


To read a simulated image and write it to fits:
```python
import rosesim

dm = rosesim.read_L3_asdf('./F158_642s.asdf')
dm.write_fits('./F158_642s.fits')
```


## Requirements
- numpy
- matplotlib
- astropy
- astroquery
- romanisim (https://github.com/AstroJacobLi/romanisim). Note that this version is required because I have modified the code to better support star injection.
- artpop
- asdf
- astrocut
- roman_datamodels

## License
MIT

## Future plans
- Add stellar population information to the `meta` of ASDF
- Enable more sky options
- Enable composite stellar population, enable user inputs on structural parameters
- Improve documentation and examples
