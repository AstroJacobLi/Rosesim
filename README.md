# RoSE-Sim: Roman Semi-resolved Galaxy Simulator

**RoSE-Sim** is a Python package designed for creating image simulations of semi-resolved dwarf galaxies for the Nancy Grace Roman Space Telescope.

<!-- insert demo.png  -->
![Concept of Rosesim](demo.png)

## Installation

### 1. Install the Package
You can install **RoSE-Sim** locally by cloning the repository and running:

```bash
git clone https://github.com/YourUsername/Rosesim.git
cd Rosesim
pip install -e .
```

### 2. Set Up Environment Variables
**RoSE-Sim** requires a dedicated data directory to store large files (e.g., isochrones, sky models). 

1. Create a directory for the data (e.g., `/scratch/gpfs/user/Rosesim_Data`).
2. Set the `ROSESIM_DATA_PATH` environment variable pointing to this directory. Add the following line to your shell configuration file (e.g., `~/.bashrc` or `~/.zshrc`):

    ```bash
    export ROSESIM_DATA_PATH="/path/to/your/data/directory"
    ```

3. Reload your shell configuration:
    ```bash
    source ~/.bashrc
    ```

### 3. Download Required Data
Once the package is installed and the environment variable is set, you can easily download the necessary data files using the built-in fetch function:

```python
import rosesim
rosesim.fetch_data()
```
This command downloads the required isochrones and sky models to your `ROSESIM_DATA_PATH`.

## Data Description

### Stellar Population Synthesis Models
**RoSE-Sim** uses [PARSEC isochrones](https://stev.oapd.inaf.it/cgi-bin/cmd) for stellar population synthesis. 

- **Isochrones:** Pre-packaged for Roman filters. Note that Roman isochrones are provided in **Vega** magnitudes.
- **Conversion:** Zeropoint offsets to convert to **AB** magnitudes are available in `rosesim.Roman_zp_AB_Vega_mist`.
- **AGB Stars:** Special attention has been given to the resolution of thermal pulse cycles (`ninTPC=20` in COLIBRI tracks) to strictly resolve AGB stars (see [Lee et al. 2025](https://ui.adsabs.harvard.edu/abs/2025ApJ...995..135L/abstract)).

**Available Grid:** (unfortunately, `artpop` doesn't support the PARSEC isochrones to be interpolated yet)
- **Filters:** F062, F087, F106, F129, F146, F158, F184, F213
- **Log(Age/yr):** [8.0, 8.1, 8.2, 8.3, 8.4, 8.5, 8.6, 8.7, 8.8, 8.9, 9.0, 9.1, 9.2, 9.3, 9.4, 9.5, 9.6, 9.7, 9.8, 9.9, 10.0, 10.1]
- **Metallicity [Fe/H]:** [-2.00, -1.75, -1.50, -1.25, -1.00, -0.75, -0.50, -0.25, 0.00, +0.25]

### Milky Way Star Model
A Milky Way star catalog is generated using the [TRILEGAL](https://stev.oapd.inaf.it/trilegal) model.
- **Default Field:** RA = 10h, Dec = 2.27 deg
- **Depth:** Limiting magnitude of 30 mag in the 4th filter.
- **Location:** Stored in `ROSESIM_DATA_PATH/TRILEGAL/`.

### Background Sky Model
To simulate realistic observations, `Rosesim` includes a background sky model comprising:
1. **Background Galaxies:** From the [JAGUAR](https://fenrir.as.arizona.edu/jaguar/download_jaguar_files.html) mock catalog.
2. **Foreground Stars:** From the TRILEGAL model.

> [!CAUTION]
> JAGUAR fluxes are currently based on JWST filters. A conversion to Roman filters is pending.

## Usage

### Command Line Interface

**RoSE-Sim** provides command-line scripts for common simulation tasks.

#### 1. Simulate a Sky Model
Generate a full sky model including background galaxies and Milky Way stars:
```bash
rosesim_sky \
  --obs_ra=150.1049 \
  --obs_dec=2.2741 \
  --size=5001 \
  --prefix='sky_jaguar_trilegal' \
  --exptime=642 \
  --filters="['F106', 'F129', 'F158']" \
  --seed=42 \
  --include_bkg=True \
  --include_star=True \
  --psf_fov_arcsec=10
```

To generate an **empty sky** (for noise-only or background-free simulations):
```bash
rosesim_sky \
  --obs_ra=150.1049 \
  --obs_dec=2.2741 \
  --size=5001 \
  --prefix='empty_sky' \
  --exptime=642 \
  --filters="['F106', 'F129', 'F158']" \
  --seed=42 \
  --include_bkg=False \
  --include_star=False
```

#### 2. Simulate a Single Dwarf Galaxy
Inject a specific dwarf galaxy into a simulation:
```bash
rosesim_gal \
  --obs_ra=150.1049 \
  --obs_dec=2.2741 \
  --distance=5 \
  --age=1.0 \
  --log_m_star=4 \
  --exptime=642
```

### Python API

You can also use the Python API to inspect or manipulate data.

**Example: inspect the sky model**
```python
import rosesim

# Read the simulated ASDF file
sky_dm = rosesim.read_L3_asdf('./F158_642s.asdf')

# Convert to FITS for inspection (e.g., with DS9)
rosesim.asdf_to_fits(sky_dm, 'F158_642s.fits', subtract_bkg=True)
```

## Requirements
The package relies on the following libraries:
- `numpy`, `matplotlib`, `astropy`, `astroquery`
- `artpop`, `asdf`, `astrocut`, `roman_datamodels`
- `romanisim` (Modified version required: [GitHub](https://github.com/AstroJacobLi/romanisim))

## Future Plans
- [ ] Add stellar population info to ASDF `meta`.
- [ ] Support more diverse sky background options.
- [ ] Support composite stellar populations and user-defined structural parameters.
- [ ] Expand documentation and examples.

## License
MIT
