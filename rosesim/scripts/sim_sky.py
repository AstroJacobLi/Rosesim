import numpy as np
import os, sys

from astropy.time import Time
from astropy import units as u
import roman_datamodels as rdm
from romanisim import ris_make_utils as ris

os.environ['NUMEXPR_MAX_THREADS'] = "30"
from rosesim.rose import RomanGalaxy, RomanSky
from rosesim import DATA_PATH, pixel_scale
from rosesim.utils import get_subpixel_dither
import fire

def simulate_sky(obs_ra=150.1049, obs_dec=2.2741, size=5001, prefix='sky_jaguar', exptime=642, nexp=6, filters=['F106', 'F129', 'F158'], seed=42, include_bkg=True, include_star=True, exclude_size_thresh=None, trilegal_file="trilegal_Roman_30mag_10h_0deg.dat", fastpointsources=False, psftype='galsim', dither_pattern=None):
    """
    Simulate an mock Roman image for the background sky.

    Parameters
    ----------
    obs_ra : float
        The right ascension of the observation center (degrees).
    obs_dec : float
        The declination of the observation center (degrees).
    size : int
        The size of the simulated image (pixels). The pixel scale of Roman WFI is 0.11 arcsec/pixel.
    prefix : str
        The prefix for the output files.
    exptime : float
        The exposure time for the observations (seconds).
    filters : list
        The list of filters to use for the observations.
    seed : int
        The random seed for reproducibility.
    include_bkg : bool
        Whether to include the background in the simulation.
    include_star : bool
        Whether to include stars in the simulation.
    exclude_size_thresh : float
        Exclude galaxies that are LARGER than this threshold (in arcsec). To only show small objects.
    trilegal_file : str
        The path to the TRILEGAL file. e.g., "trilegal_Roman_30mag_10h_0deg.dat". It needs to be in the ROSESIM_DATA_PATH/TRILEGAL directory.
    psf_fov_arcsec : float
        The field of view of the PSF (point spread function) in arcseconds.
    dither_pattern : str
        The dither pattern to use for the observations, e.g., BOXGAP4_1, SUB4, LINEGAP4_1.
    """
    if dither_pattern is not None:
        do_dither = True
        print(f"Dither following {dither_pattern}")
        pattern = get_subpixel_dither(obs_ra, obs_dec, dither_pattern, subpix=True)
    else:
        do_dither = False
        pattern = [[0.0, 0.0]]

    if not do_dither:
        sky = RomanSky(obs_ra, obs_dec, xy_dim=np.array([size, size]), prefix=prefix, suffix=None, data_dir=DATA_PATH)
        radius = sky.xy_dim[0] * 0.11 / 3600
        sky.obs_time = Time('2025-01-01T00:00:00')
        sky.load_trilegal_star(radius=radius, area=2, path=os.path.join(DATA_PATH, 'TRILEGAL', trilegal_file))
        sky.load_jaguar_bkg(radius=radius, seed=seed)
        sky.gen_catalog(include_bkg=include_bkg, include_star=include_star, exclude_size_thresh=exclude_size_thresh)
        for filt in filters:
            sky.observe(filt, exptime=exptime, nexp=nexp, fastpointsources=fastpointsources, psftype=psftype)
    else:
        for i, dither in enumerate(pattern):
            sky = RomanSky(obs_ra, obs_dec, dither=dither, xy_dim=np.array([size, size]), prefix=prefix, suffix=f"{dither_pattern}.{i}", data_dir=DATA_PATH)
            radius = sky.xy_dim[0] * 0.11 / 3600
            sky.obs_time = Time('2025-01-01T00:00:00')
            sky.load_trilegal_star(radius=radius, area=2, path=os.path.join(DATA_PATH, 'TRILEGAL', trilegal_file))
            sky.load_jaguar_bkg(radius=radius, seed=seed)
            sky.gen_catalog(include_bkg=include_bkg, include_star=include_star, exclude_size_thresh=exclude_size_thresh)
            for filt in filters:
                sky.observe(filt, exptime=exptime, nexp=nexp, fastpointsources=fastpointsources, psftype=psftype)

def main():
    fire.Fire(simulate_sky)

if __name__ == '__main__':
    main()

# rosesim_sky --obs_ra=150.1049 --obs_dec=2.2741 --size=5001 --prefix='sky_jaguar_trilegal_new' --exptime=642 --nexp=6 --filters="['F106', 'F129', 'F158']" --seed=42 --include_bkg=True --include_star=True --trilegal_file="trilegal_Roman_30mag_10h_0deg.dat" --psf_fov_arcsec=10

# If you wanna make an empty sky, you can use the following command:
# rosesim_sky --obs_ra=150.1049 --obs_dec=2.2741 --size=5001 --prefix='sky_jaguar_trilegal_new' --exptime=642 --nexp=6 --filters="['F106', 'F129', 'F158']" --seed=42 --include_bkg=False --include_star=False --trilegal_file="trilegal_Roman_30mag_10h_0deg.dat" --psf_fov_arcsec=10

# python sim_sky.py --obs_ra=201.3704160 --obs_dec=-43.0166667 --size=501 --prefix='sky_jaguar_trilegal_cena_test'   --exptime=600 --nexp=1 --filters="['F106']"   --seed=42 --include_bkg=False --include_star=True   --fastpointsources=False --psftype='epsf' --trilegal_file="trilegal_Roman_30mag_2deg2_CenA.dat" --dither_pattern="LINEGAP5_1"