import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, hstack
from astropy.io import fits
from astropy.wcs import WCS
import astropy.units as u

from astropy.stats import sigma_clipped_stats
from astropy.modeling.fitting import LevMarLSQFitter
from photutils.psf import EPSFModel, ImagePSF
from photutils.detection import DAOStarFinder
from photutils.datasets import make_model_image
from photutils.psf import PSFPhotometry, SourceGrouper
from photutils.aperture import CircularAperture, aperture_photometry

# For parallelization
from joblib import Parallel, delayed

import fire
from kuaizi.display import draw_circles
from astropy.convolution import convolve_fft

# Project paths - adjust as needed
# sys.path.append('/home/jiaxuanl/Research/SALAD/script/')
# os.chdir('/scratch/gpfs/JENNYG/jiaxuanl/Data/SBF/Rosesim/')


def inject_stars(bkg_image, epsf_model, stars):
    """
    Inject stars into a background image.

    Parameters
    ----------
    bkg_image : (ny, nx) ndarray
    epsf_model : EPSFModel (Astropy model)
    stars : astropy.table.Table with columns: x_0, y_0, flux

    Returns
    -------
    img_inj : ndarray
    """
    model_img = make_model_image(
        shape=bkg_image.shape,
        model=epsf_model,
        params_table=stars
    )
    return bkg_image + model_img


def run_photometry(img, epsf_model, fwhm=2.0, threshold_sigma=5.0, error=None, aperture_radius=2.0,
    aperture_positions="psf"):
    """
    Detect sources, perform PSF photometry, and optionally compute aperture fluxes.

    Parameters
    ----------
    img : 2D ndarray
        Image to photometer.
    epsf_model : astropy.modeling.Model
        PSF model (e.g., EPSFModel or ImagePSF) with a defined bounding box.
    fwhm : float
        FWHM (pixels) for DAOStarFinder.
    threshold_sigma : float
        Detection threshold in units of the background std.
    error : float or 2D ndarray, optional
        1-sigma per-pixel uncertainties. If provided, will be passed to PSFPhotometry and
        used to compute aperture flux errors (quadrature sum).
    aperture_radius : float or None
        Radius (pixels) for circular aperture photometry. If None, skip aperture photometry.
    aperture_positions : {"psf","det"}
        Use PSF-fitted (x_fit,y_fit) positions or detection (xcentroid,ycentroid) positions
        for aperture photometry.

    Returns
    -------
    res : astropy.table.Table
        Table including PSF-fit results plus:
            - ap_flux (if aperture_radius is not None)
            - ap_flux_err (if error is not None and aperture_radius is not None)
    """
    mean, med, std = sigma_clipped_stats(img, sigma=4.0)

    finder = DAOStarFinder(fwhm=fwhm, threshold=threshold_sigma * std)
    det = finder(img - med)
    if det is None or len(det) == 0:
        return Table()

    # Photutils expects initial guesses named x_0, y_0 in many PSF classes
    init = Table()
    init["x_0"] = det["xcentroid"]
    init["y_0"] = det["ycentroid"]
    init["flux"] = np.maximum(det["flux"], 1.0)

    fitter = LevMarLSQFitter()
    grouper = SourceGrouper(min_separation=2.0)
    phot = PSFPhotometry(
        psf_model=epsf_model,
        fit_shape=(15, 15),
        finder=None,          # we already ran a finder
        grouper=grouper,
        fitter=fitter,
        localbkg_estimator=None,
    )
    res = phot(img, init_params=init, error=error)

    # -----------------------------------------
    # Add aperture photometry (optional)
    # -----------------------------------------
    if aperture_radius is not None:
        # choose positions for apertures
        if aperture_positions == "psf":
            # be robust to different column names across versions
            if "x_fit" in res.colnames and "y_fit" in res.colnames:
                x_ap, y_ap = res["x_fit"], res["y_fit"]
            else:
                # fall back (some versions keep x_0/y_0)
                x_ap, y_ap = res["x_0"], res["y_0"]
        elif aperture_positions == "det":
            x_ap, y_ap = det["xcentroid"], det["ycentroid"]
        else:
            raise ValueError("aperture_positions must be 'psf' or 'det'")

        positions = np.transpose([np.asarray(x_ap), np.asarray(y_ap)])
        aper = CircularAperture(positions, r=aperture_radius)

        # aperture photometry on background-subtracted image
        ap = aperture_photometry(img - med, aper, error=error)

        # aperture_photometry returns 'aperture_sum' and optionally 'aperture_sum_err'
        res["ap_flux"] = np.asarray(ap["aperture_sum"])
        if error is not None and "aperture_sum_err" in ap.colnames:
            res["ap_flux_err"] = np.asarray(ap["aperture_sum_err"])

    return res


def match_injected_to_recovered(stars_in, stars_out, rmatch=1.0):
    """
    Simple nearest-neighbor match within rmatch (pixels).

    Returns: recovered_mask (len = N_in), plus matched output indices (or -1).
    """
    if stars_out is None or len(stars_out) == 0:
        return np.zeros(len(stars_in), dtype=bool), -np.ones(len(stars_in), dtype=int)

    # Try to find plausible output column names across versions
    xout_candidates = ["x_fit", "x_0", "x"]
    yout_candidates = ["y_fit", "y_0", "y"]
    
    def first_col(tbl, names):
        for n in names:
            if n in tbl.colnames:
                return np.asarray(tbl[n])
        raise KeyError(f"None of {names} found in output table columns {tbl.colnames}")

    x_out = first_col(stars_out, xout_candidates)
    y_out = first_col(stars_out, yout_candidates)

    x_in = np.asarray(stars_in["x_0"])
    y_in = np.asarray(stars_in["y_0"])

    recovered = np.zeros(len(stars_in), dtype=bool)
    match_idx = -np.ones(len(stars_in), dtype=int)

    for i in range(len(stars_in)):
        dx = x_out - x_in[i]
        dy = y_out - y_in[i]
        r = np.hypot(dx, dy)
        j = np.argmin(r)
        if r[j] <= rmatch:
            recovered[i] = True
            match_idx[i] = j

    return recovered, match_idx


def one_star_trial(bkg_imgs, epsfs, *, x0, y0, fluxes,
                   fwhm=2.0, threshold_sigma=5.0, rmatch=1.0, aperture_radius=2.0,
                   bands=['F106', 'F158'], display=False):
    """
    Inject one star into both bands and test if it is recovered in EACH band independently.
    
    Parameters
    ----------
    bkg_imgs : dict
        Dictionary of background images, keys=bands.
    epsfs : dict
        Dictionary of EPSF models, keys=bands.
    fluxes : dict
        Dictionary of fluxes, keys=bands.
    """
    
    result_row = {'x_in': x0, 'y_in': y0}

    for band in bands:
        bkg_img = bkg_imgs[band]
        epsf = epsfs[band]
        flux = fluxes[band]

        # truth table for one star
        stars_in = Table()
        stars_in["x_0"] = [float(x0)]
        stars_in["y_0"] = [float(y0)]
        stars_in["flux"] = [float(flux)]

        img_inj = inject_stars(bkg_img, epsf, stars_in)
        
        # Estimate error map
        bkg_std = sigma_clipped_stats(bkg_img, sigma=4.0)[2]
        error_map = bkg_std * np.ones_like(img_inj)

        stars_out = run_photometry(img_inj, epsf, fwhm=fwhm, threshold_sigma=threshold_sigma, 
                                   error=error_map, aperture_radius=aperture_radius)

        if stars_out is None or len(stars_out) == 0:
            rec, sep, ffit, ffit_err = False, np.nan, np.nan, np.nan
        else:
            recovered_mask, match_idx_arr = match_injected_to_recovered(stars_in, stars_out, rmatch=rmatch)
            rec = recovered_mask[0]
            
            if rec:
                j = match_idx_arr[0]
                
                # Robustly get columns
                def getcol(tbl, candidates):
                    for n in candidates:
                        if n in tbl.colnames:
                            return tbl[n][j]
                    return np.nan
                
                xfit = getcol(stars_out, ["x_fit", "x_0", "x"])
                yfit = getcol(stars_out, ["y_fit", "y_0", "y"])
                ffit = getcol(stars_out, ["flux_fit", "flux_0", "flux"])
                ffit_err = getcol(stars_out, ["flux_err"])

                sep = np.hypot(xfit - x0, yfit - y0)
            else:
                sep, ffit, ffit_err = np.nan, np.nan, np.nan

        # Store results for this band
        result_row[f'recovered_{band}'] = rec
        result_row[f'sep_{band}'] = sep
        result_row[f'flux_out_{band}'] = ffit
        result_row[f'flux_err_{band}'] = ffit_err
        result_row[f'flux_in_{band}'] = flux
        # SNR
        result_row[f'snr_{band}'] = ffit / ffit_err

        if display:
            if len(stars_out) == 0:
                draw_circles(convolve_fft(img_inj, epsf.data), stars_in, colnames=['x_0', 'y_0'], color='cyan')
            else:
                draw_circles(convolve_fft(img_inj, epsf.data), stars_out, colnames=['x_fit', 'y_fit'])

    # return result_row
    return stars_out


def run_ast_parallel(bkg_imgs, epsfs, ZPs, 
                     mag_range_f158, color_range,
                     n_trials=200, margin=10, 
                     fwhm=2.0, threshold_sigma=5.0, rmatch=1.0, 
                     bands=['F106', 'F158'], n_jobs=-1, seed=None):
    
    if seed is None:
        rng = np.random.default_rng()
    else:
        rng = np.random.default_rng(seed)

    ny, nx = bkg_imgs[bands[0]].shape

    # Pre-generate parameters
    # 1. Sample F158 magnitudes uniformly
    mags_f158 = rng.uniform(mag_range_f158[0], mag_range_f158[1], n_trials)
    
    # 2. Sample Color (F106 - F158) uniformly
    colors = rng.uniform(color_range[0], color_range[1], n_trials)
    
    # 3. Calculate F106 magnitudes
    mags_f106 = mags_f158 + colors
    
    # 4. Positions
    x0s = rng.uniform(margin, nx - margin, n_trials)
    y0s = rng.uniform(margin, ny - margin, n_trials)

    # 5. Calculate Fluxes
    fluxes_list = []
    for i in range(n_trials):
        f_dict = {}
        # F158
        f_dict['F158'] = 10 ** (-0.4 * (mags_f158[i] - ZPs['F158']))
        # F106
        f_dict['F106'] = 10 ** (-0.4 * (mags_f106[i] - ZPs['F106']))
        fluxes_list.append(f_dict)

    # Parallel Execution
    results = Parallel(n_jobs=n_jobs)(
        delayed(one_star_trial)(
            bkg_imgs, epsfs, x0=x0s[i], y0=y0s[i], fluxes=fluxes_list[i],
            fwhm=fwhm, threshold_sigma=threshold_sigma, rmatch=rmatch, bands=bands, 
            aperture_radius=2.0, aperture_positions="psf"
        ) for i in range(n_trials)
    )
    
    # Convert list of dicts to Table
    t = Table(results)
    
    # Add input magnitudes for convenience
    t['mag_in_F158'] = mags_f158
    t['mag_in_F106'] = mags_f106
    
    # Calculate output magnitudes where recovered
    for band in bands:
        zp = ZPs[band]
        f_out = t[f'flux_out_{band}']
        # Avoid log(<=0) or log(NaN)
        valid = (t[f'recovered_{band}']) & (f_out > 0)
        
        mag_out = np.full(len(t), np.nan)
        mag_out[valid] = -2.5 * np.log10(f_out[valid]) + zp
        t[f'mag_out_{band}'] = mag_out
        t[f'dmag_{band}'] = t[f'mag_out_{band}'] - t[f'mag_in_{band}']

    return t


def run(n_jobs=25, n_trials=1000000):
    # Define bands and Zero Points
    bands = ['F106', 'F158']
    zp = -2.5 * np.log10(1 * ((0.11 * u.arcsec)**2).to(u.steradian).value * 1e6 / 3631)
    ZPs = {'F106': zp, 'F158': zp}

    # PSF setup
    import romanisim
    import romanisim.parameters
    from romanisim.l3 import l3_psf
    scale = romanisim.parameters.pixel_scale
    nx_psf = ny_psf = (int(5 / scale) // 2) * 2 + 1

    epsfs = {}
    psf_imgs = {}

    # Generate PSFs for both bands
    for band in bands:
        # Use appropriate wavelength for PSF generation if available, else approximate
        # romanisim might handle string band names
        roman_psf = l3_psf(band, scale=1, stpsf=True, fov_arcsec=5) # scale might need adjustment per band wavelength if critical
        img_psf_model = roman_psf.drawImage(nx=nx_psf, ny=ny_psf, scale=scale, method="auto")
        psf_data = img_psf_model.array / img_psf_model.array.sum()
        
        psf_imgs[band] = psf_data
        epsfs[band] = ImagePSF(psf_data, oversampling=1)


    # Load Background Images
    bkg_imgs = {}

    import sep

    for band in bands:
        fname = f'/scratch/gpfs/JENNYG/jiaxuanl/Data/SBF/Rosesim/sky_jaguar_test/{band}_642s.fits'
        hdu = fits.open(fname)
        # w = WCS(hdu[0].header) # Not stricly used for injection unless RA/Dec needed
        bkg_data = hdu[0].data.astype(float)
        
        # Subtract background
        bkg_model = sep.Background(bkg_data, bh=128, bw=128)
        bkg_data -= bkg_model.globalback
        bkg_imgs[band] = bkg_data
        print(f"Loaded {band} image from {fname}")

    # Run AST
    # Define parameters
    mag_range_f158 = (24.0, 27.5)
    color_range = (-0.5, 0.8) # F106 - F158

    print(f"Running AST trials...")
    print(f"F158 Range: {mag_range_f158}")
    print(f"Color (F106-F158) Range: {color_range}")

    n_trials = 1000000
    tab = run_ast_parallel(
        bkg_imgs, epsfs, ZPs,
        mag_range_f158=mag_range_f158,
        color_range=color_range,
        n_trials=n_trials,      # Number of stars to inject
        n_jobs=25,          # Adjust based on available cores
        margin=30,
        fwhm=1.25,
        threshold_sigma=3.5,
        rmatch=2.0,
        bands=bands
    )

    # Show results sample
    print("\nResults Sample:")
    print(tab[:5])

    # Save results?
    tab.write(f'/scratch/gpfs/JENNYG/jiaxuanl/Data/SBF/Rosesim/sky_jaguar_test/ast_results_{n_trials}.fits', overwrite=True)

# use fire to make it callable
if __name__ == '__main__':
    fire.Fire(run)

# example call
# python ast_simulation.py run --n_jobs=25 --n_trials=1000000