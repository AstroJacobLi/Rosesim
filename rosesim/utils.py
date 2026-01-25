import os, sys
import numpy as np
from astropy import units as u
from astropy.table import Table, Column, vstack, hstack
from astropy.io import fits
import artpop
import asdf
import roman_datamodels as rdm
from . import pixel_scale, DATA_PATH

######### Utility functions for Roman-I-Sim ########

def read_L3_asdf(file):
    af = asdf.open(file)
    dm = rdm.open(af)
    # dm.meta.wcs = wcs.get_mosaic_wcs(dm)
    return dm

# Utility function for creating a WFI WCS to turn catalogue (x, y) into (RA, DEC)
def make_wcs(ra, dec, pa, xy_dim):
    from astropy.wcs import WCS
    wfi_scale = pixel_scale # arcsec/pixel
    w = WCS(naxis=2)
    w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    w.wcs.crpix = xy_dim // 2
    w.wcs.crval = [ra, dec]
    w.wcs.cdelt = [wfi_scale / 3600., wfi_scale / 3600.]
    pc_11 = np.cos(np.radians(pa))
    pc_12 = -np.sin(np.radians(pa))
    pc_21 = np.sin(np.radians(pa))
    pc_22 = np.cos(np.radians(pa))
    w.wcs.pc = [[pc_11, pc_12], [pc_21, pc_22]]
    return w

# Utility function for creating a sample point-source catalogue from
# image co-ordinates and fluxes
def create_point_source_catalogue(wcs, xs, ys, mag_table):
    """
    mag_table should be in maggies. 
    """
    ra_dec = wcs.all_pix2world(xs, ys, 1)
    ras = ra_dec[0]
    decs = ra_dec[1]
    t = Table()
    t['ra'] = Column(data=ras)
    t['dec'] = Column(data=decs)
    t['type'] = Column(data=['PSF']*len(ras))
    t = hstack([t, mag_table])
    return t

def create_smoooth_sersic_catalogue(wcs, src):
    models = [src.smooth_model]
    xs = [model.x_0.value for model in models]
    ys = [model.y_0.value for model in models]

    ra_dec = wcs.all_pix2world(xs, ys, 1)
    ras = ra_dec[0]
    decs = ra_dec[1]

    t = Table()
    t['ra'] = Column(data=ras)
    t['dec'] = Column(data=decs)
    t['type'] = Column(data=['SER']*len(ras))
    t['n'] = Column(data=[model.n.value for model in models])
    t['half_light_radius'] = Column(data=[model.r_eff.value * pixel_scale for model in models])
    t['pa'] = Column(data=[90 - np.rad2deg(model.theta.value) for model in models])
    t['ba'] = Column(data=[1 - model.ellip.value for model in models])
    
    filters = src.mags.colnames
    temp = Table()
    for filt in filters:
        temp["F"+filt[1:]] = [10**(-0.4 * src.sp.mag_integrated_component(filt))] # in maggies

    t = hstack([t, temp])
    return t


def make_jaguar_galaxies(coord,
                         radius=0.1,
                         bandpasses=None,
                         rng=None,
                         seed=50,
                         ):
    """
    Generate a catalog of galaxies around a given coordinate using the JAGUAR mock catalog.
    See https://fenrir.as.arizona.edu/jaguar/jades-mock-catalog_v1.2.pdf and https://fenrir.as.arizona.edu/jaguar/download_jaguar_files.html.
    
    Parameters
    ----------
    coord : astropy.coordinates.SkyCoord
        Location around which to generate sources.
    radius : float
        Radius in degrees in which to uniformly generate sources.
    bandpasses : list[str]
        List of names of bandpasses in which to generate fluxes.
    rng : galsim.BaseDeviate
        Random number generator to use.
    seed : int
        Seed to use for random numbers, only used if rng is None.

    Returns
    -------
    catalog : astropy.Table
        Table for use with table_to_catalog to generate catalog for simulation.
    """
    import galsim

    if rng is None:
        rng = galsim.UniformDeviate(seed)

    import romanisim
    from romanisim import util
    BANDPASSES = set(romanisim.bandpass.galsim2roman_bandpass.values())

    # Generate list of required filters (and Roman, if necessary)
    cos_filt = []
    if bandpasses is None:
        cos_filt = ['HST_F606W_fnu', 'HST_F814W_fnu', 'NRC_F115W_fnu', 'NRC_F150W_fnu', 'NRC_F200W_fnu']
        bandpasses = BANDPASSES
    else:
        for opt_elem in bandpasses:
            if opt_elem == "F062":
                cos_filt.append('HST_F606W_fnu')
            if opt_elem == "F087":
                cos_filt.append('HST_F814W_fnu')
            if opt_elem == "F106": # overlaps with F090W and F115W
                cos_filt.append('NRC_F115W_fnu')
            if opt_elem in ("F129", "F158", "F146"):
                cos_filt.append('NRC_F150W_fnu')
            if opt_elem in ("F213", "F184"):
                cos_filt.append('NRC_F200W_fnu')

    # Open JAGUAR file and pare to required tabs
    cat_all = Table.read(os.path.join(DATA_PATH, "JAGUAR/JADES_all_mock_r1_v1.2.fits"))
    cat_area = 11 * 11 / 3600 # 11x11 arcmin area in square degrees
    
    # remove too faint galaxies
    cat_all = cat_all[-2.5 * np.log10(cat_all['NRC_F115W_fnu'] * 1e-9 / 3631) < 30] # remove sources with F115W flux < 0.1 nJy
    
    cos_density = len(cat_all) / cat_area # number of sources per square degree

    # Calculate total sources
    sim_count = cos_density * np.pi * (radius * u.deg)**2

    # Filter for flags
    cos_filt += ["ID", "RA", "DEC", 'redshift', 'Re_maj', 'axis_ratio', 'sersic_n', 'position_angle']
    cos_filt = list(set(cos_filt))
    
    # Trim catalog
    gal_cat = cat_all[cos_filt]

    # Obtain random sources from the catalog
    rng_numpy_seed = rng.raw()
    rng_numpy = np.random.default_rng(rng_numpy_seed)
    sim_count = rng_numpy.poisson(sim_count.value)
    sim_ids = rng_numpy.integers(size=sim_count, low=0, high=len(gal_cat["ID"])).tolist()
    sim_cat = gal_cat[sim_ids]

    # Match cosmos filters to roman filters
    print('Warning -- need to figure out the coefficients between JWST and Roman filters.')
    for opt_elem in bandpasses:
        if opt_elem == "F062":
            sim_cat['FLUX_F062'] = sim_cat['HST_F606W_fnu'] * 1e-9 / 3631 # convert from nJy to maggies
        elif opt_elem == "F087":
            sim_cat['FLUX_F087'] = sim_cat['HST_F814W_fnu'] * 1e-9 / 3631
        elif opt_elem == "F106":
            F115W = -2.5 * np.log10(sim_cat['NRC_F115W_fnu'] * 1e-9 / 3631)
            F150W = -2.5 * np.log10(sim_cat['NRC_F150W_fnu'] * 1e-9 / 3631)
            Y106 = F115W + (F115W - F150W) * (0.3510699) + 0.00826781 # based on simple SSP, scatter < 0.004 mag
            sim_cat['FLUX_F106'] = 10**(-0.4 * Y106)
        elif opt_elem == "F129":
            F115W = -2.5 * np.log10(sim_cat['NRC_F115W_fnu'] * 1e-9 / 3631)
            F150W = -2.5 * np.log10(sim_cat['NRC_F150W_fnu'] * 1e-9 / 3631)
            J129 = F115W + (F115W - F150W) * (-0.42427082) + -0.05232558 # based on simple SSP, scatter < 0.002 mag
            sim_cat['FLUX_F129'] = 10**(-0.4 * J129)
        elif opt_elem == "F146":
            F115W = -2.5 * np.log10(sim_cat['NRC_F115W_fnu'] * 1e-9 / 3631)
            F150W = -2.5 * np.log10(sim_cat['NRC_F150W_fnu'] * 1e-9 / 3631)
            W146 = F150W + (F115W - F150W) * (0.27329887) + -0.03628754 # based on simple SSP, scatter < 0.004 mag
            sim_cat['FLUX_F146'] = 10**(-0.4 * W146)
        elif opt_elem == "F158":
            F115W = -2.5 * np.log10(sim_cat['NRC_F115W_fnu'] * 1e-9 / 3631)
            F150W = -2.5 * np.log10(sim_cat['NRC_F150W_fnu'] * 1e-9 / 3631)
            H158 = F150W + (F115W - F150W) * (-0.23109312) + -0.11065068 # based on simple SSP, scatter < 0.002 mag
            sim_cat['FLUX_F158'] = 10**(-0.4 * H158)
        elif opt_elem == "F184":
            sim_cat['FLUX_F184'] = sim_cat['NRC_F200W_fnu'] * 1e-9 / 3631
        elif opt_elem == "F213":
            sim_cat['FLUX_F213'] = sim_cat['NRC_F200W_fnu'] * 1e-9 / 3631
        else:
            print(f'Unknown filter {opt_elem} skipped in object catalog creation.')
                
    # Randomize positions of the sources
    locs = util.random_points_in_cap(coord, radius, len(sim_ids), rng=rng)

    # Set profile types
    types = np.zeros(len(sim_ids), dtype='U3')
    types[:] = 'SER'

    # Return Table with source parameters
    out = Table()
    out['ra'] = locs.ra.to(u.deg).value
    out['dec'] = locs.dec.to(u.deg).value
    out['type'] = types

    out['n'] = sim_cat['sersic_n'].astype('f4')
    out['half_light_radius'] = sim_cat['Re_maj'].astype('f4')
    out['pa'] = sim_cat['position_angle'].astype('f4')
    out['ba'] = sim_cat['axis_ratio'].astype('f4')

    # Perturb source fluxes by ~20%
    source_pert = np.ones(len(sim_ids))
    # source_pert += ((0.2) * rng_numpy.normal(size=len(sim_ids)))

    # Convert fluxes to maggies by converting to Jankskys and normalizing for zero-point
    for bandpass in bandpasses:
        # Perturb sources fluxes by 5% per bandwidth
        band_source_pert = 0 # ((0.05) * rng_numpy.normal(size=len(sim_ids)))

        # Convert fluxes to maggies by converting to Jankskys, normalizing for zero-point, and applying perturbations
        out[bandpass] = sim_cat[f'FLUX_{bandpass}'].value * (1 + source_pert + band_source_pert)

    return out

def asdf_to_fits(dm, output_filename, subtract_bkg=True):
    """
    Convert ASDF dataset to FITS file with proper header
    
    Parameters:
    -----------
    dm : ASDF dataset object
        The dataset with .data and .meta attributes
    output_filename : str
        Output FITS filename
    subtract_bkg : bool, optional
        If True, subtract a global background from the data before saving (default is False)
    """

    # Extract the data
    data = dm.data
    
    if subtract_bkg:
        import sep
        bkg = sep.Background(data.astype(float), bw=256, bh=256)
        data -= bkg.globalback
    
    # Create primary HDU with data
    primary_hdu = fits.PrimaryHDU(data=data)
    header = primary_hdu.header
    
    # Add basic observation info
    header['TELESCOP'] = dm.meta.get('telescope', 'UNKNOWN')
    header['INSTRUME'] = dm.meta['instrument']['name']
    header['FILTER'] = dm.meta['instrument']['optical_element']
    header['FILENAME'] = dm.meta.get('filename', 'unknown.fits')
    
    # Add calibration info
    header['CAL_VER'] = dm.meta.get('calibration_software_version', 'UNKNOWN')
    header['CAL_NAME'] = dm.meta.get('calibration_software_name', 'UNKNOWN')
    
    # Add exposure information from coadd_info
    coadd_info = dm.meta.get('coadd_info', {})
    if 'exposure_time' in coadd_info:
        header['EXPTIME'] = coadd_info['exposure_time']
        header['TEXPTIME'] = coadd_info['exposure_time']  # Total exposure time
    
    if 'max_exposure_time' in coadd_info:
        header['MAXEXP'] = coadd_info['max_exposure_time']
    
    # Add timing information
    if 'time_mean' in coadd_info:
        time_mean = coadd_info['time_mean']
        if hasattr(time_mean, 'isot'):
            header['DATE-OBS'] = time_mean.isot
            header['MJD-OBS'] = time_mean.mjd
    
    # Add program information
    program = dm.meta.get('program', {})
    if 'subcategory' in program:
        header['SUBCAT'] = program['subcategory']
    
    
    # Add resample information
    resample = dm.meta.get('resample', {})
    if 'pointings' in resample:
        header['NPOINT'] = resample['pointings']
    if 'pixel_scale_ratio' in resample:
        header['PIXRATIO'] = resample['pixel_scale_ratio']
    if 'pixfrac' in resample:
        header['PIXFRAC'] = resample['pixfrac']
    
    # Add calibration step completion status
    cal_step = dm.meta.get('cal_step', {})
    for step, status in cal_step.items():
        header[f'CAL_{step.upper()[:4]}'] = status
        
    # Add background information if available
    if subtract_bkg:
        header['BKG_SUB'] = True
        header['BKG_TYPE'] = 'GLOBAL'  # Assuming global background subtraction
        header['BKG_MEAN'] = bkg.globalback  # Background level subtracted
    else:
        header['BKG_SUB'] = False
        header['BKG_TYPE'] = 'NONE'
        header['BKG_MEAN'] = 0.0

    # Add WCS information from the WCS object
    if 'wcs' in dm.meta:
        wcs_obj = dm.meta['wcs']
        xy_dim = data.shape
        wcs_header = wcs_obj.to_fits_sip(((0, xy_dim[0]), (0, xy_dim[1])))
        
        # Add WCS keywords to main header
        for key, value in wcs_header.items():
            # FITS header keywords must be <= 8 characters
            if len(key) <= 8:
                try:
                    header[key] = value
                except ValueError:
                    # Skip problematic values
                    continue
    
    # Add additional WCS info from wcsinfo
    wcsinfo = dm.meta.get('wcsinfo', {})
    if wcsinfo:
        # Add basic coordinate info
        if 'ra' in wcsinfo:
            header['RA'] = wcsinfo['ra']
        if 'dec' in wcsinfo:
            header['DEC'] = wcsinfo['dec']
        if 'pixel_scale' in wcsinfo:
            header['PIXSCALE'] = wcsinfo['pixel_scale']
        if 'projection' in wcsinfo:
            header['PROJ'] = wcsinfo['projection']
        if 'orientation' in wcsinfo:
            header['ORIENT'] = wcsinfo['orientation']
        if 's_region' in wcsinfo:
            # S_REGION might be too long for a single header card
            s_region = wcsinfo['s_region']
            if len(s_region) <= 68:  # FITS limit minus keyword name
                header['S_REGION'] = s_region
            else:
                # Split into multiple cards if too long
                header['S_REG1'] = s_region[:68]
                if len(s_region) > 68:
                    header['S_REG2'] = s_region[68:136]
    
    # Add coordinate reference frame
    coordinates = dm.meta.get('coordinates', {})
    if 'reference_frame' in coordinates:
        header['RADESYS'] = coordinates['reference_frame']
    
    # Add file creation date
    if 'file_date' in dm.meta:
        file_date = dm.meta['file_date']
        if hasattr(file_date, 'isot'):
            header['DATE'] = file_date.isot
    
    # Add origin
    header['ORIGIN'] = dm.meta.get('origin', 'UNKNOWN')
    
    # Add model type
    header['MODELTYP'] = dm.meta.get('model_type', 'UNKNOWN')
    
    # Add some comments
    header['COMMENT'] = 'Converted from ASDF format'
    header['COMMENT'] = 'Original file: ' + dm.meta.get('filename', 'unknown')
    header['HISTORY'] = 'Converted using asdf_to_fits function'
    
    # Create HDU list and write to file
    hdul = fits.HDUList([primary_hdu])
    
    # Write to file
    hdul.writeto(output_filename, overwrite=True)
    print(f"Successfully wrote FITS file: {output_filename}")
    
    return None

############### DOLPHOT #####################
import re
from pathlib import Path
import romanisim
import romanisim.bandpass
import romanisim.parameters
from astropy.coordinates import SkyCoord

def _slug(s: str) -> str:
    s = s.strip().lower()
    s = re.sub(r"[^\w]+", "_", s)
    s = re.sub(r"_+", "_", s).strip("_")
    return s

def _parse_imgid(imgid):
    """
    Parse a DOLPHOT ImageID like 'F158_642s' or 'F106_1200sec'.

    Returns
    -------
    filt : str
        Filter name (e.g. 'F158')
    exptime : float or None
        Exposure time in seconds
    """
    m = re.match(r"(?P<filt>[^_]+)(?:_(?P<t>\d+)(?:s|sec))?$", imgid)
    if not m:
        return imgid, None

    filt = m.group("filt")
    t = m.group("t")
    return filt, (float(t) if t is not None else None)

def _parse_dolphot_columns(columns_path, fake=False):
    base = [
        "ext", "chip", "x", "y", "chi", "snr", "sharp", "round",
        "major_axis", "crowd", "obj_type", "pass_det",
    ]

    metric_map = {
        "Measured counts": "counts",
        "Measured sky level": "sky",
        "Normalized count rate": "rate",
        "Normalized count rate uncertainty": "rate_err",
        "Instrumental VEGAMAG magnitude": "mag_vega",
        "Transformed UBVRI magnitude": "mag_ubvri",
        "Magnitude uncertainty": "mag_err",
        "Chi": "chi",
        "Signal-to-noise": "snr",
        "Sharpness": "sharp",
        "Roundness": "round",
        "Crowding": "crowd",
        "Photometry quality flag": "qflag",
    }

    names = []
    filters = set()

    for line in Path(columns_path).read_text().splitlines():
        m = re.match(r"^\s*(\d+)\.\s*(.+?)\s*$", line)
        if not m:
            continue

        idx = int(m.group(1))
        desc = m.group(2)

        if idx <= len(base):
            names.append(base[idx - 1])
            continue

        # "<Metric>, <ImageID> (...)"
        m2 = re.match(r"([^,]+),\s*([^\s(]+)", desc)
        if m2:
            metric, imgid = m2.group(1).strip(), m2.group(2).strip()
            filt, exptime = _parse_imgid(imgid)
            filters.add(filt)
            short = metric_map.get(metric, _slug(metric))
            name = f"{short}_{_slug(filt)}"
        else:
            name = _slug(desc)

        names.append(name)

    filters = sorted(filters, key=lambda f: int(f[1:]))

    # ensure uniqueness
    out, seen = [], {}
    for n in names:
        seen[n] = seen.get(n, 0) + 1
        out.append(n if seen[n] == 1 else f"{n}_{seen[n]}")

    if fake:
        print('Filters:', filters)
        temp = ['ext_ast', 'chip_ast', 'x_true', 'y_true']
        for filt in filters:
            temp += [f'counts_true_{filt.lower()}', f'mag_true_{filt.lower()}']
        out = temp + out
    
    return out, sorted(filters), exptime

def read_dolphot_cat(path, column_path, fake=False):
    import romanisim
    import romanisim.bandpass

    names, filters, exptime = _parse_dolphot_columns(column_path, fake=fake)
    print('Exptime:', exptime, 'Filters:', filters)

    cat = Table.read(path, format='ascii.no_header', names=names)
    cat.remove_columns([f'mag_ubvri_{filt.lower()}' for filt in filters])
    cat.remove_columns([f'rate_{filt.lower()}' for filt in filters])
    cat.remove_columns([f'rate_err_{filt.lower()}' for filt in filters])

    for filt in filters:
        maggytoes = romanisim.bandpass.get_abflux(filt, sca=cat['chip'][0] + 1)
        cat[f'mag_vega_{filt.lower()}'] = -2.5 * np.log10(cat[f'counts_{filt.lower()}'] / exptime / maggytoes)
        if fake:
            cat[f'mag_true_{filt.lower()}'] = -2.5 * np.log10(cat[f'counts_true_{filt.lower()}'] / exptime / maggytoes)

    cat.rename_columns([f'mag_vega_{filt.lower()}' for filt in filters], [f'mag_ab_{filt.lower()}' for filt in filters])

    if fake:
        for col in cat.colnames:
            if cat[col].dtype.kind in "f":  # int or float
                cat[col][cat[col] == 99.999] = np.nan
                cat[col][cat[col] == 9.999] = np.nan
                # cat[col][cat[col] == 0.0] = np.nan
            
    return cat


def xmatch_true_meas(true_cat, cat, radius=0.4*u.arcsec,
                     true_ra='ra', true_dec='dec',
                     meas_ra='ra', meas_dec='dec',
                     true_prefix='true_', meas_prefix='meas_',
                     keep='all'):
    """
    Cross-match measured catalog 'cat' to truth 'true_cat' within 'radius'. This is used to test how much of the sources from Rosesim are detected by Dolphot, for QA purposes.

    keep:
      - 'all'   : keep all matches (many measured can map to same true)
      - 'best'  : keep at most one measured per true (smallest separation)
    """
    c_true = SkyCoord(true_cat[true_ra], true_cat[true_dec], unit='deg')
    c_meas = SkyCoord(cat[meas_ra], cat[meas_dec], unit='deg')

    # For each measured source, find nearest truth source
    idx_true, d2d, _ = c_meas.match_to_catalog_sky(c_true)

    m = d2d < radius
    meas_m = cat[m]
    true_m = true_cat[idx_true[m]]
    sep_m  = d2d[m]

    if keep == 'best':
        # enforce one-to-one on the truth side: keep the closest measured per true
        order = np.argsort(sep_m)  # closest first
        true_ids = idx_true[m][order]
        _, first = np.unique(true_ids, return_index=True)
        keep_idx = order[first]

        meas_m = meas_m[keep_idx]
        true_m = true_m[keep_idx]
        sep_m  = sep_m[keep_idx]

    # rename to avoid column collisions, then merge side-by-side
    meas_m = meas_m.copy()
    true_m = true_m.copy()
    meas_m.rename_columns(meas_m.colnames, [meas_prefix + c for c in meas_m.colnames])
    true_m.rename_columns(true_m.colnames, [true_prefix + c for c in true_m.colnames])

    out = hstack([meas_m, true_m], join_type='exact')
    out['sep_arcsec'] = sep_m.to_value(u.arcsec)

    return out, meas_m, true_m

### DOLPHOT photometric uncertainty and completeness analysis ###
from scipy.optimize import curve_fit
from astropy.stats import sigma_clip

def logistic_completeness(m, m50, w):
    return 1.0 / (1.0 + np.exp((m - m50) / w))

def mag_uncertainty_func(m, slope, intercept):
    return 10**(m * slope + intercept)

def fit_logistic_completeness(mag, comp):
    """
    Fit a logistic completeness function.

    Returns
    -------
    popt : (m50, w)
    pcov : covariance matrix
    """
    mag = np.asarray(mag)
    comp = np.asarray(comp)

    # Use only valid bins
    ok = np.isfinite(mag) & np.isfinite(comp) & (comp > 0) & (comp < 1)
    mag = mag[ok]
    comp = comp[ok]

    # Initial guesses:
    m50_init = np.interp(0.5, comp[::-1], mag[::-1])  # rough 50% estimate
    w_init = 0.3                                      # typical transition width

    popt, pcov = curve_fit(
        logistic_completeness,
        mag, comp,
        p0=[m50_init, w_init],
        bounds=([mag.min(), 0.01], [mag.max(), 5.0]),
    )
    return popt, pcov

def delta_m_scatter_vs_mag(
    mag_in,
    delta_m,
    *,
    mag_bins,
    clip_sigma=3.0,
    clip_iters=3,
    method="mad",
    min_per_bin=10,
):
    """
    Compute sigma-clipped scatter of delta_m as a function of mag_in.

    Parameters
    ----------
    mag_in : array-like
        Input magnitudes.
    delta_m : array-like
        delta_m = mag_in - mag_out (or vice versa; sign does not matter for scatter).
    mag_bins : array-like
        Magnitude bin edges.
    clip_sigma : float
        Sigma threshold for sigma clipping.
    clip_iters : int
        Number of sigma-clipping iterations.
    method : {'mad','rms'}
        Scatter estimator.
    min_per_bin : int
        Minimum number of points required to compute scatter.

    Returns
    -------
    tab : astropy.table.Table
        Columns:
          mag_center, n_used, dmag_med, scatter
    """
    mag_in = np.asarray(mag_in)
    delta_m = np.asarray(delta_m)

    idx = np.digitize(mag_in, mag_bins) - 1

    out = Table()
    out["mag_center"] = 0.5 * (mag_bins[:-1] + mag_bins[1:])
    out["n_used"] = np.zeros(len(out), dtype=int)
    out["dmag_med"] = np.nan
    out["scatter"] = np.nan

    for i in range(len(out)):
        m = (idx == i) & np.isfinite(delta_m)
        if np.sum(m) < min_per_bin:
            continue

        dm = delta_m[m]

        # iterative sigma clipping
        clipped = sigma_clip(
            dm,
            sigma=clip_sigma,
            maxiters=clip_iters,
            cenfunc="median",
            stdfunc="std",
        )

        dm_clipped = dm[~clipped.mask]
        if len(dm_clipped) < min_per_bin:
            continue

        out["n_used"][i] = len(dm_clipped)
        out["dmag_med"][i] = np.median(dm_clipped)

        if method == "mad":
            out["scatter"][i] = (
                1.4826 * np.median(np.abs(dm_clipped - np.median(dm_clipped)))
            )
        elif method == "rms":
            out["scatter"][i] = np.sqrt(np.mean((dm_clipped - np.mean(dm_clipped))**2))
        else:
            raise ValueError("method must be 'mad' or 'rms'")

    return out


def quality_cut(cat, snr=5, crowd={'f158': 0.5, 'f106': 0.5}, sharp={'f158': 0.03, 'f106': 0.03}):
    flag = cat['snr_f158'] > snr
    flag &= cat['snr_f106'] > snr
    flag &= (cat['crowd_f106'] < crowd['f106'])
    flag &= (cat['crowd_f158'] < crowd['f158'])
    flag &= (cat['qflag_f158'] <= 4)
    flag &= (cat['qflag_f106'] <= 4)
    flag &= (cat['obj_type'] <= 2)
    flag &= (cat['pass_det'] < 3)
    return flag

def point_source_cut(cat, snr=5, crowd={'f158': 0.5, 'f106': 0.5}, sharp={'f158': 0.03, 'f106': 0.03}):
    flag = quality_cut(cat, snr=snr, crowd=crowd, sharp=sharp)
    flag &= (cat['sharp_f106']**2 < sharp['f106'])
    flag &= (cat['sharp_f158']**2 < sharp['f158'])
    flag &= (cat['chi_f106'] < 3)
    flag &= (cat['chi_f158'] < 3)
    return flag