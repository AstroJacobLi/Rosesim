import numpy as np
from astropy.wcs import WCS
from astropy.table import Table, Column, vstack, hstack
from astropy.io import fits
from astropy.time import Time
from astropy.coordinates import SkyCoord

pixel_scale = 0.11 # arcsec/pixel
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

def gen_catalog_from_artpop(src, ra=180., dec=0., pa=0.):
    # check if src is loaded
    roman_wcs = make_wcs(ra, dec, pa, src.xy_dim)
    
    # In Romanisim, the catalog brightnesses are specified in units of maggies, 
    # which are defined such that one maggie is equal to the reference AB magnitude flux (3,631 Jy), 
    # i.e., maggies = 10^(-0.4 * m_AB).
    filters = src.mags.colnames
    mag_table = src.mags.copy()
    mag_table.rename_columns(mag_table.colnames, ['F' + name[1:] for name in mag_table.colnames])
    for filt in mag_table.colnames:
        mag_table[filt] = 10**(-0.4 * mag_table[filt])
    pts_table = create_point_source_catalogue(roman_wcs, src.x, src.y, mag_table)
    if len(pts_table) == 0:
        print("No point sources found in the source catalog.")
    else:
        pass
        # pts_table.write(f'./temp/pts_table.ecsv', format='ascii.ecsv', overwrite=True)

    gal_table = create_smoooth_sersic_catalogue(roman_wcs, src)
    # gal_table.write(f'./temp/gal_table.ecsv', format='ascii.ecsv', overwrite=True)
    
    full_table = gal_table.copy()
    if len(pts_table) != 0:
        full_table = vstack([full_table, pts_table])

    # full_table.write(f'./temp/full_table.ecsv', format='ascii.ecsv', overwrite=True)
    return full_table


############### DOLPHOT #####################
import re
from pathlib import Path
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

"""
Utility functions for photometry
"""
def flux_to_mag(flux):
    """
    Convert flux from maggie to AB magnitude.

    AB_mag = -2.5 * log10(flux [maggie])
    """
    flux = np.asarray(flux)
    m = np.full_like(flux, np.nan, dtype=float)
    good = flux > 0
    m[good] = -2.5 * np.log10(flux[good])
    return m

def mag_to_flux(mag):
    """Convert AB magnitude to flux in maggie.

    flux [maggie] = 10**(-0.4 * AB_mag)
    """
    mag = np.asarray(mag)
    return 10**(-0.4 * mag)

def logistic_completeness(m, m50, w):
    """
    Logistic completeness function.

    completeness = 1 / (1 + exp((m - m50) / w))
    """
    return 1.0 / (1.0 + np.exp((m - m50) / w))

def mag_uncertainty_func(m, slope, intercept):
    """
    Linear uncertainty function, used to model photometric scatter.

    uncertainty = 10**(m * slope + intercept)
    """
    return 10**(m * slope + intercept)

def apply_obs_model_two_band(
    mag1,
    mag2,
    completeness_dict,
    mag_uncertainty_dict,
    bands=['F106', 'F158'],
    depth_HLWAS=1,
    rng=None,
    detection="both",
    ref_band="F158",
    clip_sigma=None,
):
    """
    Forward-model incompleteness + photometric scatter for two-band photometry.

    Parameters
    ----------
    mag1, mag2 : array-like
        "True" magnitudes (must be >0 to yield finite magnitudes).
    completeness_dict : dict
        e.g. {'F106': [m50, w], 'F158': [m50, w]}
    mag_uncertainty_dict : dict
        e.g. {'F106': [slope, intercept], 'F158': [slope, intercept]}
    rng : np.random.Generator, optional
        If None, uses default_rng().
    detection : {"both","either","ref"}
        - "both": require detection in both bands (common for colors).
        - "either": detected if either band is detected.
        - "ref": detected if detected in ref_band only (common for detection image).
    ref_band : {"F106","F158"}
        Used only if detection="ref".
    clip_sigma : float or None
        If set, clip magnitude perturbations to +/- clip_sigma * sigma_m (rarely needed).

    Returns
    -------
    data : structured array
        Keys include:
          - 'detected' (bool mask)
          - 'm_true_1', 'm_true_2'
          - 'm_obs_1', 'm_obs_2' (NaN if not detected)
          - 'sigma_m_1', 'sigma_m_2' (evaluated at m_true)
          - 'pdet_1', 'pdet_2'
    """
    if rng is None:
        rng = np.random.default_rng()

    mag1 = np.asarray(mag1, dtype=float)
    mag2 = np.asarray(mag2, dtype=float)

    # True magnitudes
    m_true_1 = mag1
    m_true_2 = mag2

    delta_mag = 1.25 * np.log10(depth_HLWAS)
    # Detection probabilities from completeness curves (evaluated at true mags)
    pdet_1 = logistic_completeness(m_true_1 - delta_mag, *completeness_dict[bands[0]])
    pdet_2 = logistic_completeness(m_true_2 - delta_mag, *completeness_dict[bands[1]])

    # Realize detections (Bernoulli trials)
    det_1 = rng.random(size=m_true_1.size) < pdet_1
    det_2 = rng.random(size=m_true_2.size) < pdet_2

    if detection == "both":
        detected = det_1 & det_2
    elif detection == "either":
        detected = det_1 | det_2
    elif detection == "ref":
        if ref_band == "F106":
            detected = det_1
        elif ref_band == "F158":
            detected = det_2
        else:
            raise ValueError("ref_band must be 'F106' or 'F158'")
    else:
        raise ValueError("detection must be one of {'both','either','ref'}")

    # Photometric uncertainty model (evaluated at true mags)
    sigma_m_1 = mag_uncertainty_func(m_true_1 - delta_mag, *mag_uncertainty_dict[bands[0]])
    sigma_m_2 = mag_uncertainty_func(m_true_2 - delta_mag, *mag_uncertainty_dict[bands[1]])

    # Initialize observed mags as NaN (non-detections remain missing)
    m_obs_1 = np.full_like(m_true_1, np.nan, dtype=float)
    m_obs_2 = np.full_like(m_true_2, np.nan, dtype=float)

    # Draw perturbations for detected objects only
    idx = np.where(detected)[0]
    if idx.size > 0:
        dm1 = rng.normal(loc=0.0, scale=sigma_m_1[idx])
        dm2 = rng.normal(loc=0.0, scale=sigma_m_2[idx])

        if clip_sigma is not None:
            dm1 = np.clip(dm1, -clip_sigma * sigma_m_1[idx], clip_sigma * sigma_m_1[idx])
            dm2 = np.clip(dm2, -clip_sigma * sigma_m_2[idx], clip_sigma * sigma_m_2[idx])

        m_obs_1[idx] = m_true_1[idx] + dm1
        m_obs_2[idx] = m_true_2[idx] + dm2

    # Create structured array
    data = np.zeros(len(m_true_1), dtype=[
        ('RA', 'f8'), ('DEC', 'f8'),
        (f'MAG_{bands[0]}', 'f4'), (f'MAG_{bands[1]}', 'f4'),
        (f'MAG_ERR_{bands[0]}', 'f4'), (f'MAG_ERR_{bands[1]}', 'f4'),
        (f'MAG_DERED_{bands[0]}', 'f4'), (f'MAG_DERED_{bands[1]}', 'f4'),
        ('MC_SOURCE_ID', 'i8'), ('sharp', 'f4'), ('crowd', 'f4'), ('snr', 'f4'), ('quality_flag', bool), ('star_flag', bool)
    ])
    data[f'MAG_{bands[0]}'] = m_obs_1
    data[f'MAG_{bands[1]}'] = m_obs_2
    data[f'MAG_ERR_{bands[0]}'] = sigma_m_1
    data[f'MAG_ERR_{bands[1]}'] = sigma_m_2
    data[f'MAG_DERED_{bands[0]}'] = m_obs_1
    data[f'MAG_DERED_{bands[1]}'] = m_obs_2
    data['MC_SOURCE_ID'] = np.arange(len(m_true_1))
    data['sharp'] = np.full_like(m_true_1, 0, dtype=float)
    data['crowd'] = np.full_like(m_true_1, 0, dtype=float)
    data['snr'] = np.full_like(m_true_1, 100, dtype=float)
    data['quality_flag'] = np.full_like(m_true_1, True, dtype=bool)
    data['star_flag'] = np.full_like(m_true_1, True, dtype=bool)

    return data

def setup_dolphot_cat(cat, bands=['F106', 'F158']):
    # Create structured array
    data = np.zeros(len(cat), dtype=[
        ('RA', 'f8'), ('DEC', 'f8'),
        (f'MAG_{bands[0]}', 'f4'), (f'MAG_{bands[1]}', 'f4'),
        (f'MAG_ERR_{bands[0]}', 'f4'), (f'MAG_ERR_{bands[1]}', 'f4'),
        (f'MAG_DERED_{bands[0]}', 'f4'), (f'MAG_DERED_{bands[1]}', 'f4'),
        ('MC_SOURCE_ID', 'i8'), ('sharp', 'f4'), ('crowd', 'f4'), ('snr', 'f4'), ('quality_flag', bool), ('star_flag', bool)
    ])

    data['RA'] = cat['ra']
    data['DEC'] = cat['dec']
    data[f'MAG_{bands[0]}'] = cat['mag_ab_f106']
    data[f'MAG_{bands[1]}'] = cat['mag_ab_f158']
    data[f'MAG_ERR_{bands[0]}'] = cat['mag_err_f106']
    data[f'MAG_ERR_{bands[1]}'] = cat['mag_err_f158']
    data[f'MAG_DERED_{bands[0]}'] = cat['mag_ab_f106']
    data[f'MAG_DERED_{bands[1]}'] = cat['mag_ab_f158']
    data['sharp'] = cat['sharp']
    data['crowd'] = cat['crowd']
    data['snr'] = cat['snr']
    data['quality_flag'] = cat['quality_flag']
    data['star_flag'] = cat['star_flag']
    
    return data