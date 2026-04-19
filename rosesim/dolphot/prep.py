"""
Prepare the FITS files from Rosesim for running DOLPHOT.

This is inspired by the romanmask.c code from DOLPHOT.

Jiaxuan Li, 2025-12-24
"""
import fire
import sys
import numpy as np
from astropy.io import fits
import os
import glob

# filters
filters = ['F062', 'F087', 'F106', 'F129', 'F146', 'F158', 'F184', 'F213']

import romanisim
import romanisim.parameters
import romanisim.bandpass

def prep(*filenames, exptime=None, rdnoise=None, chip=1, nodark=False):
    """
    Prepare Roman FITS files for DOLPHOT, bypassing standard romanmask checks.

    Args:
        filenames: Input FITS file(s).
        exptime: Override exposure time.
        rdnoise: Override read noise.
        chip: Chip number (1-18) to assign (default: 1 -> sca01).
        nodark: Ignore dark current in noise calculation.
    """
    if not filenames:
        print("No filenames provided.")
        return

    # can we get the files from a wildcard? 
    # e.g., the first argument is a wildcard, and we expand it to a list of files
    if len(filenames) == 1 and "*" in filenames[0]:
        filenames = glob.glob(filenames[0])

    # calculate etomjysr, this is 1 e/s = X MJy/sr. Thus 1 MJy/sr = 1/X e/s
    etomjysr_dict = {}
    for filter in filters:
        etomjysr_dict[filter] = romanisim.bandpass.etomjysr(filter, sca=chip)
    
    # Constants from romanmask.c
    DARKERR = 0.00367 # electrons/sec residual

    for fname in filenames:
        print(f"Processing {fname}...")
        # make a copy of the file first, get the filename without the path
        filename = os.path.basename(fname)
        if os.path.abspath(fname) != os.path.abspath(filename):
             os.system(f"cp {fname} {filename}")
        
        try:
            with fits.open(filename, mode='update') as hdul:
                # Check directly in the primary header or extension 1 if primary is empty
                # DOLPHOT usually expects the image in the first extension or primary? 
                # romanmask.c checks ext[0].X/Y which implies the image is in the primary extension usually, 
                # or possibly extension 1 if it's a multi-extension FITS handled by its C library.
                # Standard convention for simplified files: PrimaryHDU has data.
                
                # Check for DOL_ROMN
                if 'DOL_ROMN' in hdul[0].header:
                    print(f"{fname} already run through romanmask (DOL_ROMN present)")
                    continue

                # Determine header to work on. If ext 1 has data, use it, else use 0.
                if len(hdul) > 1 and hdul[1].data is not None:
                     process_hdu = hdul[1]
                else:
                     process_hdu = hdul[0]

                header = process_hdu.header
                data = process_hdu.data

                # Convert the flux from MJy/sr to electrons/second
                etomjysr = etomjysr_dict.get(header.get('FILTER'), 1.0)
                print(f"etomjysr: {etomjysr}")
                data *= etomjysr**-1 # now in electrons/second

                # Get Exposure Time
                if exptime:
                    use_exptime = exptime
                else:
                    use_exptime = header.get('EXPTIME')
                    if use_exptime is None:
                        # Fallback for romanisim raw/standard files if not in primary
                        print("EXPTIME not found, defaulting to 1.0 (use -exptime to override)")
                        use_exptime = 1.0
                
                # Get Read Noise
                if rdnoise:
                    use_rdnoise = rdnoise
                else:
                    use_rdnoise = header.get('RDNOISE')
                    if use_rdnoise is None:
                        # Try finding generic equivalent
                        use_rdnoise = 0.0
                        print("RDNOISE not found, defaulting to 0.0 (use -rdnoise to override)")
                print(f"rdnoise: {use_rdnoise}")

                # Dark Current calculation
                dark_term = 0.0
                if not nodark:
                    dark_term = DARKERR * use_exptime
                
                # Effective Read Noise Calculation (romanmask.c line 43)
                # RN = sqrt( RN*RN*NCOMBINE + DARK*DARK)
                # Assuming NCOMBINE=1 for single image
                eff_rn = np.sqrt(use_rdnoise**2 + dark_term**2)
                print(f"Effective readout noise = {eff_rn:.2f} electrons")

                # Min/Max Calculation (before scaling? romanmask does get cards, get max min, then exptime mult)
                # Wait, romanmask.c:
                # 1. ROMANgetcards (calculates RN etc)
                # 2. ROMANgetMaxMin (calculates DMIN/DMAX on RAW data)
                # 3. ROMAN_exptime_mult (multiplies data AND DMIN/DMAX by EXP)
                # So we should calculate min/max on current data, then scale everything.

                dmin = np.min(data)
                dmax = np.max(data)
                
                if dmax > 10:
                    dmax *= 1.1
                else:
                    dmax += 1.0
                
                if dmin < -10:
                    dmin *= 1.1
                else:
                    dmin -= 1.0

                # Scale Data
                print(f"Scaling data by EXPTIME={use_exptime}...")
                data *= use_exptime
                dmin *= use_exptime
                dmax *= use_exptime

                # Update Data
                process_hdu.data = data # now everything is in electrons

                # Inject Headers
                # insertcards(fits.ext,1.,RN,EXP,DMIN,DMAX,EPOCH,0.0,EXP0);
                # insertcards arguments: gain, rn, exp, dmin, dmax, epoch, airmass, exp0
                # From fits_lib.c (inferred):
                # GAIN, RDNOISE, EXPTIME, MINVAL, MAXVAL, EPOCH, AIRMASS, EXP0?
                # Actually checking dolphot headers usually:
                
                header['GAIN'] = (1.0, 'Gain (e-/ADU)') # Roman is usually 2.0
                header['RDNOISE'] = (eff_rn, 'Effective readout noise (e-)')
                header['EXPTIME'] = (use_exptime, 'Exposure time (s)')
                header['MINVAL'] = (dmin, 'Minimum pixel value')
                header['MAXVAL'] = (dmax, 'Maximum pixel value')
                header['EPOCH'] = (0.0, 'Epoch')
                # AIRMASS is 0.0 in romanmask call
                # EXP0 is exptime/ncombine -> exptime
                
                # Tag
                # DOL_ROMN tag. Chip mapping:
                # romanmask checks detector name. Here we force it.
                # Integer 0-17. args.chip is 1-18.
                dol_romn_val = chip - 1
                header['DOL_ROMN'] = (dol_romn_val, 'DOLPHOT ROMAN tag')

                print(f"Writing {filename}...")
            
        except Exception as e:
            print(f"Error processing {fname}: {e}")
            # sys.exit(1) # Don't exit on one failure if processing multiple?
            continue

if __name__ == "__main__":
    fire.Fire(prep)


# Example:

# python /home/jiaxuanl/Research/Rosesim/rosesim/dolphot/prep.py ./raw/F158_642s.fits ./raw/F106_642s.fits --exptime 642 --rdnoise 5 --chip 2 --nodark
# calcsky F158_642s 10 25 -128 2.25 2.00 ; calcsky F106_642s 10 25 -128 2.25 2.00
# dolphot ngc253_642s.phot -p/home/jiaxuanl/Research/Rosesim/rosesim/dolphot/phot_642s.param > ngc253_642s.log
# fakelist ngc253_642s.phot Roman_F106 Roman_F158 20 29 -0.3 2.2 -nstar=100000 > fake_642s.inputlist
# time dolphot ngc253_642s.phot -p/home/jiaxuanl/Research/Rosesim/rosesim/dolphot/fake_642s.param > fake_642s.log


# python /home/jiaxuanl/Research/Rosesim/rosesim/dolphot/prep.py ./raw/F158_5136s.fits ./raw/F106_5136s.fits --exptime 5136 --chip 2 --nodark
# calcsky F158_5136s 10 25 -128 2.25 2.00 ; calcsky F106_5136s 10 25 -128 2.25 2.00
# dolphot ngc253_5136s.phot -p/home/jiaxuanl/Research/Rosesim/rosesim/dolphot/phot_5136s.param > ngc253_5136s.log
# fakelist ngc253_5136s.phot Roman_F106 Roman_F158 20 29 -0.3 2.2 -nstar=100000 > fake_5136s.inputlist
# dolphot ngc253_5136s.phot -p/home/jiaxuanl/Research/Rosesim/rosesim/dolphot/fake_5136s.param > fake_5136s.log


# python /home/jiaxuanl/Research/Rosesim/rosesim/dolphot/prep.py ./raw/F158_642s.fits ./raw/F106_642s.fits --exptime 642 --rdnoise 5 --chip 2 --nodark
# calcsky F158_642s 10 25 -128 2.25 2.00 ; calcsky F106_642s 10 25 -128 2.25 2.00
# dolphot cena_642s_sig5.phot -p/home/jiaxuanl/Research/Rosesim/rosesim/dolphot/phot_642s.param > cena_642s_sig5.log
# fakelist cena_642s.phot Roman_F106 Roman_F158 20 29 -0.3 2.2 -nstar=100000 > fake_642s.inputlist
# fakelist cena_5136s.phot Roman_F106 Roman_F158 20 29 -0.3 2.2 -nstar=100000 > fake_5136s.inputlist
# dolphot cena_642s.phot -p/home/jiaxuanl/Research/Rosesim/rosesim/dolphot/fake_642s.param > fake_642s.log
# dolphot cena_5136s.phot -p/home/jiaxuanl/Research/Rosesim/rosesim/dolphot/fake_5136s.param > fake_5136s.log