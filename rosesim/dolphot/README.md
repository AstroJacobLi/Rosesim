# DOLPHOT Preparation for Roman

This repository contains a script to prepare the FITS files from Rosesim for running DOLPHOT.

Note: should set `fitsky=1` in the param file for background fields. This will make the completeness better.

### Usage
1. Generate your image using Rosesim, and write as FITS files. 
2. Make a new directory like `/scratch/gpfs/JENNYG/jiaxuanl/Data/SBF/Rosesim/DOLPHOT/dw1e7_3.5Mpc_nap_sfh`, then copy the fits files to the `raw` directory under it.
```bash
cp /scratch/gpfs/JENNYG/jiaxuanl/Data/SBF/Rosesim/dw_1e7.0_3.5Mpc_napping_sfh/*.fits /scratch/gpfs/JENNYG/jiaxuanl/Data/SBF/Rosesim/DOLPHOT/dw1e7_3.5Mpc_nap_sfh/raw/
```
3. We now prepare the files for DOLPHOT, by running the `prep.py` script. Under the `dw1e7_3.5Mpc_nap_sfh` directory, run:
```bash
python /home/jiaxuanl/Research/Rosesim/rosesim/dolphot/prep.py ./raw/*.fits --exptime 600 --rdnoise 5 --chip 2 --nodark
```
Chip number could be 1-18, corresponding to sca01-sca18. The default is 2 (sca02), which is the default in Rosesim.

4. We now run DOLPHOT. Under the `dw1e7_3.5Mpc_nap_sfh` directory, run:
```bash
calcsky F106_600s_LINEGAP5_1.0 10 25 -64 2.25 2.00
calcsky F106_600s_LINEGAP5_1.1 10 25 -64 2.25 2.00
calcsky F106_600s_LINEGAP5_1.2 10 25 -64 2.25 2.00
calcsky F106_600s_LINEGAP5_1.3 10 25 -64 2.25 2.00
calcsky F106_600s_LINEGAP5_1.4 10 25 -64 2.25 2.00
calcsky F158_600s_LINEGAP5_1.0 10 25 -64 2.25 2.00
calcsky F158_600s_LINEGAP5_1.1 10 25 -64 2.25 2.00
calcsky F158_600s_LINEGAP5_1.2 10 25 -64 2.25 2.00
calcsky F158_600s_LINEGAP5_1.3 10 25 -64 2.25 2.00
calcsky F158_600s_LINEGAP5_1.4 10 25 -64 2.25 2.00

# Now I use `fitsky=1`. This is not optimal. In crowded fields, one should use `fitsky=2`, but that's too slow for my purpose. I can try `fitsky=2` in AST. Also set align=1 (used to be 0)
dolphot dw1e7_napping_sfh_3000s.phot -p/home/jiaxuanl/Research/Rosesim/rosesim/dolphot/phot_dw1e7_sfh.param > phot_dw1e7_sfh.log

# let's try fitsky=2
dolphot dw1e7_napping_sfh_3000s_fitsky2.phot -p/home/jiaxuanl/Research/Rosesim/rosesim/dolphot/phot_dw1e7_sfh_fitsky2.param > phot_dw1e7_sfh_fitsky2.log


fakelist dw1e7_napping_sfh_3000s_fitsky2.phot Roman_F106 Roman_F158 20 29 -0.3 2.2 -nstar=100000 > fake_dw1e7_napping_sfh_3000s_fitsky2.inputlist

time dolphot dw1e7_napping_sfh_3000s_fitsky2.phot -p/home/jiaxuanl/Research/Rosesim/rosesim/dolphot/fake_dw1e7_sfh_fitsky2.param > fake_dw1e7_napping_sfh_3000s_fitsky2.log
```


## Dolphot for background fields (CenA)

Now I have simulated 5 exposures following LINEGAP5_1 dither pattern, with 1 exposure per filter. Let's try to do photometry for these to get a catalog of background sources. This would be very useful for testing 1) UFD detection, 2) stream detection, 3) SFH recovery.

The fits files are under `/scratch/gpfs/JENNYG/jiaxuanl/Data/SBF/Rosesim/sky_jaguar_trilegal_cena`. 

1. Make a new directory like `/scratch/gpfs/JENNYG/jiaxuanl/Data/SBF/Rosesim/DOLPHOT/sky_jaguar_trilegal_cena`, then copy the fits files to the `raw` directory under it.
```bash
cp /scratch/gpfs/JENNYG/jiaxuanl/Data/SBF/Rosesim/sky_jaguar_trilegal_cena/*LINEGAP*.fits /scratch/gpfs/JENNYG/jiaxuanl/Data/SBF/Rosesim/DOLPHOT/sky_jaguar_trilegal_cena/raw/
```

2. We now prepare the files for DOLPHOT, by running the `prep.py` script. Under the `sky_jaguar_trilegal_cena` directory, run:
```bash
python /home/jiaxuanl/Research/Rosesim/rosesim/dolphot/prep.py ./raw/*.fits --exptime 600 --rdnoise 5 --chip 2 --nodark
```
Chip number could be 1-18, corresponding to sca01-sca18. The default is 2 (sca02), which is the default in Rosesim.

3. We now run DOLPHOT. Under the `sky_jaguar_trilegal_cena` directory, run:
```bash
calcsky F106_600s_LINEGAP5_1.0 10 25 -64 2.25 2.00
calcsky F106_600s_LINEGAP5_1.1 10 25 -64 2.25 2.00
calcsky F106_600s_LINEGAP5_1.2 10 25 -64 2.25 2.00
calcsky F106_600s_LINEGAP5_1.3 10 25 -64 2.25 2.00
calcsky F106_600s_LINEGAP5_1.4 10 25 -64 2.25 2.00

calcsky F158_600s_LINEGAP5_1.0 10 25 -64 2.25 2.00
calcsky F158_600s_LINEGAP5_1.1 10 25 -64 2.25 2.00
calcsky F158_600s_LINEGAP5_1.2 10 25 -64 2.25 2.00
calcsky F158_600s_LINEGAP5_1.3 10 25 -64 2.25 2.00
calcsky F158_600s_LINEGAP5_1.4 10 25 -64 2.25 2.00

# this only includes 1 exposure in each band, as a test
dolphot sky_jaguar_trilegal_cena_600s.phot -p/home/jiaxuanl/Research/Rosesim/rosesim/dolphot/phot_cena_600s.param > sky_jaguar_trilegal_cena_600s.log

# all 10 exposures, fitsky=1
dolphot sky_jaguar_trilegal_cena_3000s.phot -p/home/jiaxuanl/Research/Rosesim/rosesim/dolphot/phot_cena_3000s.param > sky_jaguar_trilegal_cena_3000s.log


fakelist sky_jaguar_trilegal_cena_3000s.phot Roman_F106 Roman_F158 20 29 -0.3 2.2 -nstar=100000 > fake_cena_3000s.inputlist

time dolphot sky_jaguar_trilegal_cena_3000s.phot -p/home/jiaxuanl/Research/Rosesim/rosesim/dolphot/fake_cena_3000s.param > fake_cena_3000s.log
```

## Dolphot for background fields (NGC253)

The fits files are under `/scratch/gpfs/JENNYG/jiaxuanl/Data/SBF/Rosesim/sky_jaguar_trilegal_n253`. 

1. Make a new directory like `/scratch/gpfs/JENNYG/jiaxuanl/Data/SBF/Rosesim/DOLPHOT/sky_jaguar_trilegal_n253`, then copy the fits files to the `raw` directory under it.
```bash
cp /scratch/gpfs/JENNYG/jiaxuanl/Data/SBF/Rosesim/sky_jaguar_trilegal_ngc253/*LINEGAP*.fits /scratch/gpfs/JENNYG/jiaxuanl/Data/SBF/Rosesim/DOLPHOT/sky_jaguar_trilegal_n253/raw/
```

2. We now prepare the files for DOLPHOT, by running the `prep.py` script. Under the `sky_jaguar_trilegal_n253` directory, run:
```bash
python /home/jiaxuanl/Research/Rosesim/rosesim/dolphot/prep.py ./raw/*.fits --exptime 600 --rdnoise 5 --chip 2 --nodark
```
Chip number could be 1-18, corresponding to sca01-sca18. The default is 2 (sca02), which is the default in Rosesim.

3. We now run DOLPHOT. Under the `sky_jaguar_trilegal_n253` directory, run:
```bash
calcsky F106_600s_LINEGAP5_1.0 10 25 -64 2.25 2.00
calcsky F106_600s_LINEGAP5_1.1 10 25 -64 2.25 2.00
calcsky F106_600s_LINEGAP5_1.2 10 25 -64 2.25 2.00
calcsky F106_600s_LINEGAP5_1.3 10 25 -64 2.25 2.00
calcsky F106_600s_LINEGAP5_1.4 10 25 -64 2.25 2.00

calcsky F158_600s_LINEGAP5_1.0 10 25 -64 2.25 2.00
calcsky F158_600s_LINEGAP5_1.1 10 25 -64 2.25 2.00
calcsky F158_600s_LINEGAP5_1.2 10 25 -64 2.25 2.00
calcsky F158_600s_LINEGAP5_1.3 10 25 -64 2.25 2.00
calcsky F158_600s_LINEGAP5_1.4 10 25 -64 2.25 2.00

# all 10 exposures. now I set "fitsky=1" to make things faster. this is only good when the field is very non-crowded.
dolphot sky_jaguar_trilegal_n253_3000s.phot -p/home/jiaxuanl/Research/Rosesim/rosesim/dolphot/phot_n253_3000s.param > sky_jaguar_trilegal_n253_3000s.log

# need to do this tomorrow.
fakelist sky_jaguar_trilegal_n253_3000s.phot Roman_F106 Roman_F158 23 30 -0.3 2.2 -nstar=100000 > fake.inputlist

time dolphot sky_jaguar_trilegal_n253_3000s.phot -p/home/jiaxuanl/Research/Rosesim/rosesim/dolphot/fake_n253_3000s.param > fake_n253_3000s.log
```


## Or batch processing all files


### Reference
- Dan Weisz's instruction on running DOLPHOT on JWST data: https://dolphot-jwst.readthedocs.io/en/latest/overview/workflow.html
- DOLPHOT documentation: http://americano.dolphinsim.com/dolphot/dolphot.pdf
- DOLPHOT documentation for the Roman module: http://americano.dolphinsim.com/dolphot/dolphotRoman.pdf