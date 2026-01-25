# DOLPHOT Preparation for Roman

This repository contains a script to prepare the FITS files from Rosesim for running DOLPHOT.

### Usage
```bash
python /home/jiaxuanl/Research/Rosesim/rosesim/dolphot/prep.py ./raw/F158_642s.fits ./raw/F106_642s.fits --exptime 642 --rdnoise 5 --chip 2 --nodark

calcsky F158_642s 10 25 -64 2.25 2.00 ; calcsky F106_642s 10 25 -64 2.25 2.00

dolphot roman_642s.phot -p/home/jiaxuanl/Research/Rosesim/rosesim/dolphot/phot.param > phot_642s.log

fakelist roman_642s.phot Roman_F106 Roman_F158 20 29 -0.3 2.2 -nstar=100000 > fake.inputlist

time dolphot roman_642s.phot -pfake.param > fake.log
```


### Reference
- Dan Weisz's instruction on running DOLPHOT on JWST data: https://dolphot-jwst.readthedocs.io/en/latest/overview/workflow.html
- DOLPHOT documentation: http://americano.dolphinsim.com/dolphot/dolphot.pdf
- DOLPHOT documentation for the Roman module: http://americano.dolphinsim.com/dolphot/dolphotRoman.pdf