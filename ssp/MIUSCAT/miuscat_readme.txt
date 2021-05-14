------------------------------------------------------------------------
MIUSCAT SSPs: extended MILES spectral coverage
(replaced by E-MILES SSPs in 2017)
------------------------------------------------------------------------

Model description: http://www.iac.es/proyecto/miles/pages/ssp-models.php
Model version: 10.0, April 2015 
Uniform resolution (FWHM = 2.51 A) from 3465-9469 A with 0.9A dispersion
Unit: The total mass of the SSP is 1Mo. 
The SSP spectra are given in units of Lo/Mo/Å, Lo = 3.826e33 erg/s

Naming convention:
	http://research.iac.es/proyecto/miles/pages/ssp-models/name-convention.php
  	Iku1.30Zm0.25T00.0300_iTp0.00_baseFe
 	Spectral Range, IMF type, IMF slope, [M/H], age (Gyr), isochrone (iP/iT)

Library downloaded from:
	http://miles.iac.es/pages/webtools/tune-ssp-models.php

Input Parameters:
	- IMF: Kroupa Universal (Ku) w/ a slope of 1.3
	- [a/Fe]: BaseFe
	- Isochrone: BaSTI (iT)

	Output: 12 M/H x 53 ages = 636 SSPs

	Z = [0.0001, 0.0003, 0.0006, 0.001, 0.002, 0.004, 0.008, 0.0100,
	0.0198, 0.0240, 0.0300, 0.0400]
	
	MH =
	[-2.27,-1.79,-1.49,-1.26,-0.96,-0.66,-0.35,-0.25,+0.06,+0.15,+0.26,+0.40]
	
	Age = [00.03, 00.04, 00.05, 00.06, 00.07, 00.08, 00.09, 00.10,
	00.15, 00.20, 00.25, 00.30, 00.35, 00.40, 00.45, 00.50, 00.60,
	00.70, 00.80, 00.90, 01.00, 01.25, 01.50, 01.75, 02.00, 02.25,
	02.50, 02.75, 03.00, 03.25, 03.50, 03.75, 04.00, 04.50, 05.00,
	05.50, 06.00, 06.50, 07.00, 07.50, 08.00, 08.50, 09.00, 09.50,
	10.00, 10.50, 11.00, 11.50, 12.00, 12.50, 13.00, 13.50, 14.00]

	- Isochrone: Padova+00 (iP)

	Output: 7 M/H x 50 ages = 350 SSPs

	Z = [0.0001, 0.0004, 0.001, 0.004, 0.008, 0.019, 0.03]

	MH = [-2.32,-1.71,-1.31,-0.71,-0.41,+0.0,+0.22]

	Age = [0.063, 0.071, 0.079, 0.089, 0.10, 0.11, 0.13, 0.14, 0.16,
	0.18, 0.20, 0.22, 0.25, 0.28, 0.32, 0.35, 0.40, 0.45, 0.50,
	0.56, 0.63, 0.71, 0.79, 0.89, 1.00, 1.12, 1.26, 1.41, 1.58,
	1.78, 2.00, 2.24, 2.51, 2.82, 3.16, 3.55, 3.98, 4.47, 5.01,
	5.62, 6.31, 7.08, 7.94, 8.91, 10.00, 11.22, 12.59, 14.13 15.85,
	17.78]

