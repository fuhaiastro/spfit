
THIS README DESCRIBES THE FILES CONTAINING THE STELLAR POPULATION MODELS OF MARASTON & STROMBACK (2011, MNRAS). MODELS ARE AVAILABLE AT WWW.MARASTON.EU/M11/.
 
These models have the same stellar energetics and atmospheric parameters as the standard Maraston (2005, M05) models. The difference is in the individual stellar spectra that are assigned to the various synthetic stars according to their effective temperatures and gravities (Teff and log). These spectra can be either empirical or theoretical, and in some cases a mixture of them. 

In particular, the empirically-based models have been extended in wavelength coverage into the UV, using the high-resolution models published in Maraston et al. (2009, 493, 425), which are based on the theoretical spectral library UVBLUE (Rodrigues-Merino et al. 2005). This extension was made at various metallicities in case of the MILES-based models, and at solar metallicity for STELIB and ELODIE. 
The UV-extended models will be indicated using the addition 'UVextended' to the file name (see below).

We also provide a version of the Pickles-based models, which extend to the UV, using the models with the theoretical UV. In this case since is not an extension we call then 'UVtheoretical'.

In general, models should be chosen according to the spectral resolution that is needed by observations, though differences induced by the use of different spectral libraries exist.
See the paper for a description or contact the authors. 


The webpage contains: 

-- Simple Stellar Population models (SSPs) as a function of age, metallicity, IMF and stellar
   library. 
FILES: SSP_M11_ELODIE.tar.gz, SSP_M11_Pickles.tar.gz, SSP_M11_STELIB.tar.gz, SSP_M11_MILES.tar.gz,    
       SSP_M11_MARCS.tar.gz

These models are extensively described in 1., see below.

-- Selected SSPs (MILES and STELIB-based) in ASCII format as in the Bruzual & Charlot models
FILE: M11_BCformat.tar.gz. A README is included

-- Composite models, i.e. models with extended star formation using the SSP highres as model unit.
   
   ** LRG best fit model. This model was published in Maraston et al. 2009, MNRAS, 394, 107. It is
      a composite population made of 97% by mass solar-metallicity stars and 3% by mass metal-poor
      stars with the same age. It was found to give a good fit to the median colours of Luminous Red
      Galaxies at redshift between 0.4 and 0.7. Here we provide the same type of model using the 
      high-resolution SSPs. FILE: CSP_LRG_highres.tar.gz

   ** More composite models will be made available at a later stage.

******************************* 1. NOTATIONS OF SSPs *********************************************

GENERAL NAMING:

------------------------------------------------------------------------------
 ssp_M11'_'library'_'any eventual extension'.'imf''zxxx'.'HB morphology': 
------------------------------------------------------------------------------

  * 'ssp_M11': Simple Stellar Populations by Maraston & Stromback  2011

  * 'library': reference to the adopted stellar library, e.g. Pickles, MILES, STELIB, ELODIE,   
     MARCS, see below for a detailed description. 					
  
  * 'any eventual extension': it refers to whether the models have been extended in wavelength coverage
                              by merging with other models or have been revised with respect to the original 
 							 library. Cases are:

     'UVextended': models extended into the UV by merging with theoretically-based SSPs
     'UVtheoretical': only for Pickles-based models. The Pickles-based UV as been substituted with the
  	                  theoretical UV
     'nearIRextended': 
     'revisednearIRslope': only for MILES-based models of solar metallicity. The slope of the models 
                           longward 6000 AA has been slightly revised as described in the Appendix of the paper.

  * 'imf': Initial Mass Funcion,IMF, namely:

    '.ss' ==> Salpeter IMF
    '.kr' ==> Kroupa IMF
    '.cha' ==> Chabrier IMF

    (Other IMFs could be obtained contacting the authors)

  * 'zxxx': chemical composition, as it follows (notation identical to M05)

     z10m4=0.0001
     z0001=0.001
     z001=0.01
     z002=0.02 (solar metallicity)
     z004=0.04 	

  *  'Horizontal Branch morphology': following the notation of M05 this can be 
      rhb=intermediate-red and bhb=intermediate-blue. Relevant to sub-solar metallicities.


COLUMNS:
------------------------------------------------------------------------------
Age (Gyr) | [Fe/H] | Wavelength (AA) | Flux (ergs /s /AA /Msun) |
------------------------------------------------------------------------------

Fluxes are normalised to 1 Msun.

--------------------------------------------------
DESCRIPTION OF MODELS BASED ON EMPIRICAL LIBRARIES
--------------------------------------------------

PICKLES-BASED MODELS:	These models are based on the stellar library by Pickles (1998).
			        * Models can be calculated only for solar metallicity.
			        * Available in a purely empirical version and in an 
			          near-IR-extended version where the M-giants are 
			          semi-empirical at wavelengths shorter than 10000 AA, 
			          and purely synthetic red-wards of this. Notation: 'Pickles-ext'.
			        * TP-AGB spectra are used at original resolution, which 
			          is higher than the Pickles resolution.
			        * O-type spectra are affected by dust attenuation in 
			          the UV, which affects the youngest ages (see paper for details).
			        * A version of the model was calculated in which the UV spectrum is
	                   based on theoretical models. Notation: 'UVtheoretical'

			          Resolution: R~500
			          5.0 AA sampling

--------------------------------------------------------------------------------------------------

ssp_M11_Pickles.'imf'z002:    58 ages, 2.5 Myr to 15 Gyr
				           1150.0 - 10620 AA
				           1895 flux points

ssp_M11_Pickles_nearIRextended.'imf'z002:	58 ages, 2.5 Myr to 15 Gyr
				                        1150.0 - 25000 AA
				                        4771 flux points

ssp_M11_Pickles_UVtheoretical.ssz002:	  55 ages, 2.5 Myr to 12 Gyr
				                     1005.0 - 10620 AA 
				                     1924 flux points
				                     UV based on theoretical models, merging around 3750 AA.

ssp_M11_Pickles_UVtheoretical_nearIRextended.'imf'z002: 55 ages, 2.5 Myr to 12 Gyr
				                                       1005.0 - 25000 AA 
				                                       4800 flux points
				                                       UV based on theoretical models 
												     merging around 3750 AA.  


--------------------------------------------------------------------------------------------------

STELIB-BASED MODELS: These models are based on the stellar library STELIB by 
                     Le Borgne et al. (2003).
			      * Cool high-resolution theoretical stars (Gustafsson et 
			        al. 2008) have been added to counter inadequacies in 
			        the STELIB library. 
                    * Many spectra in the library are not complete over the quoted wavelength 					        range, which limits the actual wavelength range of the youngest 
			        ages at solar, and all ages at non-solar metallicities.
			      * Somewhat coarse sampling of stellar parameter space at 
			        non-solar metallicities may affect the reliability of these models.

			        Resolution: 3.1-3.4 AA (fwhm)
			        0.5 AA sampling


ssp_M11_STELIB.'imf'z001: 24 ages, 200 Myr to 15 Gyr
				         3201.0 - 9296.5 AA
				         12192 flux points
				         !Usable wavelength range for all ages is
				         3201.0 - 7900.0 AA!

ssp_M11_STELIB_UVextended.ssz001: 21 ages, 200 Myr to 12 Gyr
				                 1001.0 - 9296.5 AA
				                 16952 flux points
				                 !Usable wavelength range for all ages is
				                 1001.0 - 7900.0 AA!

ssp_M11_STELIB.'imf'z002: 39 ages, 30 Myr to 15 Gyr
				         3201.0 - 9296.5 AA
				         12192 flux points

ssp_M11_STELIB_UVextended.ssz002: 36 ages, 30 Myr to 12 Gyr
				                1001.0 - 9296.5 AA 
				                16592 flux points
				                UV-extended version, models merged around 3750 AA.

ssp_M11_STELIB.'imf'z004:  22 ages, 400 Myr to 15 Gyr
				        3201.0 - 9296.5 AA
				        12192 flux points
				        !Usable wavelength range for all ages is
				        3201.0 - 7900.0 AA!

ssp_M11_STELIB_UVextended.ssz004: 18 ages, 400 Myr to 12 Gyr
				               1001.0 - 9296.5 AA
				               16952 flux points
				               !Usable wavelength range for all ages is
				               1001.0 - 7900.0 AA!

--------------------------------------------------------------------------------------------------

MILES-BASED MODELS: These models are based on the stellar library MILES by 
	                Sanchez-Blazquez et al. (2006).
			      * A few cool, high-resolution theoretical stars 
			       (Gustafsson et al. 2008) have been added for 
			       completeness.
			      * Fairly poor sampling of stellar parameter space at the 
			        lowest metallicity may affect the reliability of these models.
 	               * We make two versions of solar metallicity models, one based on the MILES
       			    library as it is, and another one in which we revised the slope longward 6000
				    based on the finding that the slope of the MILES stellar spectra at the 	
					reddest lambdas 
                      is steeper (i.e. less flux) than the one of other libraries, which we 
                      figured out to be related to the assumed temperature calibration for RGB 
                      stars in 	MILES (Maraston & Stromback 2011). This manipulated model version
                      is indicated as 'revised'. Its merit with respect to the standard model      
                      needs to be assessed through comparison with data (future work).
				  * The actual spectral resolution of MILES has been determined in
				    Beifiori et al. 2011, arXiv:1012.3428

			        Resolution: 2.54 AA (fwhm)
			        0.9 AA sampling

ssp_M11_MILES.'imf'z10m4.bhb: 11 ages, 5-15 Gyr
				             3500.0 - 7429.4 AA
				             4367 flux points
							Blue Horizontal Branch

ssp_M11_MILES.'imf'z10m4.rhb: 11 ages, 5-15 Gyr
				             3500.0 - 7429.4 AA
				             4367 flux points
						    Red Horizontal Branch	 

ssp_M11_MILES.'imf'z0001.bhb: 14 ages, 2-15 Gyr
				            3500.0 - 7429.4 AA
				            4367 flux points
						   Blue Horizontal Branch

ssp_M11_MILES.'imf'z0001.rhb: 14 ages, 2-15 Gyr
				            3500.0 - 7429.4 AA
				            4367 flux points
						   Red Horizontal Branch

ssp_M11_MILES.'imf'z001: 34 ages, 55 Myr to 15 Gyr
				       3500.0 - 7429.4 AA
				       4367 flux points

ssp_M11_MILES_UVextended.'imf'z001: 31 ages, 55 Myr to 15 Gyr
				                  1000.7 - 7429.4 AA
				       			7144 flux points

ssp_M11_MILES.'imf'z002:   50 ages, 6.5 Myr to 15 Gyr
				         3500.0 - 7429.4 AA
				         4367 flux points

ssp_M11_MILES_UVextended.'imf'z002: 47 ages, 6 Myr to 15 Gyr
				                  1000.7 - 7429.4 AA 
				                  7144 flux points 
							     UV extended (models merged around 3750 AA)

ssp_M11_MILES_revisednearIRslope.'imf'z002: 50 ages, 6.5 Myr to 15 Gyr
				                          3500.0 - 7429.4 AA
				                          4367 flux points
				                          Revised near-IR slope

ssp_M11_MILES_UVextended_revisednearIRslope.'imf'z002: 47 ages, 6.5 Myr to 15 Gyr
				                                      1000.7 - 7429.4 AA 
				   		                              7144 flux points
				                                      UV extended (models merged around 3750 AA)
												    Revised near-IR slope
ssp_M11_MILES_UVextended_revisednearIRslope_nearIRextended.ssz002: 10 ages, 3 Gyr to 15 Gyr
								      1000.7 -199948.9032 AA
									     73626 flux points
								      	UV extended (models merged around 3750 AA)
												    Revised near-IR slope
									Extended into the near-IR by merging with MARCS-based models 		         

ssp_M11_MILES.'imf'z004:   25 ages, 100 Myr to 15 Gyr
				          3500.0 - 7429.4 AA
				          4367 flux points

ssp_M11_MILES_UVextended.'imf'z004: 22 ages, 100 Myr to 12 Gyr
				                   1000.7 - 7429.4 AA
				                   7144 flux points
								  UV extended (models merged around 3750 AA)


--------------------------------------------------------------------------------------------------

ELODIE-BASED MODELS: These models are based on the stellar library ELODIE by Prugniel et al. (2007).
			* Caveat regarding the flux calibration/de-reddening of 
			  the hottest stars, which will affect the youngest 
			  ages, < 100 Myr at solar metallicity.
			* A few cool, high-resolution theoretical stars 
			  (Gustafsson et al. 2008) have been added for 
			  completeness.
			* Fairly poor sampling of stellar parameter space at the 
			  lowest metallicity affects the reliability 
			  of these models.

		      Resolution: 0.55 A (fwhm)
			  0.2 AA sampling
	

ssp_M11_ELODIE.'imf'z10m4.rhb: 10 ages, 6-15 Gyr
				              3900.0 - 6800.0 AA
				              14501 flux points

ssp_M11_ELODIE.'imf'z10m4.bhb: 6 ages, 6-11 Gyr
				              3900.0 - 6800.0 AA
				              14501 flux points

ssp_M11_ELODIE.'imf'z001: 34 ages, 55 Myr to 15 Gyr
				          3900.0 - 6800.0 AA
				          14501 flux points

ssp_M11_ELODIE.'imf'z002: 57 ages, 3 Myr to 15 Gyr
				          3900.0 - 6800.0 AA
				          14501 flux points

ssp_M11_ELODIE_UVextended.'imf'z002: 54 ages, 3 Myr to 15 Gyr
				                   1000.2 - 6800.0 AA
				                   29000 flux points
				                   UV extended version, models merged around 3950 AA.

ssp_M11_ELODIE.'imf'z004: 25 ages, 100 Myr to 15 Gyr
				         3900.0 - 6800.0 AA
				         14501 flux points


----------------------------------------------------
DESCRIPTION OF MODELS BASED ON THEORETICAL LIBRARIES
----------------------------------------------------

MARCS-BASED MODELS: These models are based on the theoretical stellar library by Gustafsson et al. 2008. 
			        * Temperatures of theoretical stars allow the calculation of only old models, 
			          namely 3 Gyr and above				    
                     * The lower metallicity Z=10^-4 ('z10m4') has [alpha/enhanced]=0.4 atmospheres.
					* There will be a dedicated follow-up publication focus on these models and more
					  under computation.

    		            Resolution: R=20,000
			        0.065 AA sampling

ssp_M11_MARCS.'imf'z002: 13 ages, 3 Gyr to 15 Gyr
				        1299.5756 - 199948.9032 AA
				        100724 flux points

ssp_M11_MARCS.'imf'z004: 13 ages, 3 Gyr to 15 Gyr
				        1299.5756 - 199948.9032 AA
				        100724 flux points

ssp_M11_MARCS.'imf'z10m4: 9 ages, 7 Gyr to 15 Gyr
				         1299.5756 - 199948.9032 AA
				         100724 flux points

