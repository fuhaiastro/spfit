------------------------------
The pPXF software distribution
------------------------------

This distribution contains an IDL (http://www.ittvis.com/idl/)
implementation of the Penalized Pixel-Fitting (pPXF) method developed by
Cappellari M., & Emsellem E., (2004, PASP, 116, 138).

The following files are included in the software distribution:

cap_any.pro                 ---> Small auxiliary routine
cap_range.pro               ---> Small auxiliary routine
image_plot.pro              ---> Display an image with axes
log_rebin.pro               ---> logarithmically rebin a spectrum
ppxf.pro                    ---> The main pPXF routine
ppxf_determine_goodPixels.pro      ---> Example for determining goodPixels vector
ppxf_emission_lines.pro     ---> generates gas emission-line templates
ppxf_kinematics_example_sauron.pro ---> pPXF kinematics usage example
ppxf_kinematics_example_sdss.pro   ---> pPXF kinematics usage example
ppxf_population_example_sdss.pro ---> Example for stellar population
ppxf_population_gas_example_sdss.pro ---> Stellar population without masking gas
ppxf_simulation_example.pro ---> Example Monte Carlo simulation
ppxf_two_components_example.pro ---> Fit two kinematic stellar components
ppxf_idl_reference_output.txt ---> Reference output from the PPXF examples
readme.txt                  ---> This README file
/spectra                    ---> directory of FITS spectra for example
/miles_models               ---> directory of FITS SSP models

------------------
pPXF usage example
------------------

To learn how to use the main program PPXF run the example programs
and read the detailed documentation at the top of the file ppxf.pro

To run the examples:

1. extract all the files in the directory ppxf
2. cd to that directory (within IDL: CD, 'ppxf')
3. type e.g. ppxf_kinematics_example at the IDL prompt

The procedure PPXF needs the following extra IDL routines,
which are not included in the distribution:

- BVLS: by M. Cappellari -> http://purl.org/cappellari/software
- MPFIT: by C.B. Markwardt -> http://purl.com/net/mpfit
- ROBUST_SIGMA: by H. Freudenreich -> http://idlastro.gfsc.nasa.gov/

The IDL Astronomy User's Library is assumed to be installed
(http://idlastro.gfsc.nasa.gov/). In particular the programs
PPXF_*_EXAMPLE use the routines FITS_READ, SXPAR, PSF_GAUSSIAN,
while PPXF uses ROBUST_SIGMA and CGPLOT.

The program was tested on IDL 5.4-8.3.

-------------------------------
IMPORTANT: Proper usage of pPXF
-------------------------------

The PPXF routine can give sensible quick results with the default BIAS
parameter, however, like in any penalized/filtered/regularized method, the
optimal amount of penalization generally depends on the problem under study.

The general rule here is that the penalty should leave the line-of-sight
velocity-distribution (LOSVD) virtually unaffected, when it is well
sampled and the signal-to-noise ratio (S/N) is sufficiently high.

EXAMPLE: If you expect an LOSVD with up to a high h4~0.2 and your
adopted penalty biases the solution towards a much lower h4~0.1 even
when the measured sigma > 3*velScale and the S/N is high, then you
are *misusing* the pPXF method!


THE RECIPE: The following is a simple practical recipe for a sensible
determination of the penalty in pPXF:

1. Choose a minimum (S/N)_min level for your kinematics extraction and
   spatially bin your data so that there are no spectra below (S/N)_min;

2. Perform a fit of your kinematics *without* penalty (PPXF keyword BIAS=0).
   The solution will be noisy and may be affected by spurious solutions,
   however this step will allow you to check the expected mean ranges in
   the Gauss-Hermite parameters [h3,h4] for the galaxy under study;

3. Perform a Monte Carlo simulation of your spectra, following e.g. the
   included ppxf_simulation_example.pro routine. Adopt as S/N in the simulation
   the chosen value (S/N)_min and as input [h3,h4] the maximum representative
   values measured in the non-penalized pPXF fit of the previous step;

4. Choose as penalty (BIAS) the *largest* value such that, for sigma > 3*velScale,
   the mean difference between the output [h3,h4] and the input [h3,h4]
   is well within the rms scatter of the simulated values
   (see e.g. Fig.2 of Emsellem et al. 2004, MNRAS, 352, 721).

-----------------------------------
Problems with your first pPXF fit ?
-----------------------------------

Common problems with your first PPXF fit are caused by incorrect wavelength
ranges or different velocity scales between galaxy and templates. To quickly
detect these problems try to overplot the (log rebinned) galaxy and the
template just before calling the PPXF procedure.

You can use something like the following IDL lines while adjusting the
smoothing window and the pixels shift. If you cannot get a rough match
by eye it means something is wrong and it is unlikely that PPXF
(or any other program) will find a good match.

  plot, galaxy
  tmp = shift(smooth(template,3),-20)
  oplot, tmp/median(tmp)*median(galaxy), COLOR=200

################

Written: Michele Cappellari, Leiden, 6 November 2003
Last updated: MC, Sydney, 5 February 2014

