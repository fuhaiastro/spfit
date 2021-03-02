## To fit a single SDSS/BOSS spectrum

```
infile = '1237664853103607848.fits'
linefile = 'emission_lines_setup.txt'
fit_sdss,infile,linefile=linefile,/verbose,$
	ssplib='miuscat-thin',mdegree=6,/blr,/ps
```

## Output files:

The `*_miuscat-thin.fits` file saves the spec-resolution-matched stellar 
population synthesis model templates. 

The `*_spfit.fits` file saves the best-fit parameters and errors. Its
format is defined by the code below. When `/BLR` is set, the first
extension saves the best-fit result without using broad Balmer lines
(which corresponds to the figures in `*_spfit.eps`), and the second
extension saves the result that include broad Balmer lines (which
corresponds to the figures in `*_spfitbr.ps`)

```idl
fit_results = create_struct($
	;-fitted emission lines-;
	'i',emission_setup.i[i_l],'name',emission_setup.name[i_l],$
	'lambda',emission_setup.lambda[i_l],$
	;-best-fit parameters---;
        'Flux',Fgas,'EW',EW,'Vel',Vgas,'Sigma',Sgas,'H3',H3gas,'H4',H4gas,$
					  ; emission line flux (1d-17 erg/s/cm2)
					  ; V, sigma (km/s)
					  ; rest-frame Equivolent Width in Ang
        'AoN',AoN,'Sigma_obs',float(sigma_obs),$ ; emission line Amplitude-to-Noise ratio
					  ; Note: noise calculated from residual spectrum
					  ; observed line width in pixels
	'EBmV',float(EBmV),$		; intrinsic reddening and err
	'kinstar',float(kinstar),$	; stellar kinematics
	'Weights',float(w_ssp),$       ; weights of SSP templates in units of 1e-30
	'M_star',float(m_star),$       ; stellar mass
	'chi2ppxf',float(chi2ppxf),$   ; Chi^2/DOF from PPXF
	'chi2nu',float(chi2spfit),$    ; Chi^2/DOF from SPFIT
	;-spectral models-------;
	'redshift',objz,$     	; redshift
	'EBmV_gal',ebv_gal,$    ; foreground Galactic extinction
	'apoly',float(apoly),'apolypars',float(apolypars),$ ; additive polynomial and parameters
	'mpoly',float(mpoly),'mpolypars',float(mpolypars),$ ; multiplicative polynomial
	'log10lam',log10lam,$   ; wavelength in rest-frame in log(A)
	'galaxy_in',galaxyin,$  ; input flux in Flam, observed-frame (1d-17 erg/s/cm2/A)
	'galaxy',galaxy,$ 	; Flam, rest-frame & Gal dust de-reddened (1d-17 erg/s/cm2/A) 
	'err', error,$ 		; 1-sigma errors of galaxy
	'good',good,$ 		; mask: good pixels are 1, bad are 0.
	'best',float(bestfit),$ ; best-fit model (stellar + emission)	
	'emis',float(emission),$ ; emission-line model spectrum
	'neat',spec_neat,$ 	; galaxy spectrum w/ emission-line model removed
	'embr',emis_br)  	; broad line only model
```
