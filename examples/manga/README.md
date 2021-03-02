## To fit a single MaNGA datacube
```idl
fit_manga,7443L,1901L,$
	infile='manga-7443-1901-LOGCUBE.fits.gz',$
	drpfile='drpall-v2_7_1.fits',$
	mdegree=6,outdir='spfit/',mtpldir='spfit_mtpl/',$
	/nobinning [,/verbose,/overwrite,/saveplot,/ps,step=100]
```
For a test run, use the options in the square brackets

## To fit a batch of datacubes 
- this mode assumes that the datacubes are downloaded to
`$MANGA_SPECTRO_REDUX/$MANGADRP_VER/`
- this mode uses N CPU-threads to fit N datacubes simultaneously using
  scripts under `pro/parallel` written by Alfred de Wijn.

```idl
; load drpall file for unique galaxy cubes
drpall = mrdfits('drpall-v2_7_1.fits',1)
; main galaxy sample: all Primary+ & full Secondary
plates = drpall.plate
ifus = long(drpall.ifudsgn)
; run 6 parallel processes
cmd = 'fit_manga'
extra = ',tag=''MPL-9'',mdegree=6,outdir=''spfit/'','+$
	'mtpldir=''spfit_mtpl/'',/nobinning,/quiet,/ignore_drp3qual'
; output files for manga_parallel to skip if alreay exist
output= 'spfit/'+strc(plates)+'-'+strc(ifus)+'.fits' 
manga_parallel,plates,ifus,cmd,extra,ncpu=6,outputs=output,/skip
```

## Output files:

The `spfit_mtpl/7443-1901_miuscat-thin.fits` file saves the spec-resolution-matched stellar 
population synthesis model templates. 

The `spfit/7443-1901.fits` file saves the best-fit parameters and errors. Its
format is defined by the code below: 
```idl
pars_str = create_struct($
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
	'chi2nu',float(chi2spfit))    ; Chi^2/DOF from SPFIT
```

The `spfit/7443-1901_spec.fits` file saves the best-fit models: 
```idl
spec_str = create_struct($
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
