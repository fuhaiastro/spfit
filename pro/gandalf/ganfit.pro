PRO GANFIT, galaxyin, errorin, log10lamin, objz, templatesin, log10lam_tpl,$
	EMISSION_SETUP_IN=emission_setup_in, BITMASK=bitmask, SRES_DATA=sres_data, $
   	DEGREE=degree, MDEGREE=mdegree, SMOMENTS=smoments, GMOMENTS=gmoments, BIAS=bias, $
	EBV_GAL=ebv_gal, REDDENING=reddening, PLOT=plot, QUIET=quiet, $ 
	FIT_RESULTS=fit_results, EMISSION_SETUP_out=emission_setup_out 
;+
; NAME
;	GANFIT
;
; PURPOSE
; 	Calls PPXF and GANDALF consecutively to fit any log10 rebinned spectrum 
;	(e.g., SDSS/MaNGA). The stellar continuum is matched with a linear 
;	combination of stellar population models, whereas emission-lines are 
; 	represented by Gaussian templates, with interdependencies regulated by 
;	the input emission-setup file
;
; INPUT
;	GALAXY 	- [N_pix] array, observed flux in f_lambda unit (1d-17 erg/s/cm2/A)
;	ERROR 	- [N_pix] array, 1-sigma error of the observed flux (1d-17 erg/s/cm2/A)
;	LOG10LAM - [N_pix] array, log10 of the observed wavelengths in observed frmae (A)
;	objz 	- constant, systemic redshift of the galaxy
;	TEMPLATES - [N_pix_tpl, N_tpl] array, L_lambda unit (erg/s/A/Msun)
;	LOG10LAM_TPL - [N_pix_tpl] array, log10 of template wavelengths
;		(rest-frame, A), must have the same spacing as LOG10LAM
;	EMISSION_SETUP_IN - structure, defined from the emission line setup file
;	BITMASK - [N_pix] array, bad pixel mask, 0: good, non-zero:bad
;	SRES_DATA - [N_pix] array, spectral resolution array (R = lambda/FWHM) of 
;		the instrument. This is used to calculate INT_DISP
;		(in km/s) array for the emission line templates (km/s) 
;
; OPTIONAL INPUT
;	DEGREE - degree of the Legendre polynomial used for correcting
;       	the template continuum shape during the fit (default: -1).
;       	This correction is ADDITIVE. Controls both pPXF and GANDALF.
;     		- IMPORTANT (GANDALF): Additive polynomials cannot be used when
;       	the REDDENING keyword is set.
;	MDEGREE - degree of the *multiplicative* Legendre polynomial (with mean of 1)
;       	used to correct the continuum shape during the fit (default: 0). The
;       	zero degree multiplicative polynomial is always included in the fit as
;       	it corresponds to the weights assigned to the templates.
;       	Note that the computation time is longer with multiplicative
;       	polynomials than with the same number of additive polynomials.
;     		- IMPORTANT (PPXF & GANDALF): Multiplicative polynomials cannot be 
;		used when the REDDENING keyword is set.
;		- Even if not set, PPXF/GANDALF will assign MDEGREE=0
;	SMOMENTS: Order of the Gauss-Hermite moments to fit for SSPs.
;		Default is set to 4 to fit [V,sigma,h3,h4], set it to 2
;		for only [V, sigma]. 
;	GMOMENTS: Order of the Gauss-Hermite moments to fit for emission
;		lines (PLACEHOLDER: currently not used in GANFIT)
;	BIAS: pPXF parameter to penalize h3, h4 ...
;	/PLOT - set to plot the best fitting solution at the end of pPXF and GANDALF fit.
;	EBV_GAL - E(B-V) of foreground Galactic extinction to be applied before any fitting.
;  	REDDENING - allows to include in the fit the effect of reddening by
;       	dust, by specifying a single E(B-V) guess for extinction, or a
;       	two-element array of E(B-V) guesses. A single guess will
;       	trigger the use of a single-screen dust model, affecting the
;       	entire spectrum, whereas passing a two-elements array of
;       	guesses will add a second dust component affecting only the
;       	nebular fluxes. This second option is recommended only when
;       	recombination lines are clearly detected, so that a
;       	temperature-dependent prior (e.g. the decrement of the Balmer
;       	lines) on their relative intensity can be use to constrain
;       	such a second, internal component. 
;	/QUIET - set this keyword to mute the output
;
; OUTPUT
;	FIT_RESULTS - a structure contains the initial and best-fit
;		parameters of gas and stellar continuum
;	EMISSION_SETUP_OUT - structure, updated with best-fit parameters
;
;
; HISTORY
;	2015/5/22 HF - Written 
;	2015/7/30 HF - now takes emission_setup_in structure instead of
;		reading the linelist. also produces emission_setup_out
;		which contains the best-fit kinematics. 
;-

c = 299792.4580d ; Speed of light in km/s

IF NOT KEYWORD_SET(DEGREE) THEN degree=-1
IF NOT KEYWORD_SET(MDEGREE) THEN mdegree=0
IF NOT KEYWORD_SET(ebv_gal) THEN ebv_gal=0.0
if ~keyword_set(smoments) then smoments = 4
if ~keyword_set(gmoments) then gmoments = 2 ; currently not used

; rename to new variables to preserve the input data; otherwise the
; input data will be alterred upon return
galaxy = galaxyin
error = errorin
log10lam = log10lamin
templates = templatesin

; De-redden the spectra for Galactic extinction
IF ebv_gal gt 0 THEN BEGIN
    l0_gal   = log10lam[0] 
    lstep_gal = mean(log10lam[1:*]-log10lam) 
    dereddening_attenuation = DUST_CALZETTI(l0_gal,lstep_gal,n_elements(galaxy),-ebv_gal,0.0d,/log10)
    galaxy = galaxy*dereddening_attenuation
    error = error*dereddening_attenuation
ENDIF

; Deredshift observed galaxy spectrum 
; lambda_rf = lambda_obs/(1+z)
; flam_rf = flam_obs * (1+z) so that lambda*flam is conserved 
log10lam = log10lam - alog10(1+objz)
galaxy = galaxy * (1+objz)
error  = error * (1+objz)

; estimate stellar mass from weights
; flam_obsved (1d-17 erg/s/A/cm) = templates (erg/s/A/Msun) * Mass / fl2ll = templates * weights
; Mass = fl2ll * weight
; flam (1d-17 erg/s/A/cm2) <-> Llam (erg/s/A)
pc = 3.0856776d18 ; pc in cm
fl2ll = 1d-17 * 4 * !dpi * (lumdist(objz,/silent)*1d6*pc)^2

; Computing the velocity offset between the starting wavelengths
; of the galaxy and template spectra. You need to include a alog(10)
; factor since we are using log10 rebinned data instead of ln rebinned.
;
; The galaxy and the template spectra do not have the same starting wavelength.
; For this reason an extra velocity shift DV has to be applied to the template
; to fit the galaxy spectrum. We remove this artificial shift by using the
; keyword VSYST in the call to PPXF below, so that all velocities are
; measured with respect to DV. 
l0_gal   = log10lam[0]     ; log10(wavelength in A) 
l0_templ = log10lam_tpl[0] ; log10(wavelength in A)
offset = -(l0_gal-l0_templ)*c*alog(10.0d) 
; dV = z*c ~= ln(lambda1/lambda0)*c = log10(lambda1/lambda0)*ln(10)*c

; Preamble to PPXF
; this is the velocity scale spanned by each log_lambda interval
; You need to include a alog(10) factor as well.
lstep_gal = mean(log10lam[1:*]-log10lam) 
velscale  = c*lstep_gal*alog(10.0d)
; Initial V and sigma guesses, in km/s. For V we use the
; previously derived artificial velocity offset between 
; templates and data. For stellar velocity dispersion we 
; assume 150 km/s
start  = [offset,150.0d]

; load original input emission setup structure, to be used for masking the
; emission-line contaminated regions and fitting the emission lines in
; GANDALF. Modify emission_setup instead of emission_setup_in so that
; the latter is preserved upon return
emission_setup = emission_setup_in 

; set action='i' (ignore) for lines outside of wavelength range
npix = n_elements(galaxy)
meml_cpix = (alog10(emission_setup.lambda)-l0_gal)/lstep_gal
ind = where(meml_cpix le 0 or meml_cpix ge npix-1,ct)
if ct gt 0 then emission_setup.action[ind] = 'i'
; output emission line setup structure
; Note: do this after changing lines outside of range to action='i' 
emission_setup_out = emission_setup

; For pPXF, mask all lines with action = 'f' or 'm'
i_fit = where(emission_setup.action eq 'f') 
emission_setup.action[i_fit] = 'm' ; toggle 'f' to 'm' for pPXF
l_rf_range = 10.^minmax(log10lam)
goodpixels = mask_emission_lines(npix,0.0,emission_setup,velscale,$
                                 l0_gal,lstep_gal,sigma=200.0,/log10,l_rf_range=l_rf_range)
emission_setup.action[i_fit] = 'f' ; toggle 'm' back to 'f'
; mask additional bad pixels using input bitmask (0 - good, 1 - bad pixel)
if n_elements(bitmask) ne 0 then begin
	ind = where(bitmask[goodpixels] eq 0,ct)
	if ct gt 0 then goodpixels = goodpixels[ind]
endif

; PPXF fit! Fit only V* and sigma*, best-fit to be supplied to GANDLAF
; Note: GANDALF does not fit for V* and sigma*
; Note: After the fit the input REDDENING is replaced with the best fitting E(B-V) value.
; 	This is not the case for GANDALF
N_TEMPL = (size(templates))[2]
if keyword_set(reddening) then ebv_ppxf = reddening[0] > 0.001 $
	else ebv_ppxf=0.001
PPXF, templates, galaxy, error, velscale, start, sol, error=sol_star_err,$
    	goodpixels=goodpixels, moments=smoments, degree=degree, mdegree=mdegree,$
	reddening=ebv_ppxf, lambda=10.^log10lam, bestfit=bestfit,$
	quiet=quiet, plot=plot, bias=bias
; compute best-fit Chi^2
chi2ppxf = robust_sigma((galaxy[goodpixels]-bestfit[goodpixels])/error[goodpixels],/zero)^2
if ~keyword_set(quiet) then print,' PPXF Chi^2/DOF = ', chi2ppxf

; For GANDALF, unmask emission lines by re-assign the goodpixels array
; masking what are trully action='m' (e.g., Na absorption lines)
goodpixels = mask_emission_lines(n_elements(galaxy),sol[0]-offset,emission_setup,velscale,$
                                 l0_gal,lstep_gal,sigma=200.0,/log10,l_rf_range=l_rf_range)
; mask bad pixels using input bitmask (0 - good, 1 - bad pixel)
if n_elements(bitmask) ne 0 then begin
	ind = where(bitmask[goodpixels] eq 0,ct)
	if ct gt 0 then goodpixels = goodpixels[ind]
endif

; Assign the stellar systemic velocity as initial guess for the gas
; kinematics, adding also the velocity offset specified in the
; emission-line setup. This makes it easier to add blue or red wings to
; the line profile.
i_na = where(emission_setup.i lt 200) ; narrow lines only
emission_setup.v[i_na] = sol[0]-offset + emission_setup.v[i_na]
; save pPXF stellar kinematics
sol_star = sol

; GANDALF_HF parameter INT_DISP: 
; set the instrumental broadening of emission lines
; R = lambda/FWHM = lambda/(2.355*sigma)
; V_disp = sigma/lambda * c = c/(2.355*R)
int_disp = interpol(c/sres_data/2.355,log10lam,alog10(emission_setup.lambda)) ; km/s

; Call Gandalf, giving it only the stellar kinematics as input
; sol. Now include reddening. Note the /log10 is turned on
GANDALF_HF, templates, galaxy, error, velscale, sol, emission_setup, $
  l0_gal, lstep_gal, degree=degree, mdegree=mdegree, $
  GOODPIXELS=goodpixels, INT_DISP=int_disp, $
  BESTFIT=bestfit, EMISSION_TEMPLATES=emission_templates, WEIGHTS=weights, $
  PLOT=plot, /LOG10, L0_TEMPL=l0_templ, $
  /FOR_ERRORS, ERROR=esol,REDDENING=reddening,quiet=quiet
; compute best-fit Chi^2
chi2gan = robust_sigma((galaxy[goodpixels]-bestfit[goodpixels])/error[goodpixels],/zero)^2
if ~keyword_set(quiet) then print,' GANDALF Chi^2/DOF = ', chi2gan

; select lines that are not doublets and action='f'
; these are the lines GANDALF actually fit and provided the solution
i_l = where(emission_setup.kind eq 'l' and emission_setup.action eq 'f') 

; Compute A/N for each line using the median noise underneath each line
resi = galaxy-bestfit
Vsys = mean(sol[2:*:4])
sigma = mean(sqrt(sol[3:*:4]^2+int_disp[i_l]^2))
sol_gas_err = noise_emission_lines(resi,error,Vsys,emission_setup,velscale,$
                                 l0_gal,lstep_gal,sigma=sigma,/log10)
sol_gas_A = sol[dindgen(n_elements(i_l))*4+1] ; amplitude
sol_gas_AoN = sol_gas_A/sol_gas_err[i_l]
; remove detected emission lines, with a constant A/N cut of 4.
spec_neat = galaxy
for i = 0,n_elements(sol_gas_A)-1 do $
    if (sol_gas_AoN[i] ge 4.0) then spec_neat -= float(emission_templates[*,i])

; pull out best-fit emission line parameters
sol_gas_F = sol[dindgen(n_elements(i_l))*4+0] ; flux
sol_gas_A = sol[dindgen(n_elements(i_l))*4+1] ; amplitude
sol_gas_V = sol[dindgen(n_elements(i_l))*4+2] ; radial velocity
sol_gas_S = sol[dindgen(n_elements(i_l))*4+3] ; intrinsic sigma
if keyword_set(reddening) then sol_EBmV  = sol[n_elements(i_l)*4:*] $
			else sol_EBmV = -99.9 ; reddening E(B-V) 
; compute Sigma as observed
sigma_obs = sqrt(sol_gas_S^2. + int_disp[i_l]^2.)/velscale ; km/s -> pix
; update emission_setup_out structure
emission_setup_out.v[i_l] = sol_gas_V
emission_setup_out.s[i_l] = sol_gas_S

; emission line only model
if ((size(emission_templates))[0] eq 1) then emission = emission_templates          
if ((size(emission_templates))[0] eq 2) then emission = total(emission_templates,2)
; 0 - bad, 1 - good, opposite to bitmask
good = intarr(n_elements(galaxy)) 
good[goodpixels] = 1
; broad lines only model
i_br = where(strpos(emission_setup.name[i_l],'_br') ge 0,ct)
if ct gt 0 then emis_br = total(emission_templates[*,i_br],2) else emis_br = galaxy*0.0 

; put best-fit models in a single array
fit_model = [[galaxy],[error],[good],[bestfit],[emission],[spec_neat],[emis_br]]

dummy = emission_setup 
if keyword_set(FOR_ERRORS) then begin
    esol_gas_F = esol[dindgen(n_elements(i_l))*4+0]
    esol_gas_A = esol[dindgen(n_elements(i_l))*4+1]
    esol_gas_V = esol[dindgen(n_elements(i_l))*4+2]
    esol_gas_S = esol[dindgen(n_elements(i_l))*4+3]
    ;
    esol_EBmV  = esol[n_elements(i_l)*4:*]

    fit_results = create_struct('i',dummy.i[i_l],'name',dummy.name[i_l],'lambda',dummy.lambda[i_l],$
                                ;'action',dummy.action[i_l],'kind',dummy.kind[i_l],'a',dummy.a[i_l],$
                                ;'v',dummy.v[i_l],'s',dummy.s[i_l],'fit',dummy.fit[i_l],$
                                'Flux',sol_gas_F,'Ampl',sol_gas_A,'Vel',sol_gas_V,'Sigma',sol_gas_S,$
                                'eFlux',esol_gas_F,'eAmpl',esol_gas_A,'eVel',esol_gas_V,'eSigma',esol_gas_S,$
                                'AoN',sol_gas_AoN,'Sigma_obs',sigma_obs,'EBmV',sol_EBmV,'eEBmV',esol_EBmV,$
				'Degree',degree,'Mdegree',mdegree,$
				'Vel_star',sol_star[0]-offset,'eVel_star',sol_star_err[0],$
                                'Sigma_star',sol_star[1],'eSigma_star',sol_star_err[1],$
				'chi2_ppxf',chi2ppxf,$ ; Chi^2/DOF from PPXF
				'Weights',weights[0:n_templ-1],$ ; weights of stellar templates
				'M_star',weights[0:n_templ-1]*fl2ll,$ ; stellar mass
				'chi2nu',chi2gan,$ 	; total chi^2
				'redshift',objz,$     	; redshift
				'EBmV_gal',ebv_gal,$    ; Galactic extinction
				'log10lam',log10lam,$   ; wavelength in rest-frame in A
				'galaxy_in',galaxyin,$  ; input flux in Flam, obs frame
				'galaxy',fit_model[*,0],$ ; deredshifted spectrum in Flam, Gal dust de-reddened 
				'err', fit_model[*,1],$ ; associated 1-sigma errors
				'good',fit_model[*,2],$ ; mask: good pixels are 1, bad are 0.
				'best',fit_model[*,3],$ ; best-fit model (stellar + emission)	
				'emis',fit_model[*,4],$ ; emission-line only spectrum
				'neat',fit_model[*,5],$ ; galaxy spectrum w/ em lines removed
				'embr',fit_model[*,6])  ; broad line only model
endif else begin
    fit_results = create_struct('i',dummy.i[i_l],'name',dummy.name[i_l],'lambda',dummy.lambda[i_l],$
                                ;'action',dummy.action[i_l],'kind',dummy.kind[i_l],'a',dummy.a[i_l],$
                                ;'v',dummy.v[i_l],'s',dummy.s[i_l],'fit',dummy.fit[i_l],$
                                'Flux',sol_gas_F,'Ampl',sol_gas_A,'Vel',sol_gas_V,'Sigma',sol_gas_S,$
                                'AoN',sol_gas_AoN,'Sigma_obs',sigma_obs,'EBmV',sol_EBmV,$
				'Degree',degree,'Mdegree',mdegree,$
				'Vel_star',sol_star[0]-offset,'eVel_star',sol_star_err[0],$
                                'Sigma_star',sol_star[1],'eSigma_star',sol_star_err[1],$
				'chi2_ppxf',chi2ppxf,$ ; Chi^2/DOF from PPXF
				'Weights',weights[0:n_templ-1],$ ; weights of stellar templates
				'M_star',weights[0:n_templ-1]*fl2ll,$ ; stellar mass
				'chi2nu',chi2gan,$ 	; total chi^2
				'redshift',objz,$     	; redshift
				'EBmV_gal',ebv_gal,$    ; Galactic extinction	
				'log10lam',log10lam,$   ; wavelength in rest-frame in A
				'galaxy_in',galaxyin,$  ; input flux in Flam, obs frame
				'galaxy',fit_model[*,0],$ ; deredshifted spectrum in Flam, Gal dust de-reddened 
				'err', fit_model[*,1],$ ; associated 1-sigma errors
				'good',fit_model[*,2],$ ; mask: good pixels are 1, bad are 0.
				'best',fit_model[*,3],$ ; best-fit model (stellar + emission)	
				'emis',fit_model[*,4],$ ; emission-line only spectrum
				'neat',fit_model[*,5],$ ; galaxy spectrum w/ em lines removed
				'embr',fit_model[*,6])  ; broad line only model

endelse

END
