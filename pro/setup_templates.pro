pro setup_templates, wave_data, ssplib, matchR=matchR, sres_data=sres_data,$
	MH=MH, age=age, ssp_matched=ssp_matched, quiet=quiet
;+
; PURPOSE
;	Convert SSP templates to the format required by GANDALF or PPXF
;	1. match spectral resolution w/ the observed spectra in 
;		rest-frame w/ MATCH_SPEC_RES.pro
;	2. log-rebin in wavelength and match the input wavelength array
;		in dlogLam (equivalent to velScale)
;
; INPUT:
; 	wave_data: [npix] array, deredshifted, log-rebined wavelength in units of AA, for the data
; 		this is used to set the velScale for log_rebin
; 
; 	ssplib: string defines the SSP library. 
;		options are 'miuscat' and 'm11_stelib', both use Kroupa IMF
;		https://trac.sdss.org/wiki/MANGA/Data/DAPDevelopment/TemplateLibraries
; 
; OPTIONAL INPUT
; 	/matchR - match spectral resolution to sres
;		if not set, match to FWHM_SDSS = 2.76 A
;
;	sres_data - [npix] array of resolution to match the SSPs to
;	MH - 2-element array giving the range of [M/H] to be included 
;	Age - 2-element array giving the range of Age (Gyr) to be included 
;
; OUTPUT
; 	ssp_matched - resolution and dlogLam matched structure of the SSPs
;		The number of pixels depend on the wavelength range of
;		the SSP template and the velScale of the input wavelength array
;		The structure contains the following tags:
;		Age (Gyr), [Fe/H], wavelength (A), Flam (erg/s/A/Msun), Source
;
; HISTORY
;	2015/5/14 - written - Hai Fu (UIowa)
;	2016/7/29 - replaced mdap_match_spectral_resolution with
;		match_spec_res.pro, huge speed gain
;
;-

; load SSP templates, SSP wavelength is linearly sampled
if ssplib eq 'miuscat' then begin
	tpldir = getenv('SPFIT_DIR')+'ssp/MIUSCAT/MIUSCAT_BaSTI_ku_1.30_fits'
	tpl = mrdfits(tpldir+'.fits',1,h_tpl,/silent)
endif
if ssplib eq 'miuscat-thin' then begin
	tpldir = getenv('SPFIT_DIR')+'ssp/MIUSCAT/MIUSCAT_BaSTI_ku_1.30_fits-thin'
	tpl = mrdfits(tpldir+'.fits',1,h_tpl,/silent)
endif
if ssplib eq 'm11_stelib' then begin
	tpldir = getenv('SPFIT_DIR')+'ssp/M11/SSP_M11_STELIB'
	tpl = mrdfits(tpldir+'_kr.fits',1,h_tpl,/silent)
endif
FWHM_ssp = sxpar(h_tpl,'FWHM') ; Template resolution in Angstrom

; restrict the range of [M/H] and Age 
if keyword_set(MH) then begin
	ind = where(tpl.MH ge MH[0] and tpl.MH le MH[1])
	tpl = tpl[ind]
endif
if keyword_set(age) then begin
	ind = where(tpl.age ge age[0] and tpl.age le age[1])
	tpl = tpl[ind]
endif

; derive velocity scale, dV = cz = c * ln(1+dlambda/lambda)
c = 299792.458d ; speed of light in km/s
log10lam = alog10(wave_data)
lstep_gal = mean(log10lam[1:*]-log10lam) 
velscale  = c*lstep_gal*alog(10.0d)

; conversion
nssp = n_elements(tpl)
if ~keyword_set(matchR) then begin
	print,'Warning: not matching spectral resolution; use /matchR to enable'
	for i=0,nssp-1 do begin
		wave_ssp = tpl[i].wave ; linearly sampled, in AA
		flux_ssp = tpl[i].flam ; erg/s/A/Msun
		; natural logarithmically bin to the same dV as the data
		lamRange = minmax(wave_ssp) 
		log_rebin,lamRange,flux_ssp,sspNew,logLam,VELSCALE=velScale
		if i eq 0 then templates = sspNew else templates = [[templates],[sspNew]]
	endfor
endif else begin
	; match spectral resolution w/ MDAP code 
	flux = transpose(tpl.flam) ; [T,S] array, T spectra, each w/ S pixels replaced upon output
	wave = tpl[0].wave 	; template wavelength, linearly sampled
	sres = wave/FWHM_ssp 	; template resolution 
	target_wave = wave_data ; data wavelength array, logarithmically sampled
	target_sres = sres_data
	match_spec_res,flux,wave,sres,target_wave,target_sres,flux_matched
	; natural logarithmically bin to the same dV as the data
	lamRange = minmax(tpl[0].wave)
	for i=0,nssp-1 do begin
		log_rebin,lamRange,reform(flux_matched[i,*]),sspNew,logLam,VELSCALE=velScale
		if i eq 0 then templates = sspNew else templates = [[templates],[sspNew]]
	endfor
endelse

; store results in a new structure
; units: Gyr, [M/H], A, erg/s/Mo/A, none
npix = (size(templates))[1]
nssp = (size(templates))[2]
tags = ['age','MH','wave','flam','source']
vals = ['0.0','0.0','fltarr('+strc(npix)+')','fltarr('+strc(npix)+')','""']
ssp_matched = mrd_struct(tags,vals,nssp)
ssp_matched.age = tpl.age
ssp_matched.MH = tpl.MH
ssp_matched.source = tpl.source
ssp_matched.wave = exp(logLam)
ssp_matched.flam = templates

return

end
