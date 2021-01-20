pro fit_sdss, infile, ssplib=ssplib, MH=MH, age=age, linefile=linefile, $
	outdir=outdir, ps=ps, verbose=verbose, BLR=BLR, GANFIT=ganfit, _EXTRA=extra
;+
; NAME
;	FIT_SDSS
;
; PURPOSE
;	Top level wrapper to fit SDSS/BOSS spectra with SPFIT/GANDALF
;
; SYNTAX
;	fit_sdss,'example/spec-0278-51900-0112.fits',/verbose,reddening=0.01
;
; INPUT
; 	infile - string giving the filename of input SDSS spectrum, 
; 		e.g., example/spec-0278-51900-0112.fits
;
; OPTIONAL INPUT
; 	ssplib - string defines the SSP library to be used. 
;		Curent options include 'miuscat', 'm11_stelib' ... 
;		both use Kroupa IMF. (default: miuscat-thin)
;	linefile - emission line setup file
;	outdir - output directory (default: the same as input file)
;	/verbose - print information at each stage of gandalf_SDSS and
;		report total elapsed time
; 	/ps - save final fit result as a PS file; if set to 0, save a PNG file
;	/BLR - include AGN broad emission lines
;	/GANFIT - use GANFIT instead of SPFIT
;	_EXTRA - SPFIT/GANFIT Keywords
;
; OUTPUT
;	infile_ssplib.fits - spectral resolution matched log10-rebinned SSP 
;		templates
;	infile_spfit.fits - best-fit emission lines + stellar continuum (ext 1)
;	infile_spfit.png/ps - summary figure illustrating the quality of the fit
;
; HISTORY
;	2015/8/6 - Written - Hai Fu (UIowa)
;
; 
;-

; start counting elapsed time
if keyword_set(verbose) then tic

; set up default parameters
IF ~keyword_set(linefile) then $
	linefile = getenv('SPFIT_DIR')+'/pro/emission_lines_setup.txt'
if ~keyword_set(ssplib) then ssplib = 'miuscat-thin' ; SSP templates
; basename for output files
basename = exfilename(infile,dirname=dirname,/noextension)
; default output results in the same directory as the input file
if ~keyword_set(outdir) then outdir = dirname 	
if ~file_test(outdir) then spawn,'mkdir '+outdir
if ~keyword_set(ganfit) then cmd = 'spfit' else cmd = 'ganfit'

; Read SDSS/BOSS spectra, which is already in constant log10-lambda steps
; Data model: http://dr12.sdss3.org/datamodel/files/BOSS_SPECTRO_REDUX/RUN2D/spectra/PLATE4/spec.html
coadd = mrdfits(infile,1,hdr,/silent)
flux = coadd.flux ; log10-rebinned, coadded calibrated flux [10-17 ergs/s/cm2/Å]
var = 1./coadd.ivar ; variance flux
mask = coadd.and_mask ; 0 - good, others - bad 
wave = 10.^coadd.loglam ; log10-rebinned vacuum wavelength in AA
; the wavelengths are also shifted such that measured velocities will be 
; relative to the solar system barycentric at the mid-point of each 15-minute exposure
;
; wdisp: wavelength dispersion in pixel=dloglam units
; wavelength array (L) has been log10-rebinned, with a constant 1e-4 dex spacing
; i.e., dlogL = dL/(ln10*L) = 0.0001, therefore, dL = dlogL*ln10*L in Angstrom
; FWHM = 2.355 * WDISP * dL = 2.355 * WDISP * dlogL*ln10*L
; R = L/FWHM = 1/(2.355 * WDISP * dlogL * ln10)
lstep_gal = mean(coadd[1:*].loglam-coadd.loglam) ; = dlogL 
sres = 1/(2.355 * alog(10) * lstep_gal * coadd.wdisp)
; get galaxy information from spAll structure
spall = mrdfits(infile,2,/silent)
objz = spall.Z 		; pipeline redshift
objvdisp = spall.vdisp 	; stellar velocity dispersion in km/s
objRA = spall.plug_ra 	; RA (deg)
objDec = spall.plug_dec ; Dec (deg)

; obtain foreground Galactic extinction if not supplied as an input
if ~keyword_set(ebv_gal) then begin
	glactc,objRA,objDec,2000.,gl,gb,1,/degree
	ebv_gal = dust_getval(gl,gb,map='Ebv',/interp)
endif

; load SSP templates
; - match the spectral resolution to the data in rest-frame
; - resample the SSPs in log10 w/ the same dlogLam as the data
if keyword_set(verbose) then print,'--> loading & spectral matching template library ...'+ssplib
mtplfile = outdir+basename+'_'+ssplib+'.fits' ; saves the matched template library
if ~file_test(mtplfile) then begin
	setup_templates,wave/(1+objz),ssplib,/matchR,sres_data=sres,$
		MH=MH,age=age,ssp_matched=ssp_matched,quiet=quiet
	mwrfits,ssp_matched,mtplfile,/create
endif else begin
	ssp_matched = mrdfits(mtplfile,1,/silent)
endelse
;; select a subset of templates
;ind = where(ssp_matched.MH eq 0.06)
;ssp_matched = ssp_matched[ind]
; extract useful arrays from the structure
flam_tpl = ssp_matched.flam ; Unit: erg/s/Mo/AA
wave_tpl = ssp_matched[0].wave ; Unit: AA

; Tailor the templates and spectra within the overlapping range between templates and data
minmax_tpl = minmax(wave_tpl)
minmax_gal = minmax(wave/(1.+objz))
if keyword_set(verbose) then begin
	print,'Rest-Frame Wavelength range (templates vs. data), original'
	print,minmax_tpl,minmax_gal
end
; compute the overlapping range
tpl_range = [0.0,0.0]  
gal_range = [0.0,0.0]
buffer = 500.0 ; km/s, buffer to allow template to shift
c_kms = 299792.4580d ; Speed of light in km/s
if minmax_tpl[0] lt minmax_gal[0] then begin
	tpl_range[0] = minmax_gal[0]*(1.-buffer/c_kms)
	gal_range[0] = minmax_gal[0]
endif else begin
	tpl_range[0] = minmax_tpl[0]
	gal_range[0] = minmax_tpl[0]*(1.+buffer/c_kms) 
endelse
if minmax_tpl[1] gt minmax_gal[1] then begin
	tpl_range[1] = minmax_gal[1]*(1.+buffer/c_kms)
	gal_range[1] = minmax_gal[1]
endif else begin
	tpl_range[1] = minmax_tpl[1]
	gal_range[1] = minmax_tpl[1]*(1.-buffer/c_kms) 
endelse
w = where(wave/(1.+objz) ge gal_range[0] and wave/(1.+objz) le gal_range[1])
wave = wave[w]
sres = sres[w]
flux = flux[w]
var  = var[w]
mask = mask[w]
w = where(wave_tpl ge tpl_range[0] and wave_tpl le tpl_range[1])
wave_tpl = wave_tpl[w]
flam_tpl = flam_tpl[w,*]
if keyword_set(verbose) then begin
	print,'Rest-Frame Wavelength range (templates vs. data), after matching'
	print,minmax(wave_tpl), minmax(wave/(1+objz))
end

;; mask out Sky lines, observed-frame 5577, 6300, 6363 A
;hw = 10. ; wavelength range to mask (A)
;ind = where( (wave gt 5577-hw and wave lt 5577+hw)) ; $
;	;or (wave gt 6300-hw and wave lt 6300+hw) $
;	;or (wave gt 6363-hw and wave lt 6363+hw))
;mask[ind] = 1

; mask out M11/STELIB interpolated regions
; see Maraston 2011: 6850–6950, 7580–7700 and 8850–9050 Å
if ssplib eq 'm11_stelib' then begin
	print,'  masking out STELIB interpolated spectral regions'
	ind = where( (wave/(1+objz) gt 6850 and wave/(1+objz) lt 6950) $
		or (wave/(1+objz) gt 7580 and wave/(1+objz) lt 7700) $
		or (wave/(1+objz) gt 8850 and wave/(1+objz) lt 9050))
	mask[ind] = 1
endif

; prepare fitting
galaxy = flux 
error = sqrt(var)
bitmask = mask
; remove NaNs
ind = where(~finite(error),ct)
if ct gt 0 then begin
	error[ind] = max(error,/nan)*1e3
	galaxy[ind] = 0.0 
	bitmask[ind] = 1 ; make sure to mask out these pixels
endif

if keyword_set(verbose) then print,'--> Running '+cmd+' ...'
; reading in the emission-line setup file, and creating the corresponding structure
readcol,linefile,eml_i,eml_name,eml_lambda,eml_action,eml_kind,eml_a,$
	eml_v,eml_s,eml_fit,f='(i,a,f,a,a,f,f,f,a)',comment='#'
; this one includes broad lines
emissionbr = create_struct('i',eml_i,'name',eml_name,'lambda',eml_lambda,$
	'action',eml_action,'kind',eml_kind,'a',eml_a,'v',eml_v,'s',eml_s,'fit',eml_fit)
; this one narrow lines only
ii = where(emissionbr.i lt 200)
emissionna = create_struct('i',eml_i[ii],'name',eml_name[ii],$
	'lambda',eml_lambda[ii],'action',eml_action[ii],'kind',eml_kind[ii],$
      'a',eml_a[ii],'v',eml_v[ii],'s',eml_s[ii],'fit',eml_fit[ii])

log10lam = alog10(wave)
log10lam_tpl = alog10(wave_tpl)
; fit w/o broad lines
xx = execute(cmd+',galaxy,error,log10lam,objz,flam_tpl,log10lam_tpl,'+$
  'EMISSION_SETUP_IN=emissionna,BITMASK=bitmask,SRES_DATA=sres,EBV_GAL=ebv_gal,'+$
  'FIT_RESULTS=fit_results,EMISSION_SETUP_out=emissionout,_EXTRA=extra')
if keyword_set(BLR) then $
; fit w/ broad lines
xx = execute(cmd+',galaxy,error,log10lam,objz,flam_tpl,log10lam_tpl,'+$
  'EMISSION_SETUP_IN=emissionbr,BITMASK=bitmask,SRES_DATA=sres,EBV_GAL=ebv_gal,'+$
  'FIT_RESULTS=fit_results_br,EMISSION_SETUP_out=emissionout,_EXTRA=extra')

; save fitting results
plotfile = outdir+basename+'_'+cmd
mylegend = ssplib+' '+basename
if keyword_set(ps) then begin
	show_spec_fit,fit_results,mylegend=mylegend,/ps,outfile=plotfile+'.eps'
endif else begin
	;device,decomp=0 
	window,0,xs=600*1.5,ys=400*1.5
	show_spec_fit,fit_results,mylegend=mylegend
	save_screen,plotfile+'.png'
endelse
; repeat for broad line fit
if keyword_set(BLR) then begin
  pars = fit_results_br
  if keyword_set(ps) then begin
  	; generate PS file to exam the quality of the fit
  	show_spec_fit,pars,mylegend=mylegend,/ps,$
  		outfile=plotfile+'br.eps'
  endif else begin
  	;device,decomp=0 
  	window,0,xs=600*1.5,ys=400*1.5
  	show_spec_fit,pars,mylegend=mylegend
  	save_screen,plotfile+'br.png'
  endelse
endif

; save result in a FITS file
mwrfits,fit_results,plotfile+'.fits',hdr,/create
if keyword_set(BLR) then mwrfits,fit_results_br,plotfile+'.fits'

if keyword_set(verbose) then toc

end
