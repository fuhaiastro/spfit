;+
; NAME:
;       MATCH_SPEC_RES
;
; PURPOSE:
;       Match the spectral resolution of a set of input spectra to a new
;       spectral resolution.  Both spectral resolutions (existing and target)
;       can both be wavelength dependent.
;
;       In the case where the existing resolution is LOWER than the target
;       resolution, match the existing resolution to the target resolution 
;       up to some constant offset that must be accounted for in subsequent analyses
;
; CALLING SEQUENCE:
;       MATCH_SPEC_RES, flux, wave, sres, target_wave, target_sres
;
; INPUTS:
;       flux dblarr[T][S]
;               Array containing T spectra, each with S spectral channels.
;
;       wave dblarr[T][S] or dblarr[S]
;               Wavelength of each spectral channel S in angstroms for each
;               spectrum T.  Default is that the spectra are linearly sampled in
;               wavelength. 
;
;       sres dblarr[T][S] or dblarr[S]
;               Spectral resolution (R=lamda/delta lambda) for of each spectral
;               channel S for each spectrum T.  Replaced upon output with the
;               matched resolution of the spectra.
;
;       target_wave dblarr[C]
;               Wavelength coordinates in angstrom of each spectral channel C 
; 		over which the target spectral resolution is sampled. 
;
;       target_sres dblarr[C]
;               Target spectral resolution (R=lamda/delta lamba) as a function
;               of wavelength.
;
; OUTPUTS:
;	flux_new dblarr[T][S]
;		The spectra array with a resolution matched
;               (as best as possible) to the target value.

; OPTIONAL KEYWORDS:
;
; REVISION HISTORY:
;
;       29 Jul 2016: (HF) Original Implementation w/
; 	gaussian_filter1d.pro
;-
;------------------------------------------------------------------------------
pro match_spec_res, flux, wave, sres, target_wave, target_sres, flux_new

; error action
ON_ERROR, 2

; interpolate target R to the same wavelength grid
sres2 = sres
ind = where(wave gt min(target_wave) and wave lt max(target_wave))
sres2[ind] = interpol(target_sres,target_wave,wave[ind])
; SRES = lambda/FWHM = R ===> FWHM = lam/sres
dlam = mean(wave[1:*]-wave) ; linear dispersion, A/pix
sig2 = wave/sres2/dlam/2.35482 ; pix
sig  = wave/sres/dlam/2.35482  ; pix
sigdiff = sqrt((sig2^2-sig^2)>0)

flux_new = flux*0
dims = (size(flux))[0:2]
if dims[0] eq 1 then begin
	flux_new = gaussian_filter1d(flux,sigdiff)
endif else if dims[0] eq 2 then begin
	for i=0,dims[1]-1 do $
		flux_new[i,*] = gaussian_filter1d(flux[i,*],sigdiff)
endif else message,'Input flux array dimension not supported!'

end
