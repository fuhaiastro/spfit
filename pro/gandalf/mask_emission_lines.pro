function mask_emission_lines,npix,Vsys,emission_setup,velscale,l0_gal,lstep_gal,$
                             sigma=sigma,l_rf_range=l_rf_range,log10=log10

; Return a list of goodpixels to fit that excludes regions potentially
; affected by gas emission and by sky lines. Unless the log10 keyword
; is specified, wavelength values are assumed to be ln-rebinned, and
; are defined by the l0_gal, lstep_gal, npix parameters. The position of
; gas and sky emission lines is set by the input emission_setup
; structure and the width of the mask by the sigma parameter. If a
; sigma value is not passed than the width of each line is taken from
; the emission_setup structure.
;
; The rest-frame fitting wavelength range can be manually restricted
; using the l_rf_range keyword to pass min and max observed
; wavelength. Typically used to exclude regions at either side of
; spectra.
;
; HF July 2015:
;	* removed sections for sky lines because the galaxy spectra
;	have been deredshifted before masking
;	* use sigma or emission_setup.s[i], whichever is greater 
;	


; speed of light
c = 299792.458d
; define good pixels array
goodpixels = range(0,npix-1) 
; if set, exclude regions at either ends of the spectra using the keyword l_rf_range
if keyword_set(l_rf_range) then begin
    pix0     = ceil((alog(l_rf_range[0])-l0_gal)/lstep_gal+Vsys/velscale)
    pix1     = ceil((alog(l_rf_range[1])-l0_gal)/lstep_gal+Vsys/velscale)
    if keyword_set(log10) then begin
        pix0     = ceil((alog10(l_rf_range[0])-l0_gal)/lstep_gal+Vsys/velscale)
        pix1     = ceil((alog10(l_rf_range[1])-l0_gal)/lstep_gal+Vsys/velscale)
    endif
    goodpixels = range(max([pix0,0]),min([pix1,npix-1])) ; NEW - V1.3 
endif

tmppixels  = goodpixels
; looping over the listed emission-lines and mask those tagged with an
; 'm' for mask. 
for i = 0,n_elements(emission_setup.i)-1 do begin
    if (emission_setup.action[i] eq 'm') then begin
        if ~keyword_set(log10) then $
          meml_cpix = ceil((alog(emission_setup.lambda[i])-l0_gal)/lstep_gal+Vsys/velscale) $
        else $
	  meml_cpix = ceil((alog10(emission_setup.lambda[i])-l0_gal)/lstep_gal+Vsys/velscale)
        ; set the width of the mask in pixels using either
        ; 3 times the sigma of each line in the emission-line setup 
        ; or the provided sigma value (only if it is greater than the
	; sigma in the emission setup) 
        if keyword_set(sigma) then $
		msigma = 3*(sigma>emission_setup.s[i])/velscale
        if not keyword_set(sigma) then msigma = 3*emission_setup.s[i]/velscale
        meml_bpix = meml_cpix - msigma
        meml_rpix = meml_cpix + msigma
        w = where(goodpixels ge meml_bpix and goodpixels le meml_rpix) 
        if (w[0] ne -1) then begin
            tmppixels[w] = -1
        endif
   endif
endfor
w = where(tmppixels ne -1)
goodpixels = goodpixels[w]

return,goodpixels
end
