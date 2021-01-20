function noise_emission_lines,resi,nois,Vsys,emission_setup,velscale,l0_gal,lstep_gal,$
                             sigma=sigma,log10=log10
;+
; Return an array of errors at the location of emission lines, computed
; using robust_sigma of the residual spec or the max of the noise
; array, whichever is larger.
;
; Modified from mask_emission_lines.pro
;-

; define pixel index array
npix = n_elements(resi)
pixind = indgen(npix)

tmperr = fltarr(n_elements(emission_setup.i))
; looping over the listed emission-lines 
for i = 0,n_elements(emission_setup.i)-1 do begin
    if ~keyword_set(log10) then $
      meml_cpix = ceil((alog(emission_setup.lambda[i])-l0_gal)/lstep_gal+Vsys/velscale) $
    else $
      meml_cpix = ceil((alog10(emission_setup.lambda[i])-l0_gal)/lstep_gal+Vsys/velscale)
    ; set the width of the mask in pixels using either
    ; +/-6 times the sigma of each line in the emission-line setup 
    ; or the provided sigma value
    ; FWHM = 2.355 sigma
    if keyword_set(sigma) then msigma = 6*sigma/velscale else $
    			   msigma = 6*emission_setup.s[i]/velscale
    meml_bpix = meml_cpix - msigma
    meml_rpix = meml_cpix + msigma
    w = where(pixind ge meml_bpix and pixind le meml_rpix,ct) 
    if ct gt 0 then $
    	; take the stddev of residual spectrum if it is larger than the input noise array
    	tmperr[i] = robust_sigma(resi[w]) > max(nois[w],/nan) 
endfor

return,tmperr
end
