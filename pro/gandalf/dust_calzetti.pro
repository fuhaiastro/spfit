function dust_calzetti,l0_gal,lstep_gal,npix,ebv,vstar,LOG10=log10
; This procedure uses the dust model of Calzetti et al. (2000, ApJ,
; 533, 682), and for a given E(B-V) value returns the flux attenuation
; array, which can be used to get reddened templates. Here the spectra
; are assumed to be binned on a ln-rebinned wavelentgh grid as defined
; by input l0_gal,lstep_gal,npix parameters. The input receiding
; velocity vstar, is used to derive the dust reddening in the galaxy
; rest-frame.
; 
; Can be used also to de-reddened the galaxy spectra by the Milky-Way
; dust extinction, using as E(B-V) the opposite of the Schlegel et
; al. values found in NED and vstar = 0.
;
; Initial version kindly provided by S. Kaviray, Oxford, 2006.
; changed lambda -> lambda0 to avoid ambiguity with IDL's LAMBDA function

;----------
; HF 20150810: following code is replaced to be consistent with pPXF
;; reconstruct the wavelength array in Anstroms, and compute rest-frame
;; values
;lambda0 = exp(dindgen(npix)*lstep_gal + l0_gal)
;if keyword_set(log10) then lambda0 = 10^(dindgen(npix)*lstep_gal + l0_gal)
;lambda0 = lambda0/exp(vstar/299792.458d)
;
;; array to hold k(lambda) values
;k = fltarr(n_elements(lambda0))           
;
;for i=0,n_elements(lambda0)-1 do begin
;     ; convert wavelength units from angstroms to micrometres
;     l = lambda0(i)/1e4                   
;     ; assign k values
;     if (l ge 0.63 and l le 2.2) then k(i) = 2.659*(-1.857+1.040/l)+4.05
;     if (l lt 0.63)              then k(i) = 2.659*(-2.156+1.509/l-0.198/l^2+0.011/l^3)+4.05
;     if (l gt 2.2)               then k(i) = 0.0
;endfor
;
;return,(10^(-0.4*ebv*k))
; this should be then multiplied by the spectrum flux array
;----------

; reconstruct the wavelength array in Anstroms, and compute rest-frame
; values
if keyword_set(log10) then lambda = 10^(dindgen(npix)*lstep_gal +l0_gal) $
	else lambda = exp(dindgen(npix)*lstep_gal + l0_gal)
lambda = lambda/exp(vstar/299792.458d)

; below copied from ppxf_reddening_curve
k1 = lambda*0
lam = 1e4/lambda ; Convert Angstrom to micrometres and take 1/lambda
rv = 4.05d ; C+00 equation (5)

w1 = where(lambda ge 6300d, m1, COMPLEMENT=w2, NCOMPLEMENT=m2)
; C+00 equation (3) but extrapolate for lam>2.2
if m1 gt 0 then k1[w1] = rv + 2.659d*(1.040d*lam[w1] - 1.857d)
; C+00 equation (4) but extrapolate for lam<0.12
if m2 gt 0 then k1[w2] = rv + $
    2.659d*(1.509d*lam[w2] - 0.198d*lam[w2]^2 + 0.011d*lam[w2]^3 - 2.156d)
fact = 10d^(-0.4d*ebv*(k1>0))  ; Calzetti+00 equation (2) with opposite sign

return, float(fact) ; The model spectrum has to be multiplied by this vector

end
