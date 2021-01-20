function gaussian_filter1d, spec, sig 
;+
; NAME
; 	GAUSSIAN_FILTER1D 
;	based on ppxf_util.py by Michele Cappellari
;
; SYNTAX
;	convolved_spec = GAUSSIAN_FILTER1D(spec,sig)
; 
; INPUT
;	spec - flux array
;	sig  - sigma array for the Gaussian kernel in pixels
;
; PURPOSE
;    	Convolve a spectrum by a Gaussian with different sigma for every
;    	pixel, given by the vector "sigma" with the same size as "spec".
;
;	The required convolution using a variable kernel is time
;	consuming. This routine solves this problem by approximating the
;	convolution with a multiplication between shifted spectra and
;	part of the Gaussian (restricted to +/-3sigma), significantly
;	reducing the computing time.
;
;    	The first/last p=max(3*sig) pixels of the convolved spec are first 
; 	filled with zeros, and later replaced with the original input values.
;
; HISTORY
;	July 28, 2016 - Translated from the same function in
;		ppxf_util.py - Jacob Isbell & Hai Fu
;	May 3, 2019 - better commented and added code to restore the
;		first/last p pixels
;
;-

; if sig is provided as a constant, make it an array
;if n_elements(sig) eq 1 then sig = spec*0+sig

; forces zero sigmas to have 0.01 pixels
sig[where(sig lt 0.01)] = 0.01
p = ceil(MAX(sig*3)) 
m = 2*p+1 ; kernal size

n = n_elements(spec)
a = FLTARR(m,n)
; for Gaussian kernel
x2 = (FINDGEN(m)-p)^2 ; make an array from -p^2 to p^2
gau = FLTARR(m,n) ; array to store Gaussian kernel
; translate spec into an (2p+1)xn array
; loop over the smaller dimension
; here the first/last p pixels are filled with zeros
for j=0,m-1 do begin 
	a[j,p:n-1-p] = spec[j:j+n-m] ; (n-1-p)-p+1 = n-2p = n-m+1 elements 
	; note that the input spec is copied and shifted by 1 pixel each row,
	; so that convolution can be done w/ a multiplication
	gau[j,*]  = exp(-1*x2[j]/(2*sig^2)) 
	; evaluate the Gaussian kernels with variable width at each row
	; each row represent an offset from the center (x2[j])
	; note that this Gaussian has not been normalized by sqrt(2pi)*sig
endfor
; normalize kernel
gau = gau / (total(gau,1) ## (fltarr(m)+1))
; gaussian area = sqrt(2*!PI)*sig, so we may divide that out instead

; multiply w/ kernel and sum
; this is equivalent to a convolution
conv_spectrum = total(a*gau,1)

; restore the first/last p pixels with original spec
conv_spectrum[0:p-1] = spec[0:p-1]
conv_spectrum[n-p:*] = spec[n-p:*]

return, conv_spectrum

end
