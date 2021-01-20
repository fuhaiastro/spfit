function emlerr,ampl_to_noise,sigma_pix,minerr=minerr,sigma=sigma
;+
; NAME
;	EMLERR
;
; PURPOSE
;	Returns fractional error in line flux, given the ampl_to_noise
;	ratio and the observed line dispersion in pixels.
;	
;	Explored parameter space: A/N = [2,20], Sigma = [1,20] pixels.
;	It would require extrapolation if the line parameter lies
;	outside. See mpfit_err.idlsrc for details. 
;
; SYNTAX
;	fractional error = EMLERR(ampl_to_noise, sigma_pix)
;
; INPUT
;	AMPL_TO_NOISE: scaler or array
;	SIGMA_PIX: scaler or array, must match the size of A/N
;	MINERR: minimum fractional error (default: 1%)
;	/sigma: returns fractional error of line width instead of flux
;
; OUTPUT
;	Fractional error in dF/F or dSigma/Sigma
;
; HISTORY
; 	written by HF, 20150730
;	added optional input MINERR, 20170914
;-
	; set minimum fractional error, default 1%
	if ~keyword_set(minerr) then minerr = 0.01

	; get polynomial parameters - kx 6x6 array
	if keyword_set(sigma) then $ ; dSigma/Sigma
		restore,'$SPFIT_DIR/pro/emlerr/mpfit_sigerr.idlsave' $
	else $ 			     ; dF/F
		restore,'$SPFIT_DIR/pro/emlerr/mpfit_err.idlsave'
	; evaluate at location
	df = poly2d(ampl_to_noise,sigma_pix,kx,/irregular)
	; return dF/F, the fractional error of line flux
	return,exp(df) > minerr 
	
end
