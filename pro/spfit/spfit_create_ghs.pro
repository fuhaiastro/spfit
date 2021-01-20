FUNCTION SPFIT_CREATE_GHS, X=x, PARS=pars,INT_DISP_PIX=int_disp_pix
;+
; Purpose
;	generate Gauss-Hermite series from input parameters	
;
; Input parameters:
;   	pars[0]: total area of the Gaussian before multiplying the Hermite
;		polynomials (unit: flux * dPix).
;   	pars[1]: centroid in pixel
;   	pars[2]: intrinsic disperison in pixel
;   	pars[3:4]: h3, h4
;   	
; 	INT_DISP_PIX: instrumental dispersion in pixels
;
; Output:
;   returns an array of a Gauss-Hermite series broadened 
;   by both intrinsic kinematics and instrument resolution. 
;   The output array has the same size as input X array
;
;-
npars=n_elements(pars)

sigma = sqrt(pars[2]^2 + int_disp_pix^2)
w = (X - pars[1])/sigma
w2 = w^2
gaus = EXP(-0.5d*w2)
if npars eq 5 then $ 
	hpoly = 1d + pars[3]/Sqrt(3d)*(w*(2d*w2-3d)) $   ; H3
          + pars[4]/Sqrt(24d)*(w2*(4d*w2-12d)+3d) $  ; H4
	else hpoly = 1d

; normalize to total integrated Gauss area and multiply by specified flux
y = gaus/total(gaus) * hpoly * pars[0]

return, y

END

