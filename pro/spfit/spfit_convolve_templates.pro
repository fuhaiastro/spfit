FUNCTION SPFIT_CONVOLVE_TEMPLATES, templates,kinstars,velscale
; HF: convolve LOSVD to the SSP templates and return the convolved SSPs

vel   = kinstars[0]/velscale   ; in pixels
sigma = kinstars[1]/velscale   ; in pixels
dx = ceil(abs(vel)+4*sigma)    ; Sample the Gaussian and GH at least to vel+4*sigma
x  = range(dx,-dx)             ; Evaluate the Gaussian using steps of 1 pixel.
w  = (x - vel)/sigma
w2 = w^2
; HF 20150807 modified to be consistent with pPXF - normalize before
; multiplying POLY
gauss = exp(-0.5d*w2)
losvd = gauss/total(gauss)

; Hermite polynomials as in van der Marel & Franx (1993).
; Coefficients are given e.g. in Appendix C of Cappellari et al. (2002)
nkins = n_elements(kinstars)
IF nkins GT 2 THEN BEGIN
    poly = 1d + kinstars[2]/Sqrt(3d)*(w*(2d*w2-3d)) $       ; H3
              + kinstars[3]/Sqrt(24d)*(w2*(4d*w2-12d)+3d)   ; H4
    IF nkins EQ 6 THEN $
        poly = poly + kinstars[4]/Sqrt(60d)*(w*(w2*(4d*w2-20d)+15d)) $      ; H5
                    + kinstars[5]/Sqrt(720d)*(w2*(w2*(8d*w2-60d)+90d)-15d)  ; H6
    losvd = losvd * poly
ENDIF 

s = size(templates)
ctemplates = dblarr(s[1],s[2])
IF s[0] EQ 2 THEN ntemp = s[2] ELSE ntemp = 1   ; Number of template spectra
FOR j=0,ntemp-1  DO BEGIN
        ctemplates[*,j] = convol(templates[*,j],losvd,/EDGE_TRUNCATE) ;
ENDFOR

return,ctemplates
END

