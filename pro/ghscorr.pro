function ghscorr, h3, h4 
;+
; Purpose
;	compute flux correction factor from Gaussian to 
;	Gauss-Hermite series from input parameters. 
;
;	Because H3 term is symmetric w.r.t. the centroid, 
;	it does not affect the correction factor (Corr).
;	Corr is controlled mostly by H4:
;		for H4 = [-0.3,0.3], Corr = [0.82, 1.18]
;	Corr only weakly depends on the line width because 
;	of finite sampling, which we donot consider here
;
; Input
;   	h3, h4
;
; Output
;	the correction factor := Flux(GHS)/Flux(Gauss)
;
; History
;	20150812 - written by HF
;-

w = range(-20,20,201) ; W := (x-centroid)/sigma
w2 = w^2
hpoly = 1d + h3/Sqrt(3d)*(w*(2d*w2-3d)) $   ; H3
  + h4/Sqrt(24d)*(w2*(4d*w2-12d)+3d)   ; H4
gaus = EXP(-0.5d*w2)

return,total(gaus*hpoly)/total(gaus)

end
