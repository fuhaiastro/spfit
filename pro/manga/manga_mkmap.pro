function manga_mkmap, binmap, val1d
;+
; NAME
;	MANGA_MKMAP
;
; PURPOSE
;	Convert the 1D measurements (linearly arranged as a function of 
;	voronoi bin number) to a 2D map on sky. Useful to illustrate the 
;	fitting results of datacubes (e.g., MaNGA).
;
; SYNTAX
;	map = manga_mkmap(binmap,val1d)
;
; INPUT
;	binmap - [nx,ny] array giving bin index vs. x & y
; 	val1d - [nbin] values to generate a map for 
;
; OUTPUT
;	val2d - [nx,ny] array to store the map, the same unit as val1d
;
; HISTORY
;	2015 May 22 - Written - HF
;
;-

val2d = binmap*0.0 ; [nx,ny] array to store the map
for ibin=0,n_elements(val1d)-1 do begin
	ind = where(binmap eq ibin)
	val2d[ind] = val1d[ibin]	
endfor
ind = where(binmap eq -1,ct)
val2d[ind] = !values.f_nan

return,val2d

end
