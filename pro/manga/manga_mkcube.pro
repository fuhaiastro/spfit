function manga_mkcube, binmap, val2d
;+
; NAME
;	MANGA_MKCUBE
;
; PURPOSE
;	Convert the 2D spectra (linearly arranged as a function of 
;	voronoi bin number) to a 3D spectra on sky. 
;
; SYNTAX
;	map = manga_mkcube(binmap,val2d)
;
; INPUT
;	binmap - [nx,ny] array giving bin index vs. x & y
; 	val2d - [nz, nbin] values to generate a cube for 
;
; OUTPUT
;	val3d - [nx,ny,nz] array to store the map, the same unit as val1d
;
; HISTORY
;	2018 Jun 27 - Written - HF
;
;-

nx = (size(binmap))[1]
ny = (size(binmap))[2]
nz = (size(val2d))[1]
nbin = (size(val2d))[2] 

val3d = fltarr(nx,ny,nz) ; [nx,ny,nz] cube
for ibin=0,nbin-1 do begin
	location = where(binmap eq ibin,ct)
	ind = ARRAY_INDICES(binmap, location)
	for i=0,ct-1 do val3d[ind[0,i],ind[1,i],*] = val2d[*,ibin]	
endfor
location = where(binmap eq -1,ct)
ind = ARRAY_INDICES(binmap, location)
for i=0,ct-1 do val3d[ind[0,i],ind[1,i],*] = !values.f_nan

return,val3d

end
