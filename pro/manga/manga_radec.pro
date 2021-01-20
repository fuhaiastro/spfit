function manga_radec,spfitfile
;+
; NAME
;	MANGA_RADEC
;
; PURPOSE
;	Calculate the on-sky coordinates (RA, Dec) of all bins in a
; 	MaNGA spfit file
;
; SYNTAX
;	radec = manga_mkcube('spfitfile')
;
; INPUT
;	spfitfile - a single string giving the filename including path
;
; OUTPUT
;	radec - [2,nbin] array to store the coordinates of all bins (degrees)
;
; HISTORY
;	2019 Apr 12 - Written - HF
;
;- 

on_error, 2

if ~file_test(spfitfile) then message,'Input file not found: '+spfitfile

vbin = mrdfits(spfitfile,2,/silent)
nbin = max(vbin.binmap)+1
radec = dblarr(2,nbin)

for i=0,nbin-1 do begin
	xy = array_indices(vbin.binmap,where(vbin.binmap eq i))
	xyad,vbin.hdr,mean(xy[0,*]),mean(xy[1,*]),ra,dec
	radec[*,i] = [ra,dec]
endfor

return,radec

end
