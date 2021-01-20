PRO mask_polygon, image, ds9reg, value=value, nreg=nreg
;+
; PURPOSE:
;	Generate a mask array the same size as the input image array
;	using polygon regions defined in ds9-format polygon files.
; 	Based on http://www.idlcoyote.com/tips/point_in_polygon.html
;
; INPUT:
;	image - 2D image array to be masked
;
;	ds9reg - string, name of the ds9 region file containing the polygon 
;		information. This file can be generated with ds9. 
;    		Or, any text file that contains lines starting with 'polygon('
;    		(no space or any other characters in front) and end with ')' 
;		can be used.  Lines not starting with 'polygon(' will be igored.  
;		The format should be:
;			
;			polygon(x1,y1,x2,y2,x3,y3,x4,y4...,xn,yn)
;		
;		for polygons with n points. The polygon loop should not be 
;		closed, i.e., xn not equal to x1, same for y.  There should
;		be no space at all in the line.  xn and yn can be integers or
;		floats.  They should be IMAGE COORDINATES.  RA/Dec will not be
;		accepted.  The points CAN be outside the image.  Multiple polygons
;		can be masked at the same time in one polygon file.
;
;	value - starting value to fill in the masked regions. Default is zero, if
;		not provided. All pixels inside a polygon will be masked 
;		with the same value. The value increments for every polygon.
;
; OUTPUT:
;	image - masked 2D image array
;
;	nreg  - integer, giving the number of polygons
;
; HISTORY:
;	Use Coyote's inside.pro to figure out if a pixel is inside or outside a
;	polygon - Aug 31, 2017 - HF
;	Following the recommendation of Rui, now use IDLanROI
;	ContainsPoints method - Aug 31, 2017 
;	Include points on edges and vertexes - Oct 15, 2017
;	Small regions may not enclose any pixels, added safe-guard to
;	prevent no-pixel-inside-mask cases - Oct 20, 2017
;
;-

; define the base value of the mask
IF n_elements(value) EQ 0 THEN value = 0.0
image[*,*] = value

; generate X & Y arrays for pixel centers
; Here we use the DS9 definition because the polygons are defined for DS9
; DS9 defines the center of the 1st pixel as 1.0
imsize = size(image,/dimen)
xx = (findgen(imsize[0])+1) # (fltarr(imsize[1])+1)
yy = (fltarr(imsize[0])+1) # (findgen(imsize[1])+1)

openr, ds9, ds9reg, /get_lun
line=''
nreg=0

WHILE NOT EOF(ds9) DO BEGIN   	
	readf, ds9, line
	; remove all blank spaces
	line = repstr(line,' ','')
	; skip lines that don't begin w/ polygon(
	IF strmid(line, 0, 8) NE 'polygon(' THEN continue   
	
	; increment the filling value by one	
	value += 1
	; increment count
	nreg += 1
	
	; extract polygon coordinates
	coors = strmid(line,8,strlen(line)-9)
	coors = float(strsplit(coors,',',/extract))
	npoint = n_elements(coors)/2
	ind = indgen(npoint)
	px = coors[ind*2]
	py = coors[ind*2+1]
	
	; call Coyote inside.pro
	;is_inside = inside(xx,yy,px,py)
	;image[where(is_inside)] = value

	; use IDLanROI ContainsPoints method
	; 0 = Exterior. The point lies strictly out of bounds of the ROI
	; 1 = Interior. The point lies strictly inside the bounds of the ROI
	; 2 = On edge. The point lies on an edge of the ROI boundary
	; 3 = On vertex. The point matches a vertex of the ROI
	object = Obj_New('IDLanROI', px, py)
	idx = where(object->ContainsPoints(xx, yy) ge 1,ct)
	; safeguard
	if ct gt 0 then begin
		image[idx] = value 
	endif else begin
		mpx = mean(px)
		mpy = mean(py)
		if mpx gt 0 and mpx lt imsize[0] and $
		mpy gt 0 and mpy lt imsize[1] then image[mpx-1,mpy-1] = value
	endelse
	Obj_Destroy, object
ENDWHILE

free_lun, ds9

END
