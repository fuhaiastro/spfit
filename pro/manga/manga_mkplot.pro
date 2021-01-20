pro manga_mkplot, map_in, zr=zr, colorbar=colorbar, ztit=ztit, $
	ztickname=ztickname, zminor=zminor, ztickint=ztickint, zsize=zsize, $
	c_map = c_map, c_levels=c_levels, $ ; c_colors=c_colors,$
	c2_map = c2_map, c2_levels=c2_levels, $ 
	sauron=sauron, position=position, reverse_ct=reverse_ct,$
	_EXTRA=EXTRA
;+
; NAME
;	MANGA_MKPLOT
;
; PURPOSE
;	Plot 2D maps with a color bar. Useful to illustrate the 
;	fitting results of datacubes (e.g., MaNGA). Has the option to
;	overlay contour plots if c_map is set
;
; SYNTAX
;	manga_mkplot, map, zr=[0,10], /colorbar, ztit='[O III]/Hbeta',$
;		c_map=c_map, c_levels=c_levels, sauron=sauron
;
; INPUT
;	map - [nx,ny] array
;	sauron - colormap index
;
; HISTORY
;	2015 May 22 - Written - HF
;	2016 Feb 26 - added PS support 
;	2017 Aug 26 - enhanced sauron keyword to support other IDL color tables
;	2017 Oct 17 - modifies map values to avoid color index 0 and 255
;	2018 Nov 12 - revert scaling ztitle font by 1.3x
;
;-

; set colormap
if ~keyword_set(sauron) then sauron = 0 ; if not set, use B&W Linear
; if set to 1 then use sauron_colormap, any other numbers use IDL
; colormaps: https://www.harrisgeospatial.com/docs/LoadingDefaultColorTables.html
if sauron eq 1 then sauron_colormap else loadct,sauron
; reverse colormap
; http://www.idlcoyote.com/color_tips/reverse_ct.html
if keyword_set(reverse_ct) then begin
	TVLCT, r, g, b, /Get
	TVLCT, Reverse(r), Reverse(g), Reverse(b)
endif
; set charsize for colorbar
if ~keyword_set(zsize) then zsize = 1.2 ; tickmark size
if ~keyword_set(tsize) then tsize = zsize ; title size

map = map_in 	; save map to another array to avoid overwriting 
		; the original input

; avoid antialiasing for small maps by enlarging them
minsize = 100
dim = (size(map))[1:2]
if min(dim) lt minsize then begin
	f = round(minsize/min(dim))
	map = rebin(map,dim[0]*f,dim[1]*f,/sample)
	if keyword_set(c_map) then $
		c_map = rebin(c_map,dim[0]*f,dim[1]*f) 
	if keyword_set(c2_map) then $
		c2_map = rebin(c2_map,dim[0]*f,dim[1]*f) 
endif

; plotting positions
if keyword_set(colorbar) then begin
	r_xpos = position[2]-position[0]
	r_ypos = position[3]-position[1]
	pos = position-[0,0,0,0.2*r_ypos]
endif else $
	pos = position 

; for PS plot - use white color for NaNs or below minimum value
if !d.name eq 'PS' then begin
	; overwrite color table
	mincolor = cgColor('White',0)
	maxcolor = cgColor('Black',!D.Table_Size-1)
	;; use mincolor for lower values
	;ind = where(map lt zr[0],ct)
	;if ct gt 0 then map[ind] = zr[0]

	; aviod white color at minimum value
	s = where(map-zr[0] lt (zr[1]-zr[0])*1./256.,ct1)
	if ct1 gt 0 then map[s] = zr[0]+(zr[1]-zr[0])*1./256
	; aviod black color at maximum value
	s = where(map-zr[0] gt (zr[1]-zr[0])*255./256.,ct1)
	if ct1 gt 0 then map[s] = zr[0]+(zr[1]-zr[0])*255./256
	
	; use mincolor for NaNs
	ind = where(~finite(map),ct)
	if ct gt 0 then map[ind] = zr[0]
endif

; show image
tvimage,bytscl(map,min=zr[0],max=zr[1],/nan,top=255),/noint,position=pos
; plot contours
if keyword_set(c_map) then begin
	if !d.name eq 'PS' then $
		c_colors=cgColor('black',!D.table_size-1) $
		else c_colors=cgColor('white',!D.table_size-1)
	contour,c_map,levels=c_levels,/noerase,c_colors=c_colors,C_THICK=1,$
		position=pos,xs=4+1,ys=4+1
endif
; plot frame
plot,[1,1],[1,1],/nodata,xs=1,ys=1,pos=pos,color=maxcolor,_extra=extra 

; plot colorbar
if keyword_set(colorbar) then begin
	cpos=[pos[0]+0.1*r_xpos,pos[3],pos[2]-0.1*r_xpos,pos[3]+0.07*r_ypos]
	cbar,vmin=zr[0],vmax=zr[1],pos=cpos,color=maxcolor,$
	   	xticklen=0.29,/top,/horizontal,xran=zr,$
		cmax=!d.table_size-1,cmin=1,$ ; avoid the first color
		charsize=zsize,xtickname=ztickname,xminor=zminor,xtickint=ztickint
	xyouts,(pos[0]+pos[2])/2,pos[3]+0.14*r_ypos*1.1,$
		ztit,/normal,align=0.5,charsize=tsize,color=maxcolor
endif

; redefining this color will mess up the color map
; so do this after drawing the colorbar
if keyword_set(c2_map) then begin
	c_colors=cgColor('blue',!D.table_size-2)
	contour,c2_map,levels=c2_levels,/noerase,c_colors=c_colors,C_THICK=1,$
		position=pos,xs=4+1,ys=4+1
endif


end
