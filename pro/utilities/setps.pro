pro setps,psfile,xsize,ysize,font=font,pthick=pthick,charsize=charsize

	if n_params() lt 1 then begin
  		print,'syntax: setps, psfile, xsize, ysize [, font=, pthick= ]'
		print,'commonly used fonts: Helvetica, Palatino-Roman'
		print,'if font is not defined, use default vector fonts'
  		retall
	endif

	set_plot,'ps' ; set graphics device
	if keyword_set(font) then begin
		!p.font=0 ; switch from vector to device fonts
		; choose a PostScript Font	
		device,xsize=xsize,ysize=ysize,bit=8,file=psfile,$
		/color,/encapsulate,set_font=font,isolatin1=1
		; /ISOLATIN1
		; Set this keyword to use Adobe ISO Latin 1 font encoding 
		; with any font that supports such coding. Use of this keyword 
		; allows access to many commonly-used foreign characters.
	endif else begin
		; if font is not set, use default vector fonts
		!p.font = -1
		device,xsize=xsize,ysize=ysize,bit=8,file=psfile,$
		/color,/encapsulate
	endelse
	
	if ~keyword_set(pthick) then pthick=4
	!x.thick=pthick
	!y.thick=pthick
	!p.thick=pthick
	if ~keyword_set(charsize) then charsize=1.2
	!p.charsize=charsize
end
