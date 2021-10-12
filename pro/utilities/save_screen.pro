;+
; For this to work properly, one must set 
;	device,decompose=0,retain=2
; before generating a plot. 
;
; retain = 2 : ensures that the window contents will be read properly
; decomposed = 0 : cause color values to be interpreted as indices into a color lookup table. 
;-
pro save_screen,outfile

	if strpos(outfile,'.png') ge 0 then begin 
		write_png,outfile,tvrd(/true) 
	endif else if strpos(outfile,'.jpg') ge 0 then begin
		write_jpeg,outfile,tvrd(/true),/true
	endif else begin
		print,'Error: unknown output file type'
	endelse
	
end
