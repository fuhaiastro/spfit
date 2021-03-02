pro save_screen,outfile

	if strpos(outfile,'.png') ge 0 then begin 
		write_png,outfile,tvrd(/true) 
	endif else if strpos(outfile,'.jpg') ge 0 then begin
		write_jpeg,outfile,tvrd(/true),/true
	endif else begin
		print,'Error: unknown output file type'
	endelse
	
end
