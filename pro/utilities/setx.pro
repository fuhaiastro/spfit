pro	setx
;
	if strupcase(!d.name) eq 'PS' then device,/close
	set_plot,'x'
	; ensure write_png works properly
	device,decompose=0,retain=2
	multiplot,/default
	!p.font=-1
	!x.thick=1
	!y.thick=1
	!p.thick=1
	!p.charthick=1 
	!p.charsize=1
;
	end
