pro	setx
;
	if strupcase(!d.name) eq 'PS' then device,/close
	set_plot,'x'
	;device,decompose=0
	multiplot,/default
	!p.font=-1
	!x.thick=1
	!y.thick=1
	!p.thick=1
	!p.charthick=1 
	!p.charsize=2
;
	end
