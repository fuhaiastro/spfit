pro burn_bridges, bridges
; destroy the bridges after we're done with the calculations
	ncpus = n_elements(bridges)
	for cpu=0,ncpus-1 do $
		obj_destroy, bridges[cpu]
end. 
