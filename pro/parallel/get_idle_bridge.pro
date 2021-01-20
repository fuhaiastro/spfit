function get_idle_bridge, bridges
; a wrapper that returns an IDL_IDLBridge object of a bridge that is currently idle.
	ncpus = n_elements(bridges)
	found = -1
	repeat begin
		for cpu=0,ncpus-1 do begin
			case (bridges[cpu])->status(error=errstr) of
				0: begin
					found = cpu
					;print,'found0 = ',cpu
				end
				2: begin
					found = cpu
					;print,'found2 = ',cpu
				end
				3: begin
					print, 'Error encountered: '+errstr
					;stop
				end
				4: begin
					print, 'Aborted execution: '+errstr
					;stop
				end
				else: ; do nothing
			endcase
			if found ne -1 then break
		endfor
		if found eq -1 then wait, 0.1 ; idle loop
	endrep until found ne -1
	bridge = (bridges[found])
	return, bridge
end
