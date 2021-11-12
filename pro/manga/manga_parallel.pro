pro manga_parallel, plates, ifus, cmd, extra, outputs=outputs, $
	skip=skip, quiet=quiet, sequential=sequential, ncpu=ncpu
;+
; PURPOSE
; 	Run parallel sessions with IDL bridge
; 	* only works when cmd and extra are scalers and plates & ifus are
; 	integer arrays of the same number of elements
;	* outputs must have the same number of elements as
;	plates and ifus
; 	* if /skip is set, then if outputs[i] file exist, skip processing i.
;	* outputs array is for manga_parallel to determine if the output
;	file already exists, it's not used to set the output name for
;	the script wrapped within (e.g., manga_fit)
;	* if /sequential is set, then start non-parallel sessions
;	(useful for debugging purposes)
;
; HISTORY:
;	written by H. Fu based on a number of bridge scripts
;-

if keyword_set(quiet) then verb = 0 else verb = 1
if keyword_set(skip) and n_elements(outputs) eq n_elements(plates) $
	then toskip = 1 else toskip = 0

if verb then begin
	tic 
	print,' # of galaxy cubes to fit: ',n_elements(plates)
	print,cmd+extra
endif

if ~keyword_set(sequential) then begin
	; build bridges
	if keyword_set(ncpu) then ncpus=ncpu else ncpus = 12 < n_elements(plates)
	nthreads = 1
	if verb then print,'use ', ncpus, ' CPUs each with',nthreads, ' threads '
	bridges = build_bridges(ncpus,nthreads)
	; looping through all datacubes 
	for i=0,n_elements(plates)-1 do begin
		; if exist then skip
		if toskip then if file_test(outputs[i]) then continue
		; grab an idle bridge
		bridge = get_idle_bridge(bridges)
		; once got an idle bridge show progress
		if verb then print,i+1,plates[i],ifus[i],f='(i,1x,i,1x,i)'
		; define input parameter
		bridge->setvar, 'plate', plates[i]
		bridge->setvar, 'ifu', ifus[i]
		; execute command asynchronously
		bridge->execute, /nowait, $
			cmd+', plate, ifu'+extra 
	endfor
	; blocks until all bridges are idle 
	barrier_bridges, bridges
	burn_bridges,bridges
endif else begin
	; non-parallel sessions
	for i=0,n_elements(plates)-1 do begin
		; if exist then skip
		if toskip then if file_test(outputs[i]) then continue
		; show progress
		if verb then print,i,plates[i],ifus[i],f='(i,1x,i,1x,i)'
		xx = execute(cmd+', plates[i], ifus[i]'+extra)
	endfor
endelse

; print ellapsed time
if verb then toc

end

