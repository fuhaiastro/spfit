;+
; NAME
;	MANGA_SAVEPLOT
;
; PURPOSE
;	A wrapper to call SHOW_SPEC_FIT repeatedly for MaNGA.
;	Given a plate number and IFU design, save PNG/PS summary 
; 	figures for all spaxels in the current directory, provided
; 	that the fitting results exist 
;
; INPUT
; 	plate, ifu - [Integers] MaNGA plate number and IFU design
;
; OPTIONAL INPUT
;	indir - input directory (default: 'spfit')
;	outdir - output directory (default = indir)
;	step - integer, produce one figure every step bins (default: 1)
; 	/ps - save summary figure as a PS file. 
;		if not set and /saveplot is set, save a PNG file
;	/BLR - set to plot the fitting results that include broad lines
;
; OUTPUT
;	Plate-IFUdesign[-br]/#.png or ps - 
;		summary figure illustrating the quality of the fit
;
; HISTORY
;	2015/7/20 - Written - HF
;	2015/8/12 - added /BLR option
; 
;-
pro manga_saveplot, plate, ifu, INDIR=indir, OUTDIR=outdir, $
	STEP=step, PS=ps, BLR=blr

if ~keyword_set(step) then step = 1
if ~keyword_set(indir) then indir = '$MANGA_DIR/hfdap/spfit'
if ~keyword_set(outdir) then outdir = indir
if ~file_test(outdir) then spawn,'mkdir '+outdir

; input directory - where the HFDAP products are saved
if ~keyword_set(BLR) then basename = strc(plate)+'-'+strc(ifu) $
	else basename = strc(plate)+'-'+strc(ifu)+'-br'
vbin_file = file_search(indir+'/'+basename+'.fits',count=ct)
if ct eq 1 then gpars = mrdfits(vbin_file,1,/silent) else return

; make a directory to save screenshots
plotdir = outdir+'/'+basename 
if ~file_test(plotdir) then spawn,'mkdir '+plotdir else spawn,'rm -f '+plotdir+'/*'

; start looping through all voronoi bins
for ibin = 0, n_elements(gpars)-1, step do begin
	fit_results = gpars[ibin]	
	; save fitting results
	plotfile = plotdir+'/'+string(ibin,f='(i04)')
	mylegend = basename+' bin#:'+strc(ibin)
	if keyword_set(ps) then begin
		; generate PS file to exam the quality of the fit
		show_spec_fit,fit_results,$
			mylegend=mylegend,/ps,outfile=plotfile+'.ps'
	endif else begin
		if ibin eq 0 then begin
			device,decomp=0 
   			window,0,xs=600*2.0,ys=400*2.0
   		endif
		show_spec_fit,fit_results,$
			mylegend=mylegend
		save_screen,plotfile+'.png'
	endelse
endfor


end
