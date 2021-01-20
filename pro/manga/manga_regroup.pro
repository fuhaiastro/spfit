;+
; Name
;	manga_regroup
;
; Purpose
; 	Organize diagnostic diagrams into subfolders
; 	by creating symbolic links or copy files
;
; Syntax 
; 	MANGA_REGROUP, plateifu, outdir, indir= , /copy
; 
; Inputs 
;	plateifu [either an array of Plate-IFUs or
;                   a filename listing Plate-IFUs]
;         outdir [string giving the output directory]
;
;	if /copy is not set
;		indir gives the relative location of the input folder
;			from the output folder
;	if /copy is set
;		indir gives the relative location of the input folder 
;			from the program running folder
;
; Options
;         /copy [option to copy instead of link files]
;         /verbose [option to report # of files]
;
;-
pro manga_regroup,plateifu,outdir,indir=indir,copy=copy,verbose=verbose

if ( N_params() lt 2 ) or ~keyword_set(indir) then begin
   print,'Syntax: MANGA_REGROUP, plateifu, outdir, indir= , /copy, /verbose'
   return
endif  

; if plateifu is a file listing Plate-IFUs
if file_test(plateifu[0]) then $
	readcol,plateifu,obj,f='a',comment='#' else $ 
	obj = strtrim(plateifu,2)

; prepare outdir
if ~file_test(outdir) then spawn,'mkdir '+outdir else $
	spawn,'rm -f '+outdir+'/*.png'

; report # of files if /verbose is set
if keyword_set(verbose) then print,outdir+': ',n_elements(obj)

; generating symbolic links or copy files
; when copy=0, indir is the relative location of the input folder from the output folder
; when copy=1, indir is the relative location of the input folder from the running folder
; get full file path for input directory
; indirfull = file_search(indir,/full)
if ~keyword_set(copy) then $
	for kk = 0,n_elements(obj)-1 do $
		spawn,'ln -s '+indir+'/'+obj[kk]+'.png '+outdir+'/' $
else 	for kk = 0,n_elements(obj)-1 do $
		spawn,   'cp '+indir+'/'+obj[kk]+'.png '+outdir+'/'

end
