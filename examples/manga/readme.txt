;------------------------------;
; 1. fit a single datacube
;------------------------------;
fit_manga,7443L,1901L,$
	infile='manga-7443-1901-LOGCUBE.fits.gz',$
	drpfile='drpall-v2_7_1.fits',$
	mdegree=6,outdir='spfit/',mtpldir='spfit_mtpl/',$
	/nobinning,/verbose,/overwrite,/saveplot

;------------------------------;
; 2. fit a batch in parallel 
; Note: this mode requires the datacubes
;	to be stored under
;	$MANGA_SPECTRO_REDUX/$MANGADRP_VER/
;------------------------------;
; load drpall file for unique galaxy cubes
drpall = mrdfits('drpall-v2_7_1.fits',1)
; main galaxy sample: all Primary+ & full Secondary
plates = drpall.plate
ifus = long(drpall.ifudsgn)
; run 6 parallel processes
cmd = 'fit_manga'
extra = ',tag=''MPL-9'',mdegree=6,outdir=''spfit/'','+$
	'mtpldir=''spfit_mtpl/'',/nobinning,/quiet,/ignore_drp3qual'
; output files for manga_parallel to skip if exist
output= 'spfit/'+strc(plates)+'-'+strc(ifus)+'.fits' 
manga_parallel,plates,ifus,cmd,extra,ncpu=6,outputs=output,/skip

