## To fit a single MaNGA datacube
```idl
fit_manga,7443L,1901L,$
	infile='manga-7443-1901-LOGCUBE.fits.gz',$
	drpfile='drpall-v2_7_1.fits',$
	mdegree=6,outdir='spfit/',mtpldir='spfit_mtpl/',$
	/nobinning,/verbose,/overwrite,/saveplot
```

## To fit a batch of datacubes 
- this mode assumes that the datacubes are downloaded to
`$MANGA_SPECTRO_REDUX/$MANGADRP_VER/`
- this mode uses N CPU-threads to fit N datacubes simultaneously using
  scripts under `pro/parallel` written by Alfred de Wijn.

```idl
; load drpall file for unique galaxy cubes
drpall = mrdfits('drpall-v2_7_1.fits',1)
; main galaxy sample: all Primary+ & full Secondary
plates = drpall.plate
ifus = long(drpall.ifudsgn)
; run 6 parallel processes
cmd = 'fit_manga'
extra = ',tag=''MPL-9'',mdegree=6,outdir=''spfit/'','+$
	'mtpldir=''spfit_mtpl/'',/nobinning,/quiet,/ignore_drp3qual'
; output files for manga_parallel to skip if alreay exist
output= 'spfit/'+strc(plates)+'-'+strc(ifus)+'.fits' 
manga_parallel,plates,ifus,cmd,extra,ncpu=6,outputs=output,/skip
```
