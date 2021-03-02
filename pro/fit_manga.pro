pro fit_manga, plate, ifu, INFILE=infile, DRPFILE=drpfile, $
	SSPLIB=ssplib, LINEFILE=linefile, $
	OUTDIR=outdir, TAG=tag, RDIR=rdir, VERBOSE=verbose, OVERWRITE=overwrite, $
	targetSN=targetSN, SAVEPLOT=saveplot, PS=ps,$
	ignore_drp3qual=ignore_drp3qual, STEP=step, $
	NOBINNING=nobinning, DS9REG=ds9reg, regnum=regnum, maskdir=maskdir,$
	BLR=blr, GANFIT=ganfit, $
	delete=delete,mtpldir=mtpldir,_EXTRA=extra
;+
; NAME
;	FIT_MANGA
;
; PURPOSE
;	Top level wrapper to fit MaNGA datacubes with SPFIT or GANFIT
;	Note that 2x2 binning to 1-arcsec spaxels is performed before any fitting
;
; INPUT
; 	plate, ifu - [Integers] MaNGA plate number and IFU design
;
; OPTIONAL INPUT
;	infile - input file name, including path
;		Default: rdir+tag+'/'+strc(plate)+'/stack/manga-'+
;			strc(plate)+'-'+strc(ifu)+'-LOGCUBE.fits*'
;	drpfile - input DRPALL file
;		Default: rdir+tag+'/drpall*.fits'
; 	ssplib - string defines the SSP library to be used. 
;		Curent options include 'miuscat', 'miuscat-thin', and 'm11_stelib', 
;		both use Kroupa IMF. (default: miuscat-thin)
;	linefile - emission line setup file
;	outdir - output directory (default: ./spfit/)
;	tag - DRP tag number (default: $MANGADRP_VER)
;	rdir - DRP redux directory (default: $MANGA_SPECTRO_REDUX)
;	/verbose - print information at each stage and report total elapsed time
;	/overwrite - if set, overwrite previous fitting results;
;		- Will not overwrite matched SSP templates even if set;
;		- If not set and fitting result exists, then quit program 
;	/saveplot - make summary figure and save it as either PNG or PS file
; 	/ps - save summary figure as a PS file. 
;		if not set and /saveplot is set, save a PNG file
;	/ignore_drp3qual - if set, proceed with fitting even if DRP3QUAL
;		indicates critical failure.
;	step - integer, giving the number of bins to skip in each FOR
;		loop iteration. Set to a high number for quick testing.
;		Default: 1 (fitting every bin)
;	/nobinning - if set, then do not perform adaptive spatial binning
;	DS9reg - the directory that keeps the DS9 polygon region file
;		regions are defined in image coordinates *after 2x2 binning*
;		e.g., 'reg1arc/' (Note: must include backslash at the end)
;	/BLR - set to fit with broad emission lines
;	/GANFIT - set to use GANFIT instead of SPFIT (Default is SPFIT)
;	/DELETE - set to remove logcube file to save space
;	MTPLDIR - directory that saves the matched SSPLIB file for each
;		plate-ifu (default = $MANGA_DIR/hfdap/TAG/spfit_mtpl/)
;	REGNUM - integer scaler giving the region number. When set, it
; 		will add the number to the plateifu for the output file.
;		It requires that (1) the DS9reg folder is set, and (2) the
;		corresponding region file exists (e.g., 8248-1902-3.reg) 
;	MASKDIR - the directory that keeps the FITS catalog for individual
;		plateifu; it must contain RA and Dec tags for the
;		objects to be masked out.
;	_EXTRA - SPFIT/GANFIT Keywords
;
; OUTPUT
;	Plate-IFUdesign.fits or _br.fits (if /BLR) - 
;		ext 1: best-fit pars
;		ext 2: RA, Dec, binmap, S/N maps, line names, FITS header w/ WCS info
;		ext 3: drpall struct
;		ext 4: spectral data and models
;	Plate-IFUdesign_ssplib.fits - spectral resolution matched 
;		log10-rebinned SSP templates
;	Plate-IFUdesign/#.png or ps (if /saveplot)- 
;		summary figure illustrating the quality of the fit
;
; HISTORY
;	2015/5/15 HF - Written
;	2015/7/31-8/5 HF - made quite a number of changes 
;		- instead of passing linefile for gandalf_fit to read,
;		  now we construct emission_setup structure from
;		  linefile when processing the 1st bin and passing it to
;		  subsequent bins.
;		- added ability to fit broad lines w/ Gaussian and test whether 
;		  the broad line is detected in the 1st bin. If detected, 
;		  hold the kinematics of the broad lines for the entire cube.
;		- Added /BLR option to decide if to fit broad lines.
;		  Useful when there is a list of objs that are known to
;		  have broad lines (Seyferts 1, 1.5, or 1.9).
;		- Narrow-line-only fits are done for all spaxels even if
;		  broad lines are detected, because sometimes including
;		  the broad lines can mess up the narrow line fits. 
;		- computes error array using Kyle's formula and
;		  modified manga_voronoi_2d_binning to bin accordingly.
;		- requires at least 80% pixels are good to be included in the SNR
;		  map, and therefore in the fitting.
;		- SNR map now computed using 97.5%tile of a large
;		  wavelength range (r.f., 3700 - 7000 A).
;		- masks bright stars by inquiry SDSS photo catalog,
;		  saved under cats/sdss/.
;		- SSP templates wavelengths converted from AIR to
;		  VACUUM when constructing the template structures.
;	2015/8/7  - converted from gandalf_manga to spfit_manga
;	2015/8/10 - unified gandalf_manga and spfit_manga to fit_manga
;		- added masking first 40 pixels, note that D.
;		  Law emailed that "the actual wavelength range of the cubes has 
;		  already been trimmed from the wavelength range covered
;		  by the detectors"
;		- calculate S/N raw map using the total S/N at Ha-NII region
;		- updated default targetSN to 20 
;		- Implement COV matrix in voronoi binning 
;		- Use COV matrix to calculate error vectors
;	2015/8/16 - for /BLR option, use 4th order Legendre polynomial
;		to approximate the AGN continuum; add the best-fit poly
;		to the templates array for subsequent bins. This
;		significantly improved the fit for the broad-line AGNs, 
;		reducing the chi^2/dof by ~2-3x. 
;		- changed masking half-width from 20 to 10 A for skylines.
;	2016/3/8 - Added DS9reg option to bin spaxels based on Polygons
;		- for /BLR option, disabled additive 4th order
;		polynomial to approximate AGN continuum
;	2016/6/29 - Utilized the sparse covariance matrix included in
;		the datacube as of MPL-5
;	2016/10/14 - Saved spec-resolution matched SSP libs in the same
;		output directory instead of a separate folder
;		- because the spec-res matching speed has improved a
;		lot by using gaussian_filter1d, we no longer need to pre-compute 
;		the matched SSPs
;		- removed manga_matchsres.pro from Pro folder
;	2016/10/14 - Retired foreground star masking using cats/sdss/
;		because MPL-5 now correctly set bitmask for regions 
;		within 4"-radius of stars
;	2016/10/17 - Saves badfrac array in the output: the fraction of
;		masked out pixels in a spectrum at (x,y)
;	2016/10/21 - Ignores foreground star masks
;		- saves middle-channel bitmask in vbin (ext [2])
;	2017/09/11 - replaced "message" statement in overwrite block
;		with a GOTO command to avoid stopping parent program
;		- redefined ds9reg as the folder name instead of region file name
;		this makes it easier to run the code in parallel mode
;		- added SNRRAW > 4 constraint for DS9-region-defined
;		bins to avoid including spaxels outside of IFU
;	2018/04/18 - added /DELETE option and a safe-guard for z < 0 cases
;	2018/07/27 - aded MTPLDIR parameter to look for SRES-matched SSP library
;		- added REGNUM parameter to be compatible with Josh's region names
; 	2018/10/15 - Josh reported a bug -- if there are DS9 regions
;		containing no valid pixels, the binmap will be messed
;		up. Fixed this by separating bin number from region
;		number and by keeping only valid bins after looping
;		through all region numbers
;	2019/07/08 - added RA & Dec tags in vbin (ext [2]) to store
;		the mean coordinates of each bin
;		- changed REGNUM parameter from string to integer and
;		added REGNUM to the output filename
;		- TODO: need to think about how to run manga_parallel with
;		the REGNUM option
;	2019/07/17 - When DS9reg is set, do not remove bins w/ S/N=0
;		before generating the 2D binmap and 2D SNR map, and do not 
;		sort bins using distance to ObjRA/Dec. This way 
;		the spectra in the output file will match the sequence of 
;		the input region file, *if all regions lie within the
;		square that snug fits the hexgon IFU. 
;		Regions outside of the square will still be skipped.
;		Low S/N bins will be skipped in
;		the SPFIT FOR loop and will have a null gpars structure.
; 	2019/07/17 - now set three default targetSN, depending on the
;		binning scheme
; 	2019/07/25 - added new NREG keywork in mask_polygon to count the
;		number of input polygons and use NREG to set NBIN (instead
;		of using max(segmap). This way, even if some polygons lie
;		outside of the image, the output still matches the input 
;		in terms of # of elements. Such outside polygons are
;		treated the same way as low S/N regions -- they are not
;		fit but will have a zero fit_results structure in return.
;		The only difference is that these outside polygons will
;		have RA/Dec = -99 while low S/N regions will still have
;		valid RA/Dec.
;	2019/07/26 - added MASKDIR option
;		- corrected how the new CRPIX is calculated (=(CRPIX1+0.5)/2)
;	2019/12/12 - save model spectra in a separate output file (*_spec.fits)
;	2021/01/19 - added optional keywords INFILE & DRPFILE 
;-

; error action
ON_ERROR, 2

; set up default parameters
IF ~keyword_set(linefile) then $
	linefile = getenv('SPFIT_DIR')+'/pro/emission_lines_setup.txt'
if ~keyword_set(ssplib) then ssplib = 'miuscat-thin' ; default SSP templates
if ~keyword_set(tag) then tag = getenv('MANGADRP_VER')	; DRP tag
if ~keyword_set(ganfit) then cmd = 'spfit' else cmd = 'ganfit'
; default output directory
if ~keyword_set(outdir) then begin
	if keyword_set(ds9reg) then outdir = ds9reg $
		else outdir = './'+cmd 
endif
if ~file_test(outdir) then spawn,'mkdir '+outdir
; default input directory
if ~keyword_set(rdir) then rdir = getenv('MANGA_SPECTRO_REDUX') ; DRP output dir
; default increment: fit every bin
if ~keyword_set(step) then step = 1
; default matched SSP templates directory
if ~keyword_set(mtpldir) then mtpldir = getenv('MANGA_DIR')+'/hfdap/'+tag+'/spfit_mtpl/' 

; set default target S/N
; skip bin if SNR < targetSN, below which we don't trust the fit anymore
if ~keyword_set(targetSN) then begin
	if keyword_set(ds9reg) then begin
		targetSN = 3.0
	endif else if keyword_set(nobinning) then begin
		targetSN = 12.0
	endif else begin
		targetSN = 17.0
	endelse
endif

; to avoid checking repeatedly if a keyword is set
if keyword_set(verbose) then verb = 1 else verb = 0

; start counting elapsed time
if verb then tic

; base filename to be used for output files and to identify the datacube
basename = strc(plate)+'-'+strc(ifu)

; catalog of observed datacubes
if ~keyword_set(drpfile) then drpfile = rdir+tag+'/drpall*.fits'
if ~file_test(drpfile) then message,basename+' File not found! '+drpfile
drpall = mrdfits(drpfile,1,/silent)
; find catalog info for the input plate & IFU
i = where(drpall.plate eq plate and strtrim(drpall.ifudsgn,2) eq strc(ifu),ct)
if ct eq 0 then message,basename+' not found!'
if ct gt 1 then message,basename+' not unique!'
; datacube file on disk?
if ~keyword_set(infile) then $ 
infile = rdir+tag+'/'+strc(plate)+'/stack/manga-'+strc(plate)+'-'+strc(ifu)+'-LOGCUBE.fits*'
if ~file_test(infile) then message,basename+' File not found! '+infile
infile = (file_search(infile))[0]

; from drpall, pull out redshift info, MANGA-ID, Galactic extinction
drpall = drpall[i]
objz = float(drpall.nsa_z)
if drpall.nsa_z le 0 and drpall.nsa_zdist gt 0 then objz = float(drpall.nsa_zdist)
if drpall.nsa_z le 0 and drpall.z gt 0 then objz = float(drpall.z)
if objz le 0 then message,'Invalid redshift, z = '+strc(objz)+', skip '+basename
if ~keyword_set(ebv_gal) then ebv_gal = float(drpall.ebvgal) ; foreground Galactic extinction
; DRP quality flag
if verb then print,'DRP3QUAL:',sdss_flagname('MANGA_DRP3QUAL',drpall.drp3qual)
if ~keyword_set(ignore_drp3qual) and $
	(drpall.drp3qual AND sdss_flagval('MANGA_DRP3QUAL', 'CRITICAL')) NE 0 then $
	message,'DRP flagged Critical Failure, skip '+basename

; check if output file already exist
if keyword_set(ds9reg) and keyword_set(regnum) then $
	outfile = outdir+basename+'-'+strc(regnum)+'.fits' $
	else outfile = outdir+basename+'.fits'
if keyword_set(BLR) then outfile = repstr(outfile,'.fits','-br.fits') 

if file_test(outfile) then begin ; if outfile already exist
	if ~keyword_set(overwrite) then begin
		print,'Already done! Skipping '+basename
		goto,theend
	endif else begin
		if verb then print,'Overwriting '+basename
	endelse
endif

; make a directory to save screenshots
if keyword_set(saveplot) then begin
	plotdir = repstr(outfile,'.fits','')
	if ~file_test(plotdir) then spawn,'mkdir '+plotdir $
	else spawn,'rm -f '+plotdir+'/*'
endif

; Read MaNGA datacube, which is already in log10-lambda steps
; https://trac.sdss.org/wiki/MANGA/TRM/TRM_MPL-3/datamodel
; make sure all input arrays are float arrays - to avoid issues w/ pPXF
if verb then print,'--> reading in datacube ' + infile
flux = mrdfits(infile,'FLUX',hdr,/silent) ; float, log10-rebinned flux cube [1E-17 erg/s/cm^2/Ang/spaxel]
ivar  = mrdfits(infile,'IVAR',/silent) ; float, log10-rebinned inverse variance cube
mask = mrdfits(infile,'MASK',/silent) ; long, corresponding MANGA_DRP3PIXMASK cube
wave = float(mrdfits(infile,'WAVE',/silent)) ; double->float, wavelength in Angstrom
sres = float(mrdfits(infile,'SPECRES',/silent)) ; double->float, spec resolution R = Wave/FWHM
; Note that the covariance matrix is not provided in DR13 and DR14
covs = mrdfits(infile,'RCORREL',hcovs,/silent) ; 'r' band sparse correlation matrix

; check Galactic extinction value
if verb then print,'--> E(B-V) Galactic: ', ebv_gal, sxpar(hdr,'EBVGAL')

; get the dimensions of the cube
dim = (size(flux))[1:3]

; evaluate covariance matrix for original cube size
; https://trac.sdss.org/wiki/MANGA/TRM/TRM_MPL-5/metadata#CorrelationMatrices
shape = repstr(repstr(sxpar(hcovs,'covshape'),'(',''),')','')
dimCOV = long(strsplit(shape,',',/ex))
covF = dblarr(dimCOV[0],dimCOV[1])
idxi = covs.indxi_c1+covs.indxi_c2*dim[0] ; 2d to 1d index
idxj = covs.indxj_c1+covs.indxj_c2*dim[0]
covF[idxi,idxj] = covs.rhoij
covF[idxj,idxi] = covs.rhoij
; transform the COV through binning (2x2 by default)
bin = 2 ; for 2x2 binning
; Binned Data (S) = T ## Input Data (F)
T = dblarr(dim[0]*dim[1],dim[0]*dim[1]/bin^2) 
for xbin=0,dim[0]/bin-1 do begin
	for ybin=0,dim[1]/bin-1 do begin
		tmp = intarr(dim[0],dim[1])
		tmp[xbin*bin:(xbin+1)*bin-1,ybin*bin:(ybin+1)*bin-1] = 1
		ind = where(tmp eq 1)
		T[ind, xbin+ybin*dim[0]/bin] = 1d/bin^2
	endfor
endfor
; compute the transformed COV matrix
; COV(S) = COV(T ## F) = T ## COV(F) ## transpose(T)
cov = T ## covF ## transpose(T)

; get VAR, which is easier to coadd than IVAR
var = ivar*0 + !values.f_nan 
ind = where(ivar ne 0)
var[ind] = 1./ivar[ind]

; object coordinates
ra = drpall.objra
dec = drpall.objdec

;;;;;;;;;;;;;;;
; RETIRED in Oct 2016, because masking f/g stars is DONE IN MPL-5
;;;;;;;;;;;;;;;
;; mask out other objects in the field
;; query SDSS catalog within 1'x1' box
;scatfile = getenv('SPFIT_DIR')+'/cats/sdss/'+basename+'.fits'
;if ~file_test(scatfile) then begin ; check if already exist
;	if verb then print,'--> querying SDSS catalog to mask out f/g stars'
;	sdsscat = QueryVizier('SDSS-DR9',[ra,dec],[1.,1.])
;	mwrfits,sdsscat,scatfile,/create
;endif else begin
;	if verb then print,'--> reading SDSS catalog to mask out f/g stars'
;	sdsscat = mrdfits(scatfile,1,/silent)
;endelse
;if size(sdsscat,/type) eq 8 then begin
;	; select primary sample to avoid duplicated sources
;	; also choose only stars brighter than 20 mags in any band and at dis > 3"
;	dis = sphdist(sdsscat.raj2000,sdsscat.dej2000,ra,dec,/deg)*3600
;	sdssmag = [[sdsscat.umag],[sdsscat.gmag],[sdsscat.rmag],$
;		[sdsscat.imag],[sdsscat.zmag]]
;	idx = where(sdsscat.mode eq 1 and sdsscat.cl eq 6 and $
;		min(sdssmag,dim=2,/nan) lt 20 and dis gt 3)
;	if idx[0] ne -1 then begin
;		sdsscat = sdsscat[idx]
;		adxy,hdr,sdsscat.raj2000,sdsscat.dej2000,x,y
;		; check if within the IFU field of view
;		dim = (size(flux))[1:3]
;		ind = where(x gt 0 and x lt dim[0]-1 and y gt 0 and y lt dim[1]-1,ct)
;		for i=0,ct-1 do begin
;			dist_ellipse,im,dim[0:1],x[ind[i]],y[ind[i]],$
;				1d,0.0,/double
;			im3d = congrid([[[im]],[[im]]],dim[0],dim[1],dim[2])
;			s = where(im3d le 4.0/0.5) ; w/i 4" radius
;			mask[s] = sdss_flagval('manga_drp3pixmask','forestar')
;		endfor
;		; delete these variables
;		;undefine,im3d,im
;	endif
;endif

; mask out objects based on object catalog that contain RA & Dec of
; objects to be masked out
if keyword_set(maskdir) then begin
	if file_test(maskdir+basename+'.fits') then begin
		maskcat = mrdfits(maskdir+basename+'.fits',1,/silent)
		adxy,hdr,maskcat.ra,maskcat.dec,x,y
		; check if within the IFU field of view
		ind = where(x gt 0 and x lt dim[0]-1 and y gt 0 and y lt dim[1]-1,ct)
		if ct gt 0 then begin
			for i=0,ct-1 do begin
				dist_ellipse,im,dim[0:1],x[ind[i]],y[ind[i]],$
					1d,0.0,/double
				im3d = congrid([[[im]],[[im]]],dim[0],dim[1],dim[2])
				; mask out pixels within 4 pixel = 2" radius
				s = where(im3d le 4.0,ct2)
				if ct2 gt 0 then ivar[s] = 0
			endfor
		endif
	endif
endif

; replace values with NaNs for pixels flagged as bad
; see idl/idlutils/data/sdss/sdssMaskbits.par for bitmask definitions
; this step is to avoid bad-pixel contamination in the 2x2 binning stage 
ind = where(ivar eq 0 $ 
	or ((mask and sdss_flagval('manga_drp3pixmask','nocov')) ne 0) $
	;or ((mask and sdss_flagval('manga_drp3pixmask','lowcov')) ne 0) $
	;or ((mask and sdss_flagval('manga_drp3pixmask','forestar')) ne 0) $
	;or ((mask and sdss_flagval('manga_drp3pixmask','donotuse')) ne 0) $
	or ((mask and sdss_flagval('manga_drp3pixmask','deadfiber')) ne 0),ct)
; save a single channel bitmask
bitmask0 = mask[*,*,dim[2]/2]
; simplify mask - either 0 (good pixel) or NaN (bad pixel)
; this step is necessary for binning the mask cubes to 1" pixels
mask *= 0.0 ; long -> float
; save original cube for later use
flux0 = flux 
var0 = var
; replace masked pixels with NaN so that they're considered missing data
if ct gt 0 then begin
	mask[ind] = !values.f_nan
	flux[ind] = !values.f_nan
	var[ind] = !values.f_nan
endif

if verb then print,'--> 2x2 rebinning '
; Rebin datacube from 0.5" spaxels to 1.0" spaxels to reduce covariance
; REBIN uses neighborhood averaging when minifying, but does not handle NAN properly
; BOXAVE3D: neighborhood averaging, treating NANs as missing data 
flux = boxave3d(flux,2,2,1)*4.0 ; x4 to preserve total flux of the datacube
var = boxave3d(var,2,2,1)*16.0	; x16 so that S/N remains the same after rebinning
				; neighboring spaxels so correlated that
				; S/N do not increase for 2x2 binning 
mask = boxave3d(mask,2,2,1)
; bin the original cubes
flux0 = boxave3d(flux0,2,2,1)*4.0
var0 = boxave3d(var0,2,2,1)*16.0
; replace !NaN back to 1 (= bad pixel)
mask[where(mask ne 0)] = 1
; mask out bad telluric lines, observed-frame, vacuum wavelength
hw = 10.0 ; wavelength range to mask (A)
skylines = [5579.,5895] ; 6300,6363
for i=0,n_elements(skylines)-1 do begin
	ind = where(wave gt skylines[i]-hw and wave lt skylines[i]+hw)
	mask[*,*,ind] = 1
endfor
; mask out first 40 pixels
dim = (size(flux))[1:3] ; reevaluate dim array after binning
mask[*,*,0:40] = 1

; update header keywords
sxaddpar,hdr,'cd1_1',sxpar(hdr,'cd1_1')*2
sxaddpar,hdr,'cd2_2',sxpar(hdr,'cd2_2')*2
; In WCS, the center of the lower-left pixel is (1.0,1.0)
; after 2x2 binning, the center of the lower-left pixel (1.0,1.0) 
; corresponds to the upper-right corner of the original first pixel (1.5,1.5)
; so we need to add 0.5 before dividing by 2
sxaddpar,hdr,'crpix1',(sxpar(hdr,'crpix1')+0.5)/2.
sxaddpar,hdr,'crpix2',(sxpar(hdr,'crpix2')+0.5)/2.

; compute tot(S/N) in Ha-NII region, over ~56 spectral channels
; (velscale=69 km/s for lstep_gal = 1e-4 in log10)
; dV = z*c ~= ln(lambda1/lambda0)*c = log10(lambda1/lambda0)*ln(10)*c
w1 = 6525 & w2 = 6610 ; r.f wave range = +/-2000 km/s or +/-28 pixels 
sigraw = fltarr(dim[0],dim[1]) ; array to keep S data
noiraw = fltarr(dim[0],dim[1]) ; array to keep N data
snrraw = fltarr(dim[0],dim[1]) ; array to keep S/N data
badfrac = total(mask,3)/dim[2] ; fraction of masked pixels
waverf = wave/(1.+objz)	 ; r.f. wavelength array
for i=0,dim[0]-1 do begin ; loop through all pixels
	for j=0,dim[1]-1 do begin
		ind = where(waverf gt w1 and waverf lt w2 and $
			mask[i,j,*] eq 0 and finite(flux[i,j,*]) and $
			finite(var[i,j,*]), ct)
		if ct gt 0 and badfrac[i,j] lt 0.2 then begin
			sigraw[i,j] = total(flux[i,j,ind])
			noiraw[i,j] = sqrt(total(var[i,j,ind]))
			snrraw[i,j] = sigraw[i,j]/noiraw[i,j]
		endif
	endfor
endfor
if max(snrraw) lt 5 then message,'Error: max(S/N) < 5'
; utilize covariance matrix for 2x2 binned cube (COV)
xind = range(0,dim[0]-1) # (fltarr(dim[1])+1) ; index of X-axis
yind = (fltarr(dim[0])+1) # range(0,dim[1]-1) ; index of Y-axis
if ~keyword_set(nobinning) and ~keyword_set(ds9reg) and $
	min(snrraw[where(snrraw gt 3)]) lt targetSN*1.5 then begin
	; do Voronoi binning only if min(S/N) < targetSN*1.5, otherwise no
	; binning is needed
	if verb then print,'--> Voronoi rebinning, targetSN = 1.5x ',targetSN
	; convert 2D -> 1D arrays
	; evaluate igood, sn, NODEs, binnum arrays, which are required in
	; subsequent analysis
	igood = where(snrraw gt 4.0) ; minimum SNR to consider
	signal = sigraw[igood]
	noise = noiraw[igood]
	x = xind[igood]
	y = yind[igood]
	; take the corresponding sub-COV-matrix 
	covsub = dblarr(n_elements(igood),n_elements(igood))
	for i=0,n_elements(igood)-1 do covsub[i,*] = cov[igood[i],igood]
	manga_voronoi_2d_binning, x, y, signal, noise, targetSN*1.5, binnum,$ ; size X
		xNode, yNode, xBar, yBar, sn, area, scale, $ ; size nBin
		COV=covsub,/QUIET
	if verb then print,'--> Number of Voronoi bins = ',n_elements(sn)
endif else if keyword_set(ds9reg) then begin
	; Use DS9 region to generate a binmap
	if keyword_set(regnum) then $ 
		regfile = ds9reg+basename+'-'+strc(regnum)+'.reg' $
	   else regfile = ds9reg+basename+'.reg'
	if verb then print,'--> Bin with DS9 region file: ',regfile
	; set a low targetSN so that all regions get fit
	if verb then print,'--> minimum SN = ',targetSN	
	; test if file exist
	if ~file_test(regfile) then message,'Region file cannot be found!'
	; 1st polygon filled w/ 1, 2nd polygon w/ 2, outside pixels w/ 0
	segmap = intarr(dim[0],dim[1])
	mask_polygon,segmap,regfile,value=0,nreg=nreg
	; evaluate igood, sn, NODEs, binnum arrays, which are required in
	; subsequent analysis
	; igood: index of spaxels in segmentation map
	igood = where(segmap ge 1,ct) 
	if ct eq 0 then message,'No poly regions defined' 
	; 1D spaxel arrays, dim = # of valid spaxels
	segmap2 = segmap[igood] 
	snrraw2 = snrraw[igood]	
	x = xind[igood]
	y = yind[igood]
	binnum = intarr(ct)-1
	; 1D bin arrays, dim = # of DS9 regions
	sn = fltarr(nreg)
	xNode = fltarr(nreg)
	yNode = fltarr(nreg)
	; loop through all DS9 regions
	for i=1,max(segmap2) do begin
		idx = where(segmap2 eq i,ct)
		; proceed only if there are spaxels in a region
		if ct gt 0 then begin
			; assign bin number
			binnum[idx] = i-1
			; unlike voronoi binning above, here we do a simple SNR 
			; calculation w/o correction for covariance
			sn[i-1] = sqrt(total(snrraw2[idx]^2))
			xNode[i-1] = mean(x[idx])
			yNode[i-1] = mean(y[idx])
		endif
	endfor
endif else begin
	; No binning
	if verb then print,'--> No binning, minimum SN = ',targetSN
	; evaluate igood, sn, NODEs, binnum arrays, which are required in
	; subsequent analysis
	igood = where(snrraw ge targetSN,ct)
	; if max(SNR) < targetSN, only fit the spaxel w/ the highest S/N
	if ct eq 0 then begin
		junk = max(snrraw,igood)
		ct = 1
	endif
	sn = snrraw[igood]
	xNode = xind[igood]
	yNode = yind[igood]
	binnum = indgen(ct)
endelse

if keyword_set(ds9reg) then begin
	; do not reorder if using DS9 region, so that the output file
	; matches that of the input DS9 region file
	sn_sorted = sn
	binnum_sorted = binnum
	nbin = n_elements(sn)
endif else begin
	; reorder bin # from closest to farthest from object center
	; note that the brightest pixel may not be the target galaxy
	adxy,hdr,ra,dec,x_ctr,y_ctr
	dis = (xNode-x_ctr)^2+(yNode-y_ctr)^2
	sidx = sort(dis)
	; apply the index
	sn_sorted = sn[sidx]
	binnum_sorted = binnum*0 ; size = # of valid spaxels
	nbin = n_elements(sn)
	for i=0,nbin-1 do begin
		idx = where(binnum eq sidx[i])
		binnum_sorted[idx] = i
	endfor
endelse

; skip this plateifu if maximum SNR is lower than target
if max(sn_sorted) lt targetSN then goto,theend

; save 2D binning map
binmap = intarr(dim[0],dim[1])-1
binmap[igood] = binnum_sorted
; save 2D binned S/N map
snrbin = fltarr(dim[0],dim[1])-1
for i=0,nbin-1 do snrbin[where(binmap eq i)] = sn_sorted[i]

; replace !NaN with original flux/var values to fill in the data of
; masked area. Use the mask array instead of NaNs to indicate bad pixels
ind = where(~finite(flux),ct)
if ct gt 0 then flux[ind] = flux0[ind]
ind = where(~finite(var),ct)
if ct gt 0 then var[ind] = var0[ind]

; load templates
; - match the spectral resolution to the data in rest-frame
; - resample the SSPs in log10 w/ the same dlogLam as the data
if verb then print,'--> loading & spectral matching template library ...'+ssplib
;mtplfile = '$MANGA_DIR/hfdap/ssp_matched/'+ssplib+'/'+basename+'.fits' ; saves the matched template library
mtplfile = mtpldir+basename+'-'+ssplib+'.fits' ; saves the matched template library
if ~file_test(mtplfile) then begin
	if verb then print,'--> Matching SSPs ...',mtplfile
	if ~file_test(mtpldir) then spawn,'mkdir '+mtpldir
	setup_templates,wave/(1.+objz),ssplib,/matchR,sres_data=sres,ssp_matched=ssp_matched
	mwrfits,ssp_matched,mtplfile,/create
endif else begin
	if verb then print,'--> Found matched SSPs ...',mtplfile
	ssp_matched = mrdfits(mtplfile,1,/silent)
endelse
if verb then print,'--> Number of templates ...',n_elements(ssp_matched)
; extract useful arrays from the structure
flam_tpl = ssp_matched.flam ; Unit: erg/s/Mo/AA
wave_tpl = ssp_matched[0].wave ; Unit: AA
; make sure SSP templates wavelength are in vacuum wavelength to avoid
; systematic velocity offset between gas and stars (Vgas - Vstar = -85 km/s)

; Tailor the wavelength ranges of templates and data to give
; maximum overlapping and minimum unnecessary data 
minmax_tpl = minmax(wave_tpl)
minmax_gal = minmax(wave/(1.+objz))
if verb then begin
	print,'Rest-Frame Wavelength range (templates vs. data), original'
	print,minmax_tpl,minmax_gal
end
; compute the overlapping range
tpl_range = [0.0,0.0]  
gal_range = [0.0,0.0]
buffer = 500.0 ; km/s, buffer to allow template to shift
c_kms = 299792.4580d ; Speed of light in km/s
if minmax_tpl[0] lt minmax_gal[0] then begin
	tpl_range[0] = minmax_gal[0]*(1.-buffer/c_kms)
	gal_range[0] = minmax_gal[0]
endif else begin
	tpl_range[0] = minmax_tpl[0]
	gal_range[0] = minmax_tpl[0]*(1.+buffer/c_kms) 
endelse
if minmax_tpl[1] gt minmax_gal[1] then begin
	tpl_range[1] = minmax_gal[1]*(1.+buffer/c_kms)
	gal_range[1] = minmax_gal[1]
endif else begin
	tpl_range[1] = minmax_tpl[1]
	gal_range[1] = minmax_tpl[1]*(1.-buffer/c_kms) 
endelse
; trim the data arrays to the overlapping wavelength range
w = where(wave/(1.+objz) ge gal_range[0] and wave/(1.+objz) le gal_range[1])
wave = wave[w]
sres = sres[w]
flux = flux[*,*,w]
var  = var[*,*,w]
mask = mask[*,*,w]
; trim template arrays to the same wavelength range
w = where(wave_tpl ge tpl_range[0] and wave_tpl le tpl_range[1])
wave_tpl = wave_tpl[w]
flam_tpl = flam_tpl[w,*]
; report the wavelength ranges after trimming
if verb then begin
	print,'Rest-Frame Wavelength range (templates vs. data), after matching'
	print,minmax(wave_tpl), minmax(wave/(1+objz))
end

; update header keywords after rebinning
dim = (size(flux))[1:3] ; update dimensions
sxaddpar,hdr,'naxis1',dim[0]
sxaddpar,hdr,'naxis2',dim[1]
sxaddpar,hdr,'crval3',alog10(wave[0])
sxaddpar,hdr,'naxis3',dim[2]

; progress report format
fmt = '(" '+basename+': ",i4,"/",i4," or ",f5.1,"% done",$,%"\r")'
; log10 wavelength arrays
log10lam = alog10(wave)
log10lam_tpl = alog10(wave_tpl)
; start looping through all voronoi bins
for ibin = 0, nbin-1, step do begin
	; find which spaxels belong to this bin
	ind = where(binmap eq ibin, n_spec)
	; skip bin if no spec in this bin or SNR < targetSN 
	if n_spec eq 0 or sn_sorted[ibin] lt targetSN then continue
	;note that sn_sorted[ibin] equals snrbin[where(binmap eq ibin)] 
	
	; find (x,y) indices
	ind2d = array_indices(binmap,ind)
	pix_x = reform(ind2d[0,*])
	pix_y = reform(ind2d[1,*])

	; combine log10 rebinned data
	galaxy2d = fltarr(n_spec,dim[2])
	var2d = fltarr(n_spec,dim[2])
	bitmask2d = fltarr(n_spec,dim[2])
	for i=0,n_spec-1 do begin
		galaxy2d[i,*] = flux[pix_x[i],pix_y[i],*]
		var2d[i,*] = var[pix_x[i],pix_y[i],*]
		bitmask2d[i,*] = mask[pix_x[i],pix_y[i],*]
	endfor
	; use MEAN to preserve surface brightness
	galaxy = mean(galaxy2d,dim=1,/nan)
	if ~keyword_set(nobinning) and n_spec gt 1 then begin
		diaCOV = total(COV[ind,ind]) ; diagnal elements only, i.e., 
				; the total err^2 when there is no covariance
		totCOV = 0.0 ; sum of all covariance, i.e., the true total err^2
		for i=0,n_elements(ind)-1 do totCOV += total(COV[ind[i],ind])
		if diaCOV eq 0 or totCOV eq 0 then corr=1.0 $
			else corr = float(sqrt(diaCOV/totCOV))
	endif else corr = 1.0
	; S/N = total(S)/sqrt(total(COV)) = <S>/sqrt<N^2> * sqrt(n_spec) * sqrt(diaCOV/totCOV)
	; so if we set S = <S> (= mean(S)), then
	; N = sqrt<N^2> / sqrt(n_spec) / sqrt(diaCOV/totCOV) 
	error = sqrt(mean(var2d,dim=1,/nan)) / sqrt(n_spec) / corr
	; bad pixel if more than half of the spaxels in a bin are bad
	bitmask = round(total(bitmask2d,1) / n_spec)
	; mask pixels w/ infinite errors
	ind = where(~finite(error),ct)
	if ct gt 0 then begin
		error[ind] = max(error,/nan)*1e3
		bitmask[ind] = 1 ; make sure to mask out these pixels
	endif

	; emission line setup structure
	if ibin eq 0 then begin
	  ; load the emission-line setup file, and creating the corresponding structure
	  readcol,linefile,eml_i,eml_name,eml_lambda,eml_action,eml_kind,eml_a,$
	  	eml_v,eml_s,eml_fit,f='(i,a,f,a,a,f,f,f,a)',comment='#'
	  ; this one includes broad lines
	  emissionbr = create_struct('i',eml_i,'name',eml_name,'lambda',eml_lambda,$
	  	'action',eml_action,'kind',eml_kind,'a',eml_a,'v',eml_v,'s',eml_s,'fit',eml_fit)
	  ; this one narrow lines only
	  ii = where(emissionbr.i lt 200)
	  emissionna = create_struct('i',eml_i[ii],'name',eml_name[ii],$
	  	'lambda',eml_lambda[ii],'action',eml_action[ii],'kind',eml_kind[ii],$
		'a',eml_a[ii],'v',eml_v[ii],'s',eml_s[ii],'fit',eml_fit[ii])
	  ; fit w/ broad lines and an additive poly (disabled 20160308) 
	  if keyword_set(BLR) then begin
	    xx = execute(cmd+',galaxy,error,log10lam,objz,flam_tpl,log10lam_tpl,'+$
	      'EMISSION_SETUP_IN=emissionbr,BITMASK=bitmask,SRES_DATA=sres,EBV_GAL=ebv_gal,'+$
	      'FIT_RESULTS=fit_results,EMISSION_SETUP_out=emissionout,degree=-1,_EXTRA=extra')
	    ; hold kinematics for broad lines for subsequent bins
	    i_br = where(strpos(emissionout.name,'_br') ge 0)
	    emissionbr.fit[i_br] = 'h'
	    emissionbr.v[i_br] = emissionout.v[i_br]
	    emissionbr.s[i_br] = emissionout.s[i_br]
	    ; hold additive poly for subsequent bins
	  endif
	endif

	if ~keyword_set(BLR) then $
	; fit w/ narrow lines only
	xx = execute(cmd+',galaxy,error,log10lam,objz,flam_tpl,log10lam_tpl,'+$
	  'EMISSION_SETUP_IN=emissionna,BITMASK=bitmask,SRES_DATA=sres,EBV_GAL=ebv_gal,'+$
	  'FIT_RESULTS=fit_results,EMISSION_SETUP_out=emissionout,_EXTRA=extra')
	if keyword_set(BLR) and ibin gt 0 then $
	; fit w/ broad lines, use line setup structure that has fixed broad lines
	xx = execute(cmd+',galaxy,error,log10lam,objz,flam_tpl,log10lam_tpl,'+$
	  'EMISSION_SETUP_IN=emissionbr,BITMASK=bitmask,SRES_DATA=sres,EBV_GAL=ebv_gal,'+$
	  'FIT_RESULTS=fit_results,EMISSION_SETUP_out=emissionout,_EXTRA=extra')
	
	; show fitting results
	if keyword_set(saveplot) then begin
		pars=fit_results
		plotfile = plotdir+'/'+string(ibin,f='(i04)')
		if keyword_set(ps) then begin
			; generate PS file to exam the quality of the fit
			mylegend = ssplib+' '+basename
			show_spec_fit,pars,mylegend=mylegend,/ps,$
				outfile=plotfile+'.ps'
		endif else begin
			mylegend = ssplib+' '+basename+' bin#:'+strc(ibin)
			if ibin eq 0 then begin
				device,decomp=0 
   				window,0,xs=1200,ys=800
   			endif
			show_spec_fit,pars,mylegend=mylegend
			save_screen,plotfile+'.png'
		endelse
	endif

	; create gpars array to store best-fit results
	if n_elements(gpars) eq 0 then begin
		tmp = fit_results
		zero_struct,tmp
		gpars = replicate(tmp,nbin)
	endif
	; save results for each bin
	gpars[ibin] = fit_results ; best-fit emission-lines + SSP models 

	; show progress
	if verb then print,f=fmt,ibin+1,nbin,100.*(ibin+1)/nbin
endfor

; split gpars into three chuncks
tags = tag_names(gpars)
i1 = where(tags eq 'LAMBDA')
i2 = where(tags eq 'CHI2NU')
pars_str = struct_selecttags(gpars,select_tags=tags[i1+1:i2])
spec_str = struct_selecttags(gpars,select_tags=tags[i2+1:*])

; compute RA, Dec, X, and Y positions of each bin
ra = dblarr(nbin)-99
dec = ra 
xbin = ra
ybin = ra
for kk=0,nbin-1 do begin
	ind = where(binmap eq kk,ct)
	if ct eq 0 then continue
	xy = array_indices(binmap,ind)
	xbin[kk] = mean(xy[0,*])
	ybin[kk] = mean(xy[1,*])
	xyad,hdr,xbin[kk],ybin[kk],bin_ra,bin_dec
	ra[kk] = bin_ra
	dec[kk] = bin_dec	
endfor
; emission line info
linename = strtrim(gpars[0].name,2)+strc(round(gpars[0].lambda))
linewave = gpars[0].lambda
; struct to save binmaps, SNR maps, etc.
vbin = {RA:ra,DEC:dec,X:xbin,Y:ybin,SNR:sn_sorted,$ ; 1D arrays, dim: # of bins
	snrraw:snrraw, bitmask:bitmask0,$ ; 2D arrays
	badfrac:badfrac, binmap:binmap, snrbin:snrbin, $
	linename:linename,linewave:linewave, $ ; fitted emission lines
	hdr:hdr} ; FITS header w/ updated WCS

; add tags to pars_str
toadd = mrd_struct(['RA','Dec','SNR','redshift'],$
	['-99d','-99d','-99.','-99.'],n_elements(gpars))
toadd.RA = ra
toadd.Dec = dec
toadd.SNR = sn_sorted
toadd.redshift = gpars.redshift
pars_str = struct_combine(toadd,pars_str)

; save best-fit pars and metadata
mwrfits,pars_str,outfile,/create ; [1] pars - emline and SSP pars
mwrfits,vbin,outfile             ; [2] vbin - linename, binmap, SNRmap, and WCS header
mwrfits,drpall,outfile 		 ; [3] drp  - drpall metadata

; save best-fit model spectra in a separate file
mwrfits,spec_str,repstr(outfile,'.fits','_spec.fits'),/create  ; [1] spec - observed and best-fit spectra

; print ellapsed time
theend:
if keyword_set(delete) then spawn,'rm -f '+infile
if verb then begin 
	print, '' ;; don't overwrite the final line
	toc
endif

end
