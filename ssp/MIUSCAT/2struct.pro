; http://www.iac.es/proyecto/miles/pages/ssp-models.php
; Model version: 10.0, April 2015 
; Uniform resolution (FWHM = 2.51 A) from 3465-9469 A with 0.9A dispersion
; Unit: The total mass of the SSP is 1Mo. 
; The SSP spectra are given in units of Lo/Mo/Å, Lo = 3.826e33 erg/s
;
; Naming convention
; Spectral Range, IMF type, IMF slope, [M/H], age (Gyr), isochrone
; ku = Kroupa Universal
; e.g., MIUSCAT_BaSTI_ku_1.30_fits/Iku1.30Zm0.25T00.0300_iTp0.00_baseFe

tpldir = 'MIUSCAT_BaSTI_ku_1.30_fits'
; 636 SSPs, 53 ages * 12 [M/H] * 1 [a/Fe]

;tpldir = 'MIUSCAT_Padova00_ku_1.30_fits'
; 350 SSPs, 50 ages * 7 [M/H]

files = file_search(tpldir+'/*.fits',count=ct)
for i=0,ct-1 do begin
	ssp = mrdfits(files[i],0,h,/silent)
	if i eq 0 then begin
		npix = sxpar(h,'naxis1')
		; units: Gyr, [M/H], A, erg/s/Mo/A, none 
		tags = ['age','MH','wave','flam','source']
		vals = ['0.0','0.0','fltarr('+strc(npix)+')','fltarr('+strc(npix)+')','""']
		struct = mrd_struct(tags,vals,ct)
	endif
	; SSP wavelengths are in Air, convert to Vacuum
	wave_air = sxpar(h,'crval1')+sxpar(h,'cdelt1')*findgen(npix)
	airtovac,wave_air,wave_vac
	struct[i].wave = wave_vac
	struct[i].flam = ssp*3.826e33 ; Lo/Mo/A -> erg/s/A/Mo
	MH2 = strmid(files[i],strpos(files[i],'Z')+1,5)
	struct[i].MH = repstr(repstr(MH2,'m','-'),'p','+')
	age2 = strmid(files[i],strpos(files[i],'Z')+7,7)
	struct[i].age = age2
	struct[i].source = repstr(files[i],tpldir+'/','')
endfor
mwrfits,struct,tpldir+'.fits',/create
; add header keywords
hdr = headfits(tpldir+'.fits',ext=1)
sxaddpar,hdr,'Units','Gyr, [M/H], A, erg/s/A/Mo'
sxaddpar,hdr,'FWHM',2.51, ' spectral resolution in Angstrom'
mwrfits,struct,tpldir+'.fits',hdr,/create

; limit the number of SSPs for miuscat-baSTI
; 26 ages from 0.04 to 13.5 Gyr
; 3 metallicities: 0.5, 1.0, 2.0 Zsun
; 78 SSPs in total
if tpldir eq 'MIUSCAT_BaSTI_ku_1.30_fits' then begin
	struct = mrdfits(tpldir+'.fits',1,hdr)
	ind = where(struct.age lt 13.7 and (struct.MH eq -0.25 or struct.MH eq 0.06 or struct.MH eq 0.4))
	struct2 = struct[ind[1:*:2]]
	mwrfits,struct2,tpldir+'-thin.fits',hdr,/create
endif 

end
