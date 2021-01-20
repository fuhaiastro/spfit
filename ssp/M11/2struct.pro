; From Maraston 2011 distribution
; STELIB-BASED MODELS: These models are based on the stellar library STELIB by 
;                      Le Borgne et al. (2003).
; * Cool high-resolution theoretical stars (Gustafsson et 
;   al. 2008) have been added to counter inadequacies in 
;   the STELIB library. 
; * Many spectra in the library are not complete over the quoted 
;   wavelength range, which limits the actual wavelength range of the youngest 
;   ages at solar, and all ages at non-solar metallicities.
; * Somewhat coarse sampling of stellar parameter space at 
;   non-solar metallicities may affect the reliability of these models.
; 
; Resolution: 3.1-3.4 AA (fwhm) MDAP adopts a value of 3.4 AA
; https://trac.sdss.org/wiki/MANGA/Data/DAPDevelopment/TemplateLibraries
;
; 0.5 AA sampling
;
; * 'zxxx': chemical composition, as it follows (notation identical to M05)
;     z10m4=0.0001
;     z0001=0.001
;     z001=0.01
;     z002=0.02 (solar metallicity)
;     z004=0.04
;
; * 'imf': Initial Mass Funcion,IMF, namely:
;    '.ss' ==> Salpeter IMF
;    '.kr' ==> Kroupa IMF
;    '.cha' ==> Chabrier IMF
;
; ssp_M11_STELIB.'imf'z001: 24 ages, 200 Myr to 15 Gyr
; 				         3201.0 - 9296.5 AA
; 				         12192 flux points
; 				         !Usable wavelength range for all ages is
; 				         3201.0 - 7900.0 AA!
; 
; ssp_M11_STELIB.'imf'z002: 39 ages, 30 Myr to 15 Gyr
; 				         3201.0 - 9296.5 AA
; 				         12192 flux points
; 
; ssp_M11_STELIB.'imf'z004: 22 ages, 400 Myr to 15 Gyr
; 				        3201.0 - 9296.5 AA
; 				        12192 flux points
; 				        !Usable wavelength range for all ages is
; 				        3201.0 - 7900.0 AA!
; 
; COLUMNS:
; ------------------------------------------------------------------------------
; Age (Gyr) | [Fe/H] | Wavelength (AA) | Flux (ergs /s /AA /Msun) |
; last column is the luminosity per AA for a 1 solar mass SSP
; ------------------------------------------------------------------------------

tpldir = 'SSP_M11_STELIB'
imf = 'kr'
m11file = file_search(tpldir+'/*.'+imf+'*',COUNT=nfiles)

; make struct array 
npix = 12192
; units: Gyr, [M/H], A, erg/s/A/Mo, none 
tags = ['age','MH','wave','flam','source']
vals = ['0.0','0.0','fltarr('+strc(npix)+')','fltarr('+strc(npix)+')','""']
struct = mrd_struct(tags,vals,24+39+22)

kk=0
for k=0,nfiles-1 do begin
	readcol,m11file[k],age,FeH,lam2,ssp_flux
	; wavelength is given in Air, convert to Vac
	airtovac,lam2,lam2vac
	; how many different ages?
	ind = rem_dup(age) 
	print,m11file[k],n_elements(ind)
	for i=0,n_elements(ind)-1 do begin
		s = where(age eq age[ind[i]])
	   	struct[kk].flam = ssp_flux[s]
		struct[kk].wave = lam2vac[s]
		struct[kk].age = age[ind[i]]
		struct[kk].MH = FeH[0]
		struct[kk].source = m11file[k]
		kk +=1
	endfor
endfor
tpldir = tpldir + '_' + imf
mwrfits,struct,tpldir+'.fits',/create
; add header keywords
hdr = headfits(tpldir+'.fits',ext=1)
sxaddpar,hdr,'Units','Gyr, [M/H], A, erg/s/A/Mo'
sxaddpar,hdr,'FWHM',3.4, ' spectral resolution in Angstrom'
mwrfits,struct,tpldir+'.fits',hdr,/create

end
