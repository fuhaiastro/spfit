;+
; Purpose
;	extract basic information from SPFIT parameter file
; input
;	spftpars - structure array, SPFIT best-fit parameters
;	vbin - structure giving linewave, linename, and binmap
;	drp  - structure giving the plateifu's DRPALL metadata
;	ssp  - structure giving SSP library info
;
; output
;	str  - structure
; optional output
;	outfile - string, filename of the output file
; history
;	8/30/2019 - written by HF
;-
function calbasic, pars, vbin, drp, ssp, outfile=outfile

; useful tags to get
tags = ['plateifu','RA','Dec','z',$ ; coordinates and abs-line redshift
	'R_in_R50','R_in_R50b',   $ ; deprojected radius in R_eff
	'logM','dlogM','Age','MH',$ ; stellar mass in log(Msun)
	'FO3','dFO3','logLo3',$	    ; [OIII]5007 line flux: 1d-17 erg/s/cm2, df/f, erg/s
	'HaHb','logLHa','logLO2',$  ; SFR estimators, not corrected for reddening
	'Vgas' ,'Sgas' ,'Vstar' ,'Sstar',$ ; gas & stellar kinematics
	'dVgas','dSgas','dVstar','dSstar'] ; errors
vals = ['""','-9.9d','-9.9d','-9.9d',replicate('-9.9',n_elements(tags)-4)] 
str = mrd_struct(tags,vals,n_elements(pars))

; line indices
linename = strtrim(vbin.linename,2)
io2 = where(linename eq 'OII3730')  ; 3727 tied to 3729 @ 1:1
io3 = where(linename eq 'OIII5008') ; 4960 tied to 5008 @ 0.35:1
ihb = where(linename eq 'Hb4863')
iha = where(linename eq 'Ha6565')

; save other useful parameters
str.plateifu = replicate(drp.plateifu,n_elements(pars))
str.RA = vbin.ra
str.dec = vbin.dec
; redshift - dz = dv/c
str.z = drp.nsa_z + pars.kinstar[0,0]/299792.46 ; c_kms
; kinematics
str.Vgas = reform(pars.vel[iha,0])
str.Sgas = reform(pars.sigma[iha,0])	
str.Vstar = pars.kinstar[0,0]
str.Sstar = pars.kinstar[1,0]
; errors
str.dVgas = pars.vel[0,1]
str.dVstar = pars.kinstar[0,1]
str.dSstar = pars.kinstar[1,1]
; sig_obs d sig_obs = sig_int d sig_int
; => d sig_int = (sig_obs^2/sig_int) * (d sig_obs / sig_obs)
velscale = 69.0 ; km/s per pixel
str.dSgas = reform((pars.sigma_obs[iha]*velscale)^2/pars.sigma[iha,0] * $
	     emlerr(pars.aon[iha],pars.sigma_obs[iha],/sigma))

; AGN luminosity estimator
; [O III] fluxes - 1d-17 erg/s/cm2
; not corrected for reddening
str.fo3 = reform(pars.flux[io3,0])
; fractional error - df/f
str.dfo3 = emlerr(pars.aon[io3],pars.sigma_obs[io3])
; 1d-17 erg/s/cm2 * Mpc^2 -> erg/s
; 1d-17 * 4 * !PI * (3.086d24)^2 = 1.1967453e+33
f2l = alog10((lumdist(drp.nsa_z,/silent))^2*1.197d33)
str.logLo3 = alog10(str.fo3)+f2l

; SFR estimators
str.HaHb = reform(pars.flux[iha,0]/pars.flux[ihb,0])
str.logLHa = alog10(reform(pars.flux[iha,0]))+f2l
str.logLO2 = alog10(reform(pars.flux[io2,0])*2)+f2l

; stellar mass (logMsun), age (Gyr), and [M/H]
str.logM = alog10(total(pars.m_star[*,0],1))
str.dlogM = sqrt(total(pars.m_star[*,1]^2))/total(pars.m_star[*,0],1)/alog(10)
; SSP parameters
str.age = total(ssp.age*pars.m_star[*,0])/total(pars.m_star[*,0])
str.MH = alog10(total(10.^ssp.mh*pars.m_star[*,0])/total(pars.m_star[*,0]))

; deprojected radius
adxy,vbin.hdr,drp.objra,drp.objdec,xc,yc ; object position
BA = drp.NSA_ELPETRO_BA
q = 0.13 ; oblateness at face-on (Giovanelli94)
if BA gt q then BA_int = sqrt((BA^2-q^2)/(1-q^2)) else BA_int = BA ; correct for oblateness
pos_ang = drp.NSA_ELPETRO_PHI    ; deg
r50     = drp.NSA_ELPETRO_TH50_R ; arcsec
scale   = 1.0 ; arcsec/pixel, scale of the spaxel
dist_ellipse,mask,(size(vbin.binmap))[1:2],xc,yc,1./BA_int,pos_ang
; no oblateness correction
dist_ellipse,mas2,(size(vbin.binmap))[1:2],xc,yc,1./BA,pos_ang
nbin = n_elements(pars)
for kk=0,nbin-1 do begin
	ind = where(vbin.binmap eq kk,ct)
	if ct gt 0 then begin
		str[kk].r_in_r50 = mean(mask(ind))*scale/r50
		str[kk].r_in_r50b= mean(mas2(ind))*scale/r50
	endif
endfor

if keyword_set(outfile) then mwrfits,str,outfile,/create
return,str

end

