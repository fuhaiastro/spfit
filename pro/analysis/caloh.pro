;+
; Purpose
; 	Measure gas-phase 12+log(O/H) with various calibrators
;	- calculate joint BPT+WHAN classifications
;	- keep log line ratios, and EW(Ha) 
;	- for SF+Comp objs, deredden line fluxes, and calculate
;	12+log(O/H) using various calibrators
; input
;	spftpars - structure array, SPFIT best-fit parameters
;	vbin - structure, stores linewave and linename
;
; output
;	outfile - string, filename of the output file
; history
;	8/16/2019 - written by HF
;	11/21/2019 - added Marino13 O3N2 calibrator
;-
function caloh, spftpars, vbin, outfile=outfile

; Ha/Hb-derived intrinsic reddening, 12+log(O/H)
tags = ['BPT','EBV',$
	'WHa','dWHa',$     	; de-redshifted H-alpha EW in Angstrom
	'O3Hb' ,'N2Ha' ,'O1Ha' ,'S2Ha' ,$  ; log(line ratios)
	'dO3Hb','dN2Ha','dO1Ha','dS2Ha',$  ; errors
	'R23','O3N2',$
	'OH_D16','OH_D16_err',$
	'OH_T04','OH_T04_err',$
	'OH_PP04','OH_PP04_err',$
	'OH_M13','OH_M13_err',$
	'OH_M08R23','OH_M08N2']
vals = ['-9',replicate('-9.9',n_elements(tags)-1)] 
nobj = n_elements(spftpars)
str = mrd_struct(tags,vals,nobj)

; get emission line indices 
linename = strtrim(vbin.linename,2)
linewave = vbin.linewave
io2 = where(linename eq 'OII3730')  ; 3727 tied to 3729 @ 1:1
io3 = where(linename eq 'OIII5008') ; 4960 tied to 5008 @ 0.35:1
ihb = where(linename eq 'Hb4863')
in2 = where(linename eq 'NII6585')  ; 6549 tied to 6585 @ 0.34:1
iha = where(linename eq 'Ha6565')
io1 = where(linename eq 'OI6302')
is2a= where(linename eq 'SII6718') 
is2b= where(linename eq 'SII6733')

; BPT+WHAN classification
; emission line fluxes
fn2 = reform(spftpars.flux[in2,0])
fha = reform(spftpars.flux[iha,0])
fo3 = reform(spftpars.flux[io3,0])
fhb = reform(spftpars.flux[ihb,0])
fo1 = reform(spftpars.flux[io1,0])
fs2a= reform(spftpars.flux[is2a,0])
fs2b= reform(spftpars.flux[is2b,0])
; fractional errors, df/f
dfn2 = emlerr(spftpars.aon[in2], spftpars.sigma_obs[in2])
dfha = emlerr(spftpars.aon[iha], spftpars.sigma_obs[iha])
dfo3 = emlerr(spftpars.aon[io3], spftpars.sigma_obs[io3])
dfhb = emlerr(spftpars.aon[ihb], spftpars.sigma_obs[ihb])
dfo1 = emlerr(spftpars.aon[io1], spftpars.sigma_obs[io1])
dfs2a= emlerr(spftpars.aon[is2a],spftpars.sigma_obs[is2a])
dfs2b= emlerr(spftpars.aon[is2b],spftpars.sigma_obs[is2b])
dfs2 = sqrt((dfs2a*fs2a)^2+(dfs2b*fs2b)^2)/(fs2a+fs2b)
; line ratios
str.n2ha = alog10(fn2/fha)
str.o3hb = alog10(fo3/fhb)
str.o1ha = alog10(fo1/fha)
str.s2ha = alog10((fs2a+fs2b)/fha)
; error propagation: d(log a) = (1/ln10) (da/a)
str.dn2ha = sqrt(dfn2^2+dfha^2)/alog(10)
str.do3hb = sqrt(dfo3^2+dfhb^2)/alog(10)
str.do1ha = sqrt(dfo1^2+dfha^2)/alog(10)
str.ds2ha = sqrt(dfs2^2+dfha^2)/alog(10)
; EW
str.wha  = reform(spftpars.EW[iha])
str.dwha = reform(spftpars.EW[iha]) * dfha
; NII/Ha BPT classification
; BPT	CODE - Seyfert=4, LINER=3, Comp=2, SF=1, RG=0, Invalid=-2
bptclass = bpt_k06_n2(str.o3hb,str.n2ha,ew=str.wha,code=bptcode)
str.bpt = bptcode

;;;;;;;;;;;;;;;;;;;;;
; 12+log(O/H) estimates for SF+Comp sources
;;;;;;;;;;;;;;;;;;;;;
sc = where(str.bpt eq 1 or str.bpt eq 2,ct)
if ct gt 0 then pars = spftpars[sc] else goto,theend 

; reddening correction 
hahb_obs = reform(pars.flux[iha,0]/pars.flux[ihb,0])
; intrinsic line ratio (Kewley06)
hahb_int = 2.85
; for Av = 1.0, calculate A_lambda at H-alpha and H-beta?
; ccm_UNRED,[6564.61,4862.68],[1.0,1.0],1.0/3.1,funred,R_V=3.1,A_Lambda=Al
; => Al[0] = A_6565/Av = 0.812214, Al[1] = A_4863/Av = 1.18113
; => Av = 2.5/(Al[1]-Al[0]) * log[(Ha/Hb)/2.85] = 6.7766447 log[(Ha/Hb)/2.85]
; => Ha,int/Ha,obs = [(Ha/Hb)/2.85]^[Al[0]/(Al[1]-Al[0])] = [Ha/Hb)/2.85]^2.2016335
; fast way to compute Av
av = alog10(hahb_obs/hahb_int)*6.7766447
ebv = av/3.1 ; Rv = 3.1
; unredden all lines
for i=0,n_elements(ebv)-1 do begin
	if ebv[i] le 0 then continue
	ccm_unred,linewave,pars[i].flux[*,0],ebv[i],funred
	pars[i].flux[*,0] = funred
endfor
str[sc].ebv = ebv

; take reddening-corrected emission line fluxes
; correct fluxes for doublets
fo2 = reform(pars.flux[io2,0])*2.0  ; 3727 tied to 3729 @ 1:1
fn2 = reform(pars.flux[in2,0])
fha = reform(pars.flux[iha,0])
fo3 = reform(pars.flux[io3,0])*1.35 ; 4960 tied to 5008 @ 0.35:1
fhb = reform(pars.flux[ihb,0])
fs2a= reform(pars.flux[is2a,0])
fs2b= reform(pars.flux[is2b,0])
; fractional errors, df/f
dfo2 = emlerr(pars.aon[io2],pars.sigma_obs[io2])
dfn2 = emlerr(pars.aon[in2],pars.sigma_obs[in2])
dfha = emlerr(pars.aon[iha],pars.sigma_obs[iha])
dfo3 = emlerr(pars.aon[io3],pars.sigma_obs[io3])
dfhb = emlerr(pars.aon[ihb],pars.sigma_obs[ihb])
dfs2a= emlerr(pars.aon[is2a],pars.sigma_obs[is2a])
dfs2b= emlerr(pars.aon[is2b],pars.sigma_obs[is2b])

; Dopita+16 calibrator using NII/Ha and NII/SII
; [S II]6718 + 6733
fs2  = fs2a+fs2b
dfs2 = sqrt((dfs2a*fs2a)^2+(dfs2b*fs2b)^2)/(fs2a+fs2b)
n2ha = alog10(fn2/fha)
n2s2 = alog10(fn2/fs2)
; error propagation: d(log a) = (1/ln10) (da/a)
dn2ha = sqrt(dfn2^2+dfha^2)/alog(10)
dn2s2 = sqrt(dfn2^2+dfs2^2)/alog(10)
; Dopita+16 calibrator
y = n2s2+0.264*n2ha
dy = sqrt(dn2s2^2+(0.264*dn2ha)^2)
str[sc].OH_D16 = 8.77+y
str[sc].OH_D16_err = dy

; Tremonti+04 calibrator using R23
; FO23 = [O II]3726,3729 + [O III]4959,5007
fo23  = fo2+fo3
dfo23 = sqrt((dfo2*fo2)^2+(dfo3*fo3)^2)/(fo2+fo3)
R23 = alog10(fo23/fhb)
dR23 = sqrt(dfo23^2+dfhb^2)/alog(10)
str[sc].R23 = R23
str[sc].OH_T04 = 9.185-0.313*R23-0.264*R23^2-0.321*R23^3
str[sc].OH_T04_err = 0.313*dR23

; O3N2 calibrators: O3N2 = log(OIII5007/Hb / NII6584/Ha)
; Pettini & Pagel 2004
O3N2 = alog10((fo3/1.35/fhb)/(fn2/fha))
str[sc].O3N2 = O3N2
dO3N2 = sqrt(dfo3^2+dfhb^2+dfn2^2+dfha^2)/alog(10)
str[sc].OH_PP04 = 8.73-0.32*O3N2
str[sc].OH_PP04_err = 0.32*dO3N2
; Marino+13 
str[sc].OH_M13 = 8.533-0.214*O3N2
str[sc].OH_M13_err = 0.214*dO3N2

; Maiolino+08 M08N2
nobj = n_elements(sc)
OH_M08N2 = fltarr(nobj)
coeffs = [-0.7732,1.2357,-0.2811,-0.7201,-0.3330]
for i=0,nobj-1 do begin
	coeffs[0] = -0.7732-n2ha[i]
	roots = fz_roots(coeffs)
	; take only real roots within the calibrated range
	ind = where(imaginary(roots) eq 0 and $
		real_part(roots)+8.69 gt 7.0 and $
		real_part(roots)+8.69 lt 9.5,ct)
	if ct gt 0 then	begin
		x = (real_part(roots[ind]))[0]
		OH_M08N2[i] = x+8.69
	endif ; else print,i,n2ha[i],' M08N2 no real solution'
endfor
str[sc].OH_M08N2 = OH_M08N2

; Maiolino+08 R23 calibrator
OH_M08r23 = fltarr(nobj)
coeffs = [0.7462,-0.7149,-0.9401,-0.6154,-0.2524]
for i=0,nobj-1 do begin
	coeffs[0] = 0.7462-R23[i]
	roots = fz_roots(coeffs)
	; take only real roots within the calibrated range
	ind = where(imaginary(roots) eq 0 and $
		real_part(roots)+8.69 gt 7.0 and $
		real_part(roots)+8.69 lt 9.5,ct)
	; use NII/Ha-derived O/H to choose between the upper and lower O/H branches
	; 8.08 corresponds to the peak of the polynomial
	if ct gt 0 then begin
		if OH_M08N2[i] lt 8.08 and OH_M08N2[i] gt 0 then $
			x = min(real_part(roots[ind])) else $
			x = max(real_part(roots[ind]))
		OH_M08r23[i] = x+8.69
	endif ; else print,i,R23[i],' M08R23 no real solution'
endfor
str[sc].OH_M08r23 = OH_M08r23

theend:
if keyword_set(outfile) then mwrfits,str,outfile,/create
return,str

end
