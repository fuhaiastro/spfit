pro show_spec_fit,pars,mylegend=mylegend,ps=ps,outfile=outfile

if ~keyword_set(mylegend) then mylegend=''
if keyword_set(ps) and ~keyword_set(outfile) then outfile='tmp.eps'

; load color map
loadct,0
; load object and model spectra
l_rf = 10.^pars.log10lam ; r.f. wavelength in A
obj = pars.galaxy ; flux array
err = pars.err    ; error array
fit = pars.best   ; best-fit model
emi = pars.emis   ; narrow+broad lines model
embr = pars.embr  ; broad-lines-only model
cfit = fit-emi    ; continuum-only model
res = obj-fit     ; residual 
good = pars.good ; goodpixels mask
chi2nu = pars.chi2nu ; reduced chi^2

; get rest-frame wavelength array
c = 299792.458d
;; deredshift, dV = cz ~= c*ln(1+z) => exp(dV/c) = 1+z
;l_rf = l/exp(pars.vel_star/c) 

; lambda plotting range
mn_l = min(l_rf,max=mx_l)
lrange = mx_l - mn_l 
lrange_rf_0 = [mn_l-0.02*lrange,mx_l+0.02*lrange]
; lambda plotting range zoom 1
lrange_rf_1 = [3700.1,4250.0]
region_1=textoidl('[OII]+[NeIII]+H\delta')
; lambda plotting range zoom 2
lrange_rf_2 = [4750.1,5100] ;5300.0]
region_2=textoidl('H\beta+[OIII]')
;lrange_rf_2 = [4750.1,5300.0]+870
;region_2='NaD region'
; lambda plotting range zoom 3
lrange_rf_3 = [6450,6800] ; [6250.0,6799.9]
region_3=textoidl('H\alpha+[NII]')
; flux plotting range
ind = where(good)
mx_f = max(obj[ind])
mn_f = -0.1*mx_f
frange = [mn_f,mx_f]
; to plot good pixels
good[where(good eq 0)] = -1
good = good*1e6
; plot positions
pos0=[0.070,0.59,0.970,0.98] ; main
posr=[0.070,0.49,0.970,0.59] ; residual
pos1=[0.070,0.08,0.337,0.42] ; zoom 1
pos2=[0.387,0.08,0.653,0.42] ; zoom 2
pos3=[0.703,0.08,0.970,0.42] ; zoom 3

; open PS file if necessary
if keyword_set(ps) then setps,outfile,29,21 ;,font='helvetica'

; set background color for legends
if keyword_set(ps) then begin
	bkgcolor=cgcolor('white')
	foncolor=cgcolor('black')
endif else begin
	bkgcolor=cgcolor('black')
	foncolor=cgcolor('white')
endelse

; plot entire wavelength range
xtitle=textoidl('\lambda_{r.f.} (A)')
ytitle=textoidl('f_\lambda (10^{-17} erg/s/cm^2/A)')
 plot,l_rf,obj,psym=10,xr=lrange_rf_0,/xs,yr=frange,/ys,xtitle='',ytitle=ytitle,pos=pos0,$
 	xtickname=replicate(' ',10)
oplot,l_rf,fit,color=cgcolor('red')
oplot,l_rf,err,color=cgcolor('sky blue')
oplot,l_rf,good,lines=1,color=cgcolor('green')
; useful legends
str = ['z = '+string(pars.redshift,f='(f6.4)'),$
	textoidl('\chi^2/DOF = ')+string(chi2nu,f='(f7.2)'),$
	'log(M*) = '+string(alog10(total(pars.m_star[*,0])),f='(f6.2)'),$
	textoidl('\sigma^* = ')+string(pars.kinstar[1,0],f='(f5.1)'),$
	'E(B-V)='+string(strmid(strcompress(round((pars.EBmV)[0]*1000.)/1000.,/remove_all),8,6,/reverse))]
al_legend,str,/top,/right,/normal,outline_color=foncolor,textcolor=foncolor,background=bkgcolor
al_legend,[mylegend],/top,/left,/normal,outline_color=foncolor,textcolor=foncolor,background=bkgcolor

; show residuals in subpanel
plot,l_rf,res/err,psym=3,xr=lrange_rf_0,yminor=-1,xticklen=0.08,$
	xtitle='',ytit=textoidl('\Delta/\sigma'),yr=[-10,10],/xs,/ys,pos=posr,/noerase
; show bad pixels
ind = where(good le 1e-6,ct)
if ct gt 0 then $
	oplot,l_rf[ind],res[ind]/err[ind],psym=4,syms=0.2,color=cgcolor('red')

; zoom region 1
ind = where(good and l_rf gt lrange_rf_1[0] and l_rf lt lrange_rf_1[1])
mx_f = max(obj[ind])*1.1
mn_f = -0.1*mx_f
frange = [mn_f,mx_f]
 plot,l_rf,obj,psym=10,xr=lrange_rf_1,/xs,yr=frange,/ys,xtickinterval=100,$
 	tit=region_1,xtitle=xtitle,ytitle=ytitle,pos=pos1,/noerase
oplot,l_rf,cfit,col=cgcolor('green')
oplot,l_rf,fit,col=cgcolor('red')
oplot,l_rf,emi,col=cgcolor('sky blue')
oplot,l_rf,res,psym=3
ind = where(good le 1e-6,ct)
if ct gt 0 then oplot,l_rf[ind],res[ind],psym=3,color=cgcolor('red')
oplot,l_rf,good,lines=1,color=cgcolor('green')

ytitle=''
; zoom region 2
ind = where(good and l_rf gt lrange_rf_2[0] and l_rf lt lrange_rf_2[1])
mx_f = max(obj[ind])*1.1
mn_f = -0.1*mx_f
frange = [mn_f,mx_f]
 plot,l_rf,obj,psym=10,xr=lrange_rf_2,/xs,yr=frange,/ys,xtickinterval=100,$
 	tit=region_2,xtitle=xtitle,ytitle=ytitle,pos=pos2,/noerase
oplot,l_rf,cfit,col=cgcolor('green')
oplot,l_rf,fit,col=cgcolor('red')
oplot,l_rf,emi,col=cgcolor('sky blue')
oplot,l_rf,embr,col=cgcolor('red'),lines=2
oplot,l_rf,res,psym=3
ind = where(good le 1e-6,ct)
if ct gt 0 then oplot,l_rf[ind],res[ind],psym=3,color=cgcolor('red')
oplot,l_rf,good,lines=1,color=cgcolor('green')

; zoom region 3
ind = where(good and l_rf gt lrange_rf_3[0] and l_rf lt lrange_rf_3[1])
mx_f = max(obj[ind])*1.1
mn_f = -0.1*mx_f
frange = [mn_f,mx_f]
 plot,l_rf,obj,psym=10,xr=lrange_rf_3,/xs,yr=frange,/ys,xtickinterval=100,$
 	tit=region_3,xtitle=xtitle,ytitle=ytitle,pos=pos3,/noerase
oplot,l_rf,cfit,col=cgcolor('green')
oplot,l_rf,fit,col=cgcolor('red')
oplot,l_rf,emi,col=cgcolor('sky blue')
oplot,l_rf,embr,col=cgcolor('red'),lines=2
oplot,l_rf,res,psym=3
ind = where(good le 1e-6,ct)
if ct gt 0 then oplot,l_rf[ind],res[ind],psym=3,color=cgcolor('red')
oplot,l_rf,good,lines=1,color=cgcolor('green')

; close PS file if necessary
if keyword_set(ps) then device,/close

end
