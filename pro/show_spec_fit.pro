pro show_spec_fit,pars,mylegend=mylegend,ps=ps,$
	outfile=outfile,lowc=lowc,obsframe=obsframe
; /lowc - low contrast mode, setting y-range with continuum model
; 	default is to set y-range with continuum+emission-line model

; default parameters
if ~keyword_set(mylegend) then mylegend=''
if keyword_set(ps) and ~keyword_set(outfile) then outfile='show_spec_fit.eps'

; Angstrom symbol
Angsym  = '!3'+String(197B)+'!X'

; The following ensures write_png works properly
; decomposed=0 to cause color values to be interpreted as indices into a color lookup table. 
device,decompose=0
; Having IDL provide the backing store (RETAIN=2) ensures that the window
; contents will be read properly.
if ~keyword_set(ps) then window,0,xs=600*2,ys=400*2,retain=2 $
	else setps,outfile,29,21 ;,font='helvetica'

; load color map
loadct,0
; load object and model spectra
l_rf = 10.^pars.log10lam ; r.f. wavelength in A
obj = pars.galaxy ; flux array
err = pars.err    ; error array
fit = pars.best   ; best-fit model
emi = pars.emis   ; narrow+broad lines model
embr= pars.embr  ; broad-lines-only model
cfit= fit-emi    ; continuum-only model
res = obj-fit     ; residual 
good= pars.good  ; goodpixels mask

; choose the array that sets the flux plotting range
if keyword_set(lowc) then yarr = cfit*2 else yarr = fit*1.5

; wavelength plotting ranges
; full spec
mn_l = min(l_rf,max=mx_l)
lrange = mx_l - mn_l 
lrange0 = [mn_l-0.02*lrange,mx_l+0.02*lrange]
; zoom regions
titles = textoidl(['[OII]+[NeIII]+H\delta','H\beta+[OIII]','H\alpha+[NII]'])
xtitle=textoidl('\lambda_{rest}')+'('+Angsym+')'
lrange = [[3700,4250],[4750,5100],[6450,6800]]

; observed frame wavelength vs. rest-frame wavelength
if keyword_set(obsframe) then begin
	l_rf *= (1+pars.redshift)
	lrange *= (1+pars.redshift)
	lrange0 *= (1+pars.redshift)
	xtitle=textoidl('\lambda_{obs}')+'('+Angsym+')'
endif

; modify good pixel mask to generate vertical dotted lines
good[where(good eq 0)] = -1
good = good*1e6

; plot positions
; full sepc
pos0=[0.070,0.59,0.970,0.98] ; main
posr=[0.070,0.49,0.970,0.59] ; residual
; zoom regions
pos1=[0.070,0.08,0.337,0.42] ; zoom 1
pos2=[0.387,0.08,0.653,0.42] ; zoom 2
pos3=[0.703,0.08,0.970,0.42] ; zoom 3
pzoom = [[pos1],[pos2],[pos3]]

; set background color for legends
if keyword_set(ps) then begin
	bkgcolor=cgcolor('white')
	foncolor=cgcolor('black')
	charsize=1.0
endif else begin
	bkgcolor=cgcolor('black')
	foncolor=cgcolor('white')
	charsize=1.5
endelse

; plot entire wavelength range
ytitle=textoidl('f_\lambda (10^{-17} erg/s/cm^2/')+Angsym+')'
; y range
ind = where(good)
yr = [-0.1,1]*max(yarr[ind])
plot,l_rf,obj,psym=10,xr=lrange0,/xs,yr=yr,/ys,$
	xtitle='',ytitle=ytitle,pos=pos0,$
 	xtickname=replicate(' ',10),chars=charsize
oplot,l_rf,fit,color=cgcolor('red')
oplot,l_rf,err,color=cgcolor('sky blue')
oplot,l_rf,good,lines=1,color=cgcolor('green')

; useful legends from pars file
str = ['z_in = '+string(pars.redshift,f='(f0.4)'),$
	textoidl('\chi^2/DOF = ')+string(pars.chi2nu,f='(f0.2)'),$
	'log(M*) = '+string(alog10(total(pars.m_star[*,0])),f='(f0.1)'),$
	textoidl('\sigma^* = ')+string(pars.kinstar[1,0],f='(i0)'),$
	textoidl('V^* = ')+string(pars.kinstar[0,0],f='(i0)'),$ ; ]
	'E(B-V)='+string(pars.EBmV[0],f='(f0.2)')]
al_legend,str,/top,/right,/normal,chars=charsize,$
	outline_color=foncolor,textcolor=foncolor,background=bkgcolor
al_legend,[mylegend],/top,/left,/normal,chars=charsize,$
	outline_color=foncolor,textcolor=foncolor,background=bkgcolor

; show residuals in subpanel
plot,l_rf,res/err,psym=3,xr=lrange0,yminor=-1,xticklen=0.08,chars=charsize,$
	xtitle='',ytit=textoidl('\Delta/\sigma'),yr=[-10,10],/xs,/ys,pos=posr,/noerase
; show bad pixels
ind = where(good le 1e-6,ct)
if ct gt 0 then $
	oplot,l_rf[ind],res[ind]/err[ind],psym=4,syms=0.2,color=cgcolor('red')

; zoom regions
for i=0,2 do begin
   ind = where(good and l_rf gt lrange[0,i] and l_rf lt lrange[1,i])
   frange = [-0.1,1]*max(yarr[ind])
    plot,l_rf,obj,psym=10,xr=lrange[*,i],/xs,yr=frange,/ys,chars=charsize,$
    	xtickinterval=200,tit=titles[i],xtitle=xtitle,$
	ytitle=(i eq 0)?ytitle:'',pos=pzoom[*,i],/noerase
   oplot,l_rf,cfit,col=cgcolor('green')
   oplot,l_rf,fit,col=cgcolor('red')
   oplot,l_rf,emi,col=cgcolor('sky blue')
   oplot,l_rf,res,psym=3
   ind = where(good le 1e-6,ct)
   if ct gt 0 then oplot,l_rf[ind],res[ind],psym=3,color=cgcolor('red')
   oplot,l_rf,good,lines=1,color=cgcolor('green')
endfor

; close PS file if necessary
if keyword_set(ps) then device,/close else save_screen,outfile

end
