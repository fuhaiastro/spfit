pro show_spec_fit,pars,mylegend=mylegend,ps=ps,$
	outfile=outfile,lowc=lowc
; /lowc - low contrast mode, setting yrange by using continuum model

if ~keyword_set(mylegend) then mylegend=''
if keyword_set(ps) and ~keyword_set(outfile) then outfile='show_spec_fit.eps'

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

; choose the array that sets the flux plotting range
if keyword_set(lowc) then yarr = cfit*2 else yarr = obj

; wavelength plotting ranges
; full spec
mn_l = min(l_rf,max=mx_l)
lrange = mx_l - mn_l 
lrange0 = [mn_l-0.02*lrange,mx_l+0.02*lrange]
; zoom regions
titles = textoidl(['[OII]+[NeIII]+H\delta','H\beta+[OIII]','H\alpha+[NII]'])
lrange = [[3700,4250],[4750,5100],[6450,6800]]

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
; y range
ind = where(good)
yr = [-0.1,1]*max(yarr[ind])
plot,l_rf,obj,psym=10,xr=lrange0,/xs,yr=yr,/ys,$
	xtitle='',ytitle=ytitle,pos=pos0,$
 	xtickname=replicate(' ',10)
oplot,l_rf,fit,color=cgcolor('red')
oplot,l_rf,err,color=cgcolor('sky blue')
oplot,l_rf,good,lines=1,color=cgcolor('green')

; useful legends from pars file
str = ['z = '+string(pars.redshift,f='(f6.4)'),$
	textoidl('\chi^2/DOF = ')+string(pars.chi2nu,f='(f7.2)'),$
	'log(M*) = '+string(alog10(total(pars.m_star[*,0])),f='(f6.2)'),$
	textoidl('\sigma^* = ')+string(pars.kinstar[1,0],f='(f5.1)'),$
	textoidl('V^* = ')+string(pars.kinstar[0,0],f='(f5.1)')]
	;'E(B-V)='+string(strmid(strcompress(round((pars.EBmV)[0]*1000.)/1000.,/remove_all),8,6,/reverse))]
al_legend,str,/top,/right,/normal,outline_color=foncolor,textcolor=foncolor,background=bkgcolor
al_legend,[mylegend],/top,/left,/normal,outline_color=foncolor,textcolor=foncolor,background=bkgcolor

; show residuals in subpanel
plot,l_rf,res/err,psym=3,xr=lrange0,yminor=-1,xticklen=0.08,$
	xtitle='',ytit=textoidl('\Delta/\sigma'),yr=[-10,10],/xs,/ys,pos=posr,/noerase
; show bad pixels
ind = where(good le 1e-6,ct)
if ct gt 0 then $
	oplot,l_rf[ind],res[ind]/err[ind],psym=4,syms=0.2,color=cgcolor('red')

; zoom regions
for i=0,2 do begin
   ind = where(good and l_rf gt lrange[0,i] and l_rf lt lrange[1,i])
   frange = [-0.1,1]*max(yarr[ind])
    plot,l_rf,obj,psym=10,xr=lrange[*,i],/xs,yr=frange,/ys,xtickinterval=100,$
    	tit=titles[i],xtitle=xtitle,ytitle=(i eq 0)?ytitle:'',pos=pzoom[*,i],/noerase
   oplot,l_rf,cfit,col=cgcolor('green')
   oplot,l_rf,fit,col=cgcolor('red')
   oplot,l_rf,emi,col=cgcolor('sky blue')
   oplot,l_rf,res,psym=3
   ind = where(good le 1e-6,ct)
   if ct gt 0 then oplot,l_rf[ind],res[ind],psym=3,color=cgcolor('red')
   oplot,l_rf,good,lines=1,color=cgcolor('green')
endfor

; close PS file if necessary
if keyword_set(ps) then device,/close

end
