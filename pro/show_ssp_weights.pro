; given the fitting results from SPFIT and the SSP library 
; show the weight distribution in Age, [M/H] plane

pro show_ssp_weights,infile

basename = exfilename(infile,dirname=outdir,/noextension)
ssplib = 'miuscat' ; SSP templates
mtplfile = outdir+basename+'_'+ssplib+'.fits'
ssp_matched = mrdfits(mtplfile,1,/silent)
spfitfile = outdir+basename+'_spfit.fits'
fit_results = mrdfits(spfitfile,1,/silent)

; start plotting
outfile = repstr(mtplfile,'.fits','.eps')
setps,outfile,10,10,font='Helvetica'

tit = textoidl(' \sigma^* = ')+$
	string(fit_results.kinstar[1,0],f='(f5.1)')+$
	textoidl('\pm')+$
	string(fit_results.kinstar[1,1],f='(f5.1)')+' km/s'

x = ssp_matched.age
y = ssp_matched.MH
z = fit_results.m_star[*,0]/total(fit_results.m_star[*,0])
plot,x,y,psym=3,xtit='Age (Gyr)',ytit='[M/H]',tit=tit,/xlog,/xs,/ys,$
	xr=minmax(x)*[0.5,2],yr=minmax(y)+[-0.2,0.2],position=[0.2,0.15,0.98,0.90]
ind = where(z gt 0,ct)
plots,x[ind],y[ind],psym=6,syms=0.5,color=cgcolor('sky blue')
for i=0,ct-1 do $
plots,x[ind[i]],y[ind[i]],psym=8,syms=5*z[ind[i]]/max(z)+1,color=cgcolor('red')
;contour,z,x,y,/irregular,/overplot

device,/close

end

