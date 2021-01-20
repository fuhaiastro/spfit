pro stack, xdata, ydata, xgrid, xctr=xctr, xlerr=xlerr, xherr=xherr, $
	ymean=ymean, yerr=yerr, npbin=npbin

; bin center X coordinates
xctr = xgrid[sort(xgrid)]

; lower and upper bounds for each bin
nbin = n_elements(xctr)
delx = xctr[1:*]-xctr
xhig = xgrid+[delx,delx[nbin-2]]/2
xlow = xgrid-[delx[0],delx]/2

; bin X range bars
xherr = xhig-xctr
xlerr = xctr-xlow

; compute mean and stddev per bin
ymean = fltarr(nbin)+!values.f_nan
yerr = ymean
npbin = intarr(nbin)
for i=0,nbin-1 do begin
	s = where(xdata ge xlow[i] and xdata lt xhig[i],ct)
	if ct gt 0 then begin
		npbin[i] = ct
		ymean[i]= median(ydata[s])
		if ct gt 10 then begin
			yerr[i] = robust_sigma(ydata[s],GOODVEC=good)
		endif else begin
			yerr[i] = stddev(ydata[s])
		endelse
	endif
endfor

return

end
