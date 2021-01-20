PRO SPFIT_SET_CONSTRAINTS, GALAXY=galaxy, NOISE=noise,L0_GAL=l0_gal, LSTEP_GAL=lstep_gal,$
	templates=templates, EMISSION_SETUP=emission_setup, L0_TEMPL=l0_templ, $
        GOODPIXELS=goodpixels, INT_DISP=int_disp, velscale=velscale, $
	START_PARS=start_pars, PARINFO=parinfo, FUNCTARGS=functargs, $
	SMOMENTS=smoments, GMOMENTS=gmoments, BIAS=bias

; This procedure sets up the constraints and boundaries for the
; variables to be fitted, preparing and returning the PARINFO and
; FUNCTARG structure for MPFIT
;
; Total number of parameters that MPFIT will deal with equal to 2x/3x the
; number of main emission lines (counting multiplets only once)
; plus the order of the multiplicative polynomials and, if
; needed, the reddening parameters. PARINFO is set up so that it
; properly account for the hold and tied parameters in the EMISSION_SETUP
; structure

i_lines = where(emission_setup.kind eq 'l')
nlines  = n_elements(i_lines)
nssp = (size(templates))[2]
npars = n_elements(start_pars)

; Setup the PARINFO structure that will allow us to control the limits
; of our parameters with MPFIT, decide whether we want to hold them at
; the input values, or tie some of them to others.
; set step following pPXF
parinfo = REPLICATE({step:1d-3,limits:[0d,0d],limited:[1,1],fixed:0, tied:' '}, npars)

; A) First of all, fix kinematics to their input values for
; the lines with a 'h' (for hold) as their FIT tag
for i = 0,nlines-1 do begin
    j = i_lines[i]
    if (emission_setup.fit[j] eq 'h') then begin
        parinfo[5*i+1:5*i+4].fixed = 1
    endif
endfor

; B) set the limits, and steps
; i) for centroid (pix), sigma (pix), h3, h4
for i=0,nlines*5 -1, 5 do begin
    parinfo[i+1].limits = start_pars[i+1] + [-1d3,1d3]/velscale ; Limits for centroid 
    parinfo[i+2].limits = [0,8d2]/velscale ; Limits for sigma_pix of narrow lines 
    parinfo[[1,2]+i].step = 1d-2 ; centroid/sigma step in pixel
    parinfo[i+3].limits = [-0.3, 0.3] ; h3
    parinfo[i+4].limits = [-0.3, 0.3] ; h4
    if gmoments eq 2 then begin
    	parinfo[i+[3,4]].fixed = [1,1]
	start_pars[i+[3,4]] = [0,0]
    endif
    ; for lines falling on bad pixels
    ind = where(goodpixels gt parinfo[i+1].limits[0] and $
    	goodpixels lt parinfo[i+1].limits[1], ngood)
    if ngood lt 1d3/velscale then begin
    	start_pars[i+[0,2,3,4]] = 0 ; set flux,sigma,h3,h4 to 0
	parinfo[i+[0,1,2,3,4]].fixed = 1 ; fix all pars
    endif
endfor
; for broad lines limit 5000 > sigma > 800 km/s, and fix h3, h4 at 0
ind = where(strpos(emission_setup.name[i_lines],'_br') ge 0,ct)
if ct gt 0 then begin
    parinfo[ind*5+2].limits = [8d2,5d3]/velscale
    parinfo[ind*5+3].fixed = 1
    parinfo[ind*5+4].fixed = 1
    start_pars[ind*5+3] = 0
    start_pars[ind*5+4] = 0
endif

; iv) constrain the line fluxes to be
; non-negative, unless we have specified negative amplitudes in the
; emission-line setup. In this case we want non-positive values.
for i = 0,nlines-1 do begin
    j = i_lines[i]
    if (emission_setup.a[j] gt 0) then begin
        parinfo[i*5].limited = [1,0]
        parinfo[i*5].limits(0) = 0
    endif else begin
        parinfo[i*5].limited = [0,1]
        parinfo[i*5].limits(1) = 0
    endelse
    parinfo[i*5].step = 1d-2
endfor

; v) for the SSP templates
; kinematics: V*, sigma* (km/s), h3*, h4*
parinfo[5*nlines].limits = start_pars[5*nlines] + [-1d3,1d3]
parinfo[5*nlines+1].limits = [0,8d2]
parinfo[5*nlines+[0,1]].step = 1d-2*velscale
parinfo[5*nlines+2].limits = [-0.3,0.3]
parinfo[5*nlines+3].limits = [-0.3,0.3]
if smoments eq 2 then begin
	parinfo[5*nlines+[2,3]].fixed = [1,1]
	start_pars[5*nlines+[2,3]] = [0,0]
endif

; SSP weights
for i=0,nssp-1 do begin
    if start_pars[5*nlines+4+i] gt 0 then begin
	; force positive
	parinfo[5*nlines+4+i].limited = [1,0]
	parinfo[5*nlines+4+i].limits(0) = 0
	parinfo[5*nlines+4+i].step = 1d-2
    endif else begin ; force zero
	parinfo[5*nlines+4+i].fixed = 1
    endelse
endfor

; for the reddening parameters (if needed). These will follow the
; emission-line parameters in the parameter array. 
if npars gt 5*nlines+4+nssp then begin
    parinfo[nlines*5+4+nssp].limits = [0d,5d]
    parinfo[nlines*5+4+nssp].step = 1d-3
endif

; C) find the lines for which the kinematics needs to be tied
; to that of other lines and set the tied parinfo (fit = 't#')
for i = 0,nlines-1 do begin
    j = i_lines[i]    
    ; check if we have a 't' tag
    if (strmid(emission_setup.fit[j],0,1) eq 't') then begin
        ; find the reference line reference index, as given in the
        ; 'fit' tag of the emission setup structure
        k_refline = fix(strmid(emission_setup.fit[j],1)) 
        ; which correspond to the following position in emission setup
        ; structure that was passed to this function
        j_refline = where(emission_setup.i eq k_refline)
        if (j_refline eq -1) then $
          message,'Hey, you tied '+emission_setup.name[j]+' to a line you are not even fitting...'
        ; and to the following position in the list of the emission
        ; lines to fit
        i_refline = where(emission_setup.i[i_lines] eq k_refline)
        ; wavelengths of the line and the ref. line 
        l_line    = emission_setup.lambda[j]         & str_l_line    = string(l_line)
        l_refline = emission_setup.lambda[j_refline] & str_l_refline = string(l_refline)
        ; centroid are tied (in pixel of log10 sampled spectra)
	parinfo[5*i+1].tied = $
          strcompress('P['+string(5*i_refline+1)+']-alog10('+str_l_refline+'/'+str_l_line+')/')+ $
          strtrim(string(lstep_gal),2)
	; intrinsic sigma (in pixel) and h3 and h4 are also tied
        parinfo[5*i+2].tied =  strcompress('P['+string(5*i_refline+2)+']')        
        parinfo[5*i+3].tied =  strcompress('P['+string(5*i_refline+3)+']')        
        parinfo[5*i+4].tied =  strcompress('P['+string(5*i_refline+4)+']')        
    endif
endfor

; D) fix conflicts between parinfo boundaries and start_pars
; if the parameter is not fixed but limited.
for i=0,npars-1 do begin
        if (~parinfo[i].fixed and parinfo[i].limited[0] and $
		start_pars[i] lt parinfo[i].limits[0]) then $
		start_pars[i] = mean(parinfo[i].limits)
        if (~parinfo[i].fixed and parinfo[i].limited[1] and $
		start_pars[i] gt parinfo[i].limits[1]) then $
		start_pars[i] =  mean(parinfo[i].limits)
        if parinfo[i].fixed then $ ; for fixed pars, do not use limits
		parinfo[i].limited = [0,0]
endfor

; E) setting the FUNCTARGS structure for MPFIT
functargs = {GALAXY:galaxy, NOISE:noise, L0_GAL:l0_gal, LSTEP_GAL:lstep_gal,$
	templates:templates, EMISSION_SETUP:emission_setup, L0_TEMPL:l0_templ, $
        GOODPIXELS:goodpixels, INT_DISP:int_disp, velscale:velscale,$
	BIAS:bias}


END
