FUNCTION SPFIT_GHTEMPL, EMISSION_SETUP=emission_setup, PARS=pars, NPIX=npix,$
	LSTEP_GAL=lstep_gal, INT_DISP_PIX=int_disp_pix

; Take the emission-setup structure and the input pars parameter array
; to make emission-line single or multi-Gaussian-Hermite (GH) templates.
;
; Note that emission_setup is used to determine the line type,
; wavelength, and the flux ratios. The kinematics stored in
; emission_setup are not used; instead, the kinematics are specified by pars
;
; 20150812 - fixed bug on the flux of satellite line f_sline

i_lines = where(emission_setup.kind eq 'l') ; index of non-doublet lines
nlines  = n_elements(i_lines)
n_pars = nlines*5 
if (n_pars ne n_elements(pars)) then message,'Hey, this is not the right emission-line parameter array'
; array that will contain the emission-line templates, including multiplets
gaus = dblarr(npix,nlines)

; a) First create the emission-line templates corresponding to each
; single line (like Hb) or to each main line of a multiplet (like
; [OIII]5007)
for i = 0,nlines-1 do begin
   ; Make Gaussian templates amplitudes specified by the input pars array
   gaus[*,i]=spfit_create_ghs(X=findgen(npix),PARS=[pars[5*i:5*i+4]],INT_DISP_PIX=int_disp_pix[i_lines[i]])
endfor

; b) Then find all the satellite lines belonging to multiplets
; (like [OIII]4959), and add them to the main emission-line template 
; we just created in a)
i_slines = where(strmid(emission_setup.kind,0,1) eq 'd')
n_slines = n_elements(i_slines)

if i_slines[0] ne -1 then begin
    ; loop over the satellite lines
    for i = 0,n_slines-1 do begin
        ; Current index in the emission-line setup structure for 
        ; the current satellite line (e.g. 1 for [OIII]4959)
        j = i_slines[i]
        ; Find the reference line index, as given in the "kind" tag of the
        ; emission setup structure, which points to the main line in the
        ; present multiplet (e.g. 2 for [OIII]5007 if kind was d2 for [OIII]4959)
        k_mline = fix(strmid(emission_setup.kind[j],1)) 
        ; which correspond to the following position in the emission setup
        ; structure that was passed to this function (e.g. still 2, if
        ; the indices of the lines in the emission setup start at 0 and
        ; increase by 1)
        j_mline = where(emission_setup.i eq k_mline)
        ; Get the wavelengths of both satellite and main lines
        ; and compute the offset (in pix) that we need to apply in order
        ; to correctly place the satellite line.     
        l_sline = emission_setup.lambda[j]         
        l_mline = emission_setup.lambda[j_mline] 
        ; offset in pixels of a log10 sampled spectrum 
	offset  = alog10(l_mline/l_sline)/lstep_gal
        ; Get the index in the array of the lines to fit, corresponding
        ; to the main line of the present multiplet, so that we can add the 
        ; satellite emission-line template in the right place in the
        ; gaussian templates array
        i_mline = where(emission_setup.i[i_lines] eq k_mline)
        ; Finally, create the satellite template, and add it to that of
        ; the corresponding main line of this multiplet.
	; emission_setup.a gives the flux ratio between satellite and main lines
        f_sline = emission_setup.a[j]*pars[i_mline*5]
        gaus_sline = spfit_create_ghs(X=findgen(npix), $
                     PARS=[f_sline,pars[i_mline*5+1]-offset,pars[i_mline*5+2:i_mline*5+4]], $
                     INT_DISP_PIX=int_disp_pix[j])
        gaus[*,i_mline] = gaus[*,i_mline] + gaus_sline 
    endfor
endif

return, gaus
END

