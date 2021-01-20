FUNCTION SPFIT_FITFUNC, pars, GALAXY=galaxy, NOISE=noise, L0_GAL=l0_gal, LSTEP_GAL=lstep_gal,$
	templates=templates, EMISSION_SETUP=emission_setup, L0_TEMPL=l0_templ, $
        GOODPIXELS=goodpixels, INT_DISP=int_disp, velscale=velscale, BIAS=bias, $
	BESTFIT=bestfit, EMISSION_TEMPLATES=emission_templates
; HF: This is the deviates function for MPFIT. 
compile_opt idl2, hidden

npix = n_elements(galaxy) ; # of data points
npars = n_elements(pars) ; # of parameters
; Count the number of single lines or multiplets
nlines = n_elements(where(emission_setup.kind eq 'l')) 
; Count the number of SSP templates
nssp = (size(templates))[2]

; start_pars include all parameters (5x nlines + 4+ntempl + E(B-V))
; for each line, we have [ampl, centroid_pix, sigma_pix, h3, h4]
; for the SSPs, we have [v*, sigma*, h3*, h4*, and weights]

; Gauss-Hermite templates convolved by velocity and sigma given in input pars
; instrumental broadening is also convolved
int_disp_pix = int_disp/velscale
eml_pars = pars[0:nlines*5-1]
gaus = spfit_ghtempl(EMISSION_SETUP=emission_setup,LSTEP_GAL=lstep_gal, NPIX=npix, PARS=eml_pars, $
                        INT_DISP_PIX=int_disp_pix)
eml_model = total(gaus,2)
; Output array containing the best matching emission-line templates
emission_templates = gaus

; SSP templates
kinstars = pars[nlines*5:nlines*5+3] ; V*,sigma*,h3*,h4*
wstars = pars[nlines*5+4:nlines*5+4+nssp-1] ; weights for each SSP
cstar = spfit_convolve_templates(templates,kinstars,velscale)
cont_model = cstar # wstars
; Here the size of the convolved SSPs is still different from the
; input galaxy spectrum. This difference can be removed by multiplying the
; reddening curve, which has the same size as the galaxy array.
; But we can add a line to make sure cont_model matches galaxy array.
cont_model = cont_model[0:npix-1]

; redden the SSP templates
; do not redden the emission line templates
if npars-5*nlines-4-nssp gt 0 then begin
	ebv = pars[nlines*5+4+nssp]
	Vstar =  kinstars[0] + (l0_gal-l0_templ)*299792.458d*alog(10.0d)
    	reddening_attenuation = DUST_CALZETTI(l0_gal,lstep_gal,npix,ebv,Vstar,/log10)
	cont_model *= reddening_attenuation
endif

; compute deviates, which is the output of the function called by MPFIT
bestfit = cont_model + eml_model
err = (galaxy[goodpixels]-bestfit[goodpixels])/noise[goodpixels]

; penalize the deviates
tmp = total(pars[3:4]^2) + total(kinstars[2:3]^2) ; h3^2+h4^2
err = err + bias*robust_sigma(err,/zero)*sqrt(tmp)

return, err
END

