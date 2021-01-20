;######################################################################
;+
; NAME:
;     PLOT_VELFIELD
;
; AUTHOR:
;       Michele Cappellari, University of Oxford
;       cappellari_at_astro.ox.ac.uk
;
; PURPOSE:
;       Plot a linearly interpolated field, given a set of (X,Y) coordinates
;       and of corresponding values to display. This routine may be used
;       to visualize the output from the VORONOI_2D_BINNING program.
;
; CALLING SEQUENCE:
;       PLOT_VELFIELD, xBin, yBin, velBin, $
;               FLUX=fluxBin, NCOLORS=nColors, /NODOTS, RANGE=[minVel,maxVel]
;
; INPUTS:
;       XBIN: Vector containing the X coordinate of the centroid of each bin.
;           The bin *centroid* should be preferably used here for plotting
;           instead of the bin generator, although the two values are generally
;           quite similar.
;       YBIN: same as XBIN for the Y coordinate of each bin.
;       VELBIN: Vector containing the quantity to plot (e.g. velocity)
;           associated to each bin of coordinates XBIN and YBIN.
;
; KEYWORDS:
;       FLUX: vector with the same size as XBIN and YBIN, containing the
;           average surface brightness (in linear units) of each bin.
;           If this vector is passed, it will be overplotted as
;           contours spaced by 1 mag intervals.
;       NCOLORS: number of colors to display (default 128, maximum 256).
;       /NODOTS: set this keyword to prevent the coordinates of the
;           bin centroids from being overplotted to the contour plot.
;       RANGE: two elements vector [velMin,velMax] defining the range of
;           VELBIN values to plot. Any value outside this range will
;           be visualized with the maximum or minimum color in the colormap.
;       _EXTRA: any additional keyword can be passed to plot via the _EXTRA
;           mechanism (e.g. TITLE, XRANGE,...).
;
; MODIFICATION HISTORY:
;       V1.0: Michele Cappellari, Vicenza, 7 December 2002
;       V1.01: Overplot contours in magnitudes. MC, Leiden, 1 August 2003
;       V1.02: Written documentation. MC, Leiden, 4 December 2004
;       V1.03: Added _EXTRA keyword to OPLOT. MC, Vicenza, 30 December 2010
;-
;----------------------------------------------------------------------------
pro plot_velfield, xBin, yBin, velBin, $
     FLUX=flux, ISO=iso, NCOLORS=ncolors, NODOTS=nodots, RANGE=range, _EXTRA=extra
compile_opt idl2
on_error, 2

nx = n_elements(xBin)
if nx ne n_elements(yBin) or nx ne n_elements(velBin) then $
    message, 'The vectors (XBIN, YBIN, VELBIN) must have the same size'

if n_elements(range) eq 0 then begin
    mx = max(velBin, MIN=mn)
    range = [mn,mx]
endif
if n_elements(ncolors) eq 0 then ncolors = 128 else ncolors = ncolors<256
if n_elements(iso) eq 0 then iso=1 else iso=0

levels = range[0] + (range[1]-range[0])/(ncolors-1.0)*findgen(ncolors)
contour, velBin>range[0]<range[1], xBin, yBin, /FILL, ISO=iso, LEVELS=levels, $
    /XSTYLE, /YSTYLE, /IRREGULAR, XTITLE='arcsec', YTITLE='arcsec', _EXTRA=extra
if not keyword_set(nodots) then oplot, xBin, yBin, PSYM=3, COLOR=0, _EXTRA=extra

nf = n_elements(flux)
if nf gt 0 then begin
    if nf ne nx then $
        message, 'The vector FLUX must have the same size as (XBIN, YBIN)'
    contour, 2.5*alog10(flux/max(flux)), xbin, ybin, /OVERPLOT, /IRREGULAR, $
        LEVELS=-reverse(findgen(10)), COLOR=0 ; Contours spaced by 1 mag
endif

end
;----------------------------------------------------------------------------
