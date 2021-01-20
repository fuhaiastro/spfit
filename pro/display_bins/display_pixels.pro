;######################################################################
;+
; NAME:
;     DISPLAY_PIXELS
;
; AUTHOR:
;       Michele Cappellari, University of Oxford
;       cappellari_at_astro.ox.ac.uk
;
; PURPOSE:
;       Plot a set of square pixels, given a set of (X,Y) coordinates
;       and of corresponding values to display. This is uesful when the pixels
;       are not arranged in an full array and cannot be displayed as an image.
;       This is often the case for integral-field spectroscopic data.
;
; CALLING SEQUENCE:
;       DISPLAY_PIXELS, x, y, vel, PA=pa, PIXELSIZE=pixelSize, RANGE=[minVel,maxVel]
;
; INPUTS:
;       X: Vector containing the X coordinate of the center of each pixel.
;       Y: same as X for the Y coordinate of each pixel.
;       VEL: Vector containing the quantity to plot (e.g. velocity)
;           associated to each pixel of coordinates (X,Y).
;
; KEYWORDS:
;       PA: position agle of the side of each pixel.
;       PIXELSIZE: size of each pixel, in the same units as the coordinates.
;           If this is not given then its value is computed automatically.
;       RANGE: two elements vector [minVel,maxVel] defining the range of
;           VEL values to plot. Any value outside this range will be
;           visualized with the maximum or minimum color in the colormap.
;           By default RANGE=[min(VEL),max(VEL)].
;       FLUX: vector of the same size as VEL, containing the flux associated
;           to each pixel. Overplot the isophotes on top of the field, spaced
;           by one magnitude intervals.
;       _EXTRA: any additional keyword can be passed to PLOT via the _EXTRA
;           mechanism (e.g. TITLE, XTITLE, XRANGE,...).
;
; MODIFICATION HISTORY:
;       V1.0: Michele Cappellari, Leiden, January 2000
;       V1.1: added RANGE keyword and cuts for numbers out of range.
;           MC, Leiden, 3 December 2002
;       V1.2: Added more input checking. MC, Leiden, 12 December 2004
;       V1.21: Shortened loop. MC, Leiden, 10 September 2005
;       V1.22: Perform robust but slow pixel size calculation by default.
;           Prevent tick marks to be covered by pixels. Updated documentation.
;           MC, Leiden, 29 September 2005
;       V1.23: added _EXTRA keyword to PLOT call to allow for the
;           TITLE keyword to be passed. MC, Leiden, 9 October 2005
;       V1.24: added FLUX keyword. MC, Leiden, 4 December 2005
;       V1.25: Use intrinsic functions BYTSCL to scale colors and 
;           DISTANCE_MEASURE for pixelSize. MC, Oxford, 10 November 2011
;-
;----------------------------------------------------------------------------
pro display_pixels, x, y, vel, $
    PA=pa, PIXELSIZE=pixelSize, RANGE=range, FLUX=flux, _EXTRA=ex
compile_opt idl2
on_error, 2

n = N_ELEMENTS(vel)
if n ne n_elements(x) or n ne n_elements(y) then $
    message, 'The vectors (X, Y, VEL) must have the same size'

; For each point, find the distance to all other points and select the minimum.
; This is a robust but slow way of determining the pixel size of unbinned data.
;
if n_elements(pixelSize) eq 0 then $
    pixelSize = min(distance_measure( transpose([[x],[y]]) ))
if n_elements(pa) eq 0 then begin
    cpa = 1.0
    spa = 0.0
endif else begin
    cpa = COS(pa/!RADEG)
    spa = SIN(pa/!RADEG)
endelse
if n_elements(range) eq 0 then $
    color = bytscl(vel) $
else $
    color = bytscl(vel, MIN=range[0], MAX=range[1])

maxx = max(x, MIN=minx)
maxy = max(y, MIN=miny)
xaper  = [-0.5, -0.5, +0.5, +0.5, -0.5] * pixelSize
yaper  = [+0.5, -0.5, -0.5, +0.5, +0.5] * pixelSize
x1 = xaper*cpa - yaper*spa
y1 = xaper*spa + yaper*cpa

; Plot an empty frame first and then overplot the axes on top of it
; to prevent the axes tick marks from being covered by the pixels.
;
dx = pixelSize/sqrt(2) ; Make sure pixels always lie inside the plot
PLOT, [minx-dx, maxx+dx], [miny-dx, maxy+dx], /NODATA, /ISO, XSTYLE=5, YSTYLE=5, _EXTRA=ex
FOR j=0L, n-1L DO POLYFILL, x[j]+x1, y[j]+y1, COLOR=color[j]
AXIS, XAXIS=0, XRANGE=[minx-dx, maxx+dx], /XSTYLE, XTITLE='arcsec', _EXTRA=ex
AXIS, XAXIS=1, XRANGE=[minx-dx, maxx+dx], /XSTYLE, XTICKFORMAT='(A1)'
AXIS, YAXIS=0, XRANGE=[miny-dx, maxy+dx], /YSTYLE, YTITLE='arcsec', _EXTRA=ex
AXIS, YAXIS=1, XRANGE=[miny-dx, maxy+dx], /YSTYLE, YTICKFORMAT='(A1)'

nf = n_elements(flux)
if nf gt 0 then begin
    if nf ne n then message, 'The vector FLUX must have the same size as (X, Y)'
    contour, 2.5*alog10(flux/max(flux)), x, y, /OVERPLOT, /IRREGULAR, $
        LEVELS=-reverse(findgen(10)), COLOR=0 ; Contours spaced by 1 mag
endif

END
;----------------------------------------------------------------------
