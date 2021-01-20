;######################################################################
;+
; NAME:
;     DISPLAY_BINS
;
; AUTHOR:
;       Michele Cappellari, University of Oxford
;       cappellari_at_astro.ox.ac.uk
;
; PURPOSE:
;       Plot a Voronoi 2D-binned field, given a set of coordinates (XBIN,YBIN)
;       for the generators of the Voronoi tessellation, the corresponding values
;       to display for each bin, and the coordinates (X,Y) of the original unbinned
;       pixels. This routine is meant to visualize the output from the
;       VORONOI_2D_BINNING program and uses the same notations.
;
; CALLING SEQUENCE:
;       DISPLAY_BINS, xBin, yBin, velBin, x, y, $
;               PA=pa, PIXELSIZE=pixelSize, RANGE=[minVel,maxVel]
;
; INPUTS:
;       XBIN: Vector containing the X coordinate of the generator of each bin.
;           Note that the *generator* has to be used here, not the centroid!
;       YBIN: same as XBIN for the Y coordinate of each bin.
;       VELBIN: Vector containing the quantity to plot (e.g. velocity)
;           associated to each bin of coordinates XBIN and YBIN.
;       X: Vector containing the X coordinate of the original unbinned
;           pixels, from which the bins where constructed.
;       Y: Same as X, but for the Y coordinate.
;
; KEYWORDS:
;       PA: position agle of the side of each pixel.
;       PIXELSIZE: size of each pixel, in the same units as the coordinates.
;           If this is not given then its value is computed automatically.
;       RANGE: two elements vector [minVel,maxVel] defining the range of
;           VEL values to plot. Any value outside this range will be
;           visualized with the maximum or minimum color in the colormap.
;           By default RANGE=[min(VEL),max(VEL)].
;       _EXTRA: any additional keyword can be passed to PLOT via the _EXTRA
;           mechanism (e.g. TITLE, XTITLE, XRANGE,...).
;
; PROCEDURES REQUIRED:
;       This procedure requires the routine DISPLAY_PIXELS to work.
;
; MODIFICATION HISTORY:
;       V1.0: Written by Michele Cappellari, Leiden, 18 February 2003
;       V1.1: Documented and included error checking. MC, Leiden, 30 September 2005
;-
;----------------------------------------------------------------------------
pro display_bins, xbin, ybin, vel, x, y, out, _EXTRA=ex
compile_opt idl2
on_error, 2

nbin = N_ELEMENTS(vel)
if nbin ne n_elements(xbin) or nbin ne n_elements(ybin) then $
    message, 'The vectors (XBIN, YBIN, VEL) must have the same size'
npix = N_ELEMENTS(x)
if npix ne n_elements(y) then message, 'The vectors (X, Y) must have the same size'
if npix lt nbin then message, 'The vectors (X, Y) cannot be smaller than (XBIN, YBIN)'

; Perform a Voronoi tessellation starting from the coordinates
; of the generators and the coordinates of the original pixels
;
out = fltarr(npix,/NOZERO)
FOR j=0,npix-1 DO BEGIN
    tmp = min((x[j]-xbin)^2+(y[j]-ybin)^2,index)
    out[j] = vel[index]
ENDFOR
display_pixels, x, y, out, _EXTRA=ex

end
;----------------------------------------------------------------------
