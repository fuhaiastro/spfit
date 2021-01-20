function boxave3d, array, xbin, ybin, zbin
;+
; Name
;	boxave3d
;
; Purpose
; 	Properly handle !NaNs as missing data with mean() and dimensional
; 	juggling (http://www.idlcoyote.com/code_tips/rebinmissing.php)
; 
; Similar IDL Functions: 
;	HBOXAVE, HCONGRID, HREBIN 
; 
; To compare with rebin:
;	IDL > a = findgen(10,10,3)        
;	IDL > a[0,0,0] = !values.f_nan
;	IDL > a_boxave = boxave3d(a,2,2,1)
;	IDL > a_rebin = rebin(a,5,5,3)    
;	IDL > print,a_boxave[*,*,0],a_rebin[*,*,0]
;
; 2015 May, Hai Fu, UIowa
;-
; cause control to be returned to the caller of a procedure in the event of an error.
ON_ERROR, 2  

dim = (size(array))[0:3]
if xbin*(dim[1]/xbin)*ybin*(dim[2]/ybin)*zbin*(dim[3]/zbin) ne n_elements(array) then $
	message,'The supplied bins must be integral factors of the original dimensions.'
array1 = reform(array, xbin, dim[1]/xbin, ybin, dim[2]/ybin, zbin, dim[3]/zbin)
array2 = transpose(array1,[0,2,4,1,3,5])
array3 = reform(array2,xbin*ybin*zbin,dim[1]/xbin,dim[2]/ybin,dim[3]/zbin)
resampled = mean(array3,dim=1,/nan)
return,resampled

end
