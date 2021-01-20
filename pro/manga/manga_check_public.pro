;+
; Check if a list of MaNGA IDs are in public DR14 or latest internal release
;-
PRO manga_check_public,dir=dir,inlist=inlist

if keyword_set(dir) then begin
	inlist = exfilename(file_search(dir+'/*.png'),/noext)
	outfile = dir+'.txt'
endif else if keyword_set(inlist) then begin
	outfile = 'check_result.txt'
endif else begin
	print,'Either DIR or INLIST has to be set!'
	print,'e.g., manga_check_public,dir="classify/agn/AGN_Cone"'
	print,'e.g., manga_check_public,inlist=["7977-9101","8465-12705"]'
	return
endelse

openw,5,outfile

public = mrdfits('$MANGA_SPECTRO_REDUX/drpall-v2_1_2.fits',1,/silent)
intern = mrdfits('$MANGA_SPECTRO_REDUX/drpall-v2_2_0.fits',1,/silent)
pubver = public[0].versdrp3
intver = intern[0].versdrp3

fmt = '(a-10,1x,a-10,1x,a)'
for i=0,n_elements(inlist)-1 do begin
	ind1 = where(public.plateifu eq inlist[i],ct1)
	ind2 = where(intern.plateifu eq inlist[i],ct2)
	if ct1 gt 0 then printf,5,inlist[i],'public',pubver,f=fmt
	if ct1 eq 0 and ct2 gt 0 then printf,5,inlist[i],'internal',intver,f=fmt
	if ct1 eq 0 and ct2 eq 0 then printf,5,inlist[i],'not found',f=fmt
endfor

close,5

end

