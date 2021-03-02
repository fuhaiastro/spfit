;+
;
; NAME:
; ZERO_STRUCT
;
; PURPOSE:
; "Zero" all the elements in a structure. Numbers are set to zero,
; strings to '', pointers and objects references to NULL
;
; CALLING SEQUENCE:
; zero_struct, struct, idx
;
; INPUTS:
; struct: Structure to be zeroed. Can be an array of structures.
;
; REVISION HISTORY:
; Created 26-OCT-2000 Erin Scott Sheldon
; Better way using struct_assign. 2006-July-21 E.S.
;-
;
;
;
; Copyright (C) 2005 Erin Sheldon, NYU. erin dot sheldon at gmail dot com
;
; This program is free software; you can redistribute it and/or modify
; it under the terms of the GNU General Public License as published by
; the Free Software Foundation; either version 2 of the License, or
; (at your option) any later version.
;
; This program is distributed in the hope that it will be useful,
; but WITHOUT ANY WARRANTY; without even the implied warranty of
; MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
; GNU General Public License for more details.
;
; You should have received a copy of the GNU General Public License
; along with this program; if not, write to the Free Software
; Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA
;
;


PRO zero_struct, struct, idx

IF N_params() EQ 0 THEN BEGIN
print,'-Syntax: zero_struct, struct, [idx]'
print,''
print,'Use doc_library,"zero_struct" for more help.'
return
ENDIF

IF size(struct,/tname) NE 'STRUCT' THEN BEGIN
message,'Input value must be a structure. Input variable unchanged',/inf
return
ENDIF

;; Make a structure with a random variable name and use struct_assign
;; to "copy" between structures. By default, fields that do not match
;; are zeroed

tagname = 'randomFront'
numstr = $
strtrim(string(ulong64(1000*systime(1))), 2) + 'moreRandom' + $
strtrim(string(long(1000000*randomu(seed))), 2)

tagname = tagname + numstr
cst = create_struct(tagname, 0)

if n_elements(idx) gt 0 then begin
	if idx[0] eq -1 then begin
		print,'Warning: no elements to be zeroed!'
		return
	endif
	struct2 = struct[idx]
	struct_assign, cst, struct2
	struct[idx] = struct2
endif else begin
	struct_assign, cst, struct
endelse

return
END


