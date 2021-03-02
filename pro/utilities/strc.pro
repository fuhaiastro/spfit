function Strc, thing

;+
; removes the extra spaces and unnecessary zeros
; will choke on arrays, only works on scalars
; June 1994, M. Liu (UCB)
;
; added chopping of unneeded zeros on the end of floats/doubles
; (beware if numbers are too long, will be chopped off - this is
; and IDL feature of the 'string' command, not due to this function)
; 11/05/96 MCL
;-

if n_params() eq 0 then return, ''

sz = size(thing)
typ = sz(sz(0)+1)

if sz(0) eq 0 then begin 

    ; round any floats of their excess zeros,
    ; putting a zero at the end if there's an exposed decimal pt
    ; -> only works for scalars, not vectors
    ; -> this screws up sci notation, e.g. '1e10' becomes '1e1'
    typ = sz(0)+1
    ss = strcompress(string(thing), /remove_all)
    if (sz(typ) eq 5) or (sz(typ) eq 4) then begin
        while (strmid(ss,strlen(ss)-1,1) eq '0') do begin
            ss = strmid(ss,0,strlen(ss)-1)
        endwhile
        if (strmid(ss,strlen(ss)-1,1) eq '.') then ss = ss+'0'
    endif 

endif else $
  ss =  strcompress(thing, /remove_all)



return,ss   
;return, strcompress(string(thing), /remove_all) 


end
