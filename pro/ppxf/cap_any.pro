;----------------------------------------------------------------------
function cap_any, a
on_error, 2
;
; Emulate MATLAB function any()
; V1.0: M. Cappellari, Leiden, December 2002
; V1.1: Uses logical NOT operator. MC, Oxford, 23 October 2007
; V1.2: Compare to BYTE type. MC, Windhoek, 3 October 2008
; V1.21: Renamed CAP_ANY to avoid potential naming conflicts.
;    MC, Paranal, 8 November 2013

    return, ~array_equal(a,0b)
end
;----------------------------------------------------------------------
