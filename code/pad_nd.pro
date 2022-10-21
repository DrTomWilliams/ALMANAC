
;+
;
; NAME:
; pad_4d
;
; PURPOSE:
; Pad edges of ND array "A" with "n0" zeroes. So if A is array of size [30,30,30,30], and n0=5,
; pad_nd(A,n0) will create an array [40,40,40,40] with A positioned at [5:34,5:34,5:34,5:34], and zeroes elsewhere.
;
; CALLING SEQUENCE:
;   b=pad_nd(a,n0)
;
; INPUTS:
;   a = n-dimensional array
;
;   n0 = required size of margin
;
; OUTPUTS:
;   
;   b = padded array
;
; OPTIONAL KEYWORDS
;  (None)
;  
; OPTIONAL INPUTS
;   (NONE)
;
; OPTIONAL OUTPUTS
;   (none)
;
; PROCEDURE:
;
;
; USE & PERMISSIONS
; Any problems/queries, or suggestions for improvements, please email Huw Morgan, hmorgan@aber.ac.uk
;
; ACKNOWLEDGMENTS:
;  This code was developed with the financial support of:
;  STFC and Coleg Cymraeg Cenedlaethol Studentship to Aberystwyth University (Humphries)
;  STFC Consolidated grant to Aberystwyth University (Morgan)
;
; MODIFICATION HISTORY:
; Created at Aberystwyth University 2021 - Huw Morgan hmorgan@aber.ac.uk
;
;
;-


function pad_nd,a,n0

sz=size(a,/dim)
ndim=(size(a))[0]
szout=sz+n0*2
type=size(a,/type)
a2=make_array(dimension=szout,type=type)
arg=strarr(ndim)
arg[*]='n0'
if ndim ge 2 then arg[0:ndim-2]=arg[0:ndim-2]+','
arg='a2['+totalstring(arg)+']=a'
void=execute(arg)
return,a2

end