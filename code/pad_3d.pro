
;+
;
; NAME:
; pad_3d
;
; PURPOSE:
; Pad edges of 3D image "a" with "n" zeroes. So if a is array of size [30,30,30], and n=5,
; pad_3d(a,n) will create an array [40,40,40] with a positioned at [5:34,5:34,5:34], and zeroes elsewhere.
;
; CALLING SEQUENCE:
;   b=pad_3d(a,n)
;
; INPUTS:
;   a = 3-dimensional array
;
;   n = required size of margin
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


function pad_3d,a,n

sz=size(a,/dim)
ndim=(size(a))[0]
szout=sz+n*2
type=size(a,/type)
a2=make_array(dimension=szout,type=type)
arg=strarr(ndim)
arg[*]='n'
if ndim ge 2 then arg[0:ndim-2]=arg[0:ndim-2]+','
arg='a2['+totalstring(arg)+']=a'
void=execute(arg)
return,a2

end