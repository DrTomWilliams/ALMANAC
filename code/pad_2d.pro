
;+
;
; NAME:
; pad_2d
;
; PURPOSE:
; Pad edges of 2D image "a" with "n" zeroes. So if a is array of size [30,30], and n=5,
; pad_2d(a,n) will create an array [40,40] with a positioned at [5:34,5:34], and zeroes elsewhere.
;
; CALLING SEQUENCE:
;   b=pad_2d(a,n)
;
; INPUTS:
;   a = 2-dimensional array
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

function pad_2d,a,n,reflect=reflect,ramp=ramp

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

if keyword_set(reflect) then begin
  a2[0:n-1,n:n+sz[1]-1]=reverse(a[0:n-1,*],1)
  a2[n+sz[0]:*,n:n+sz[1]-1]=reverse(a[sz[0]-n:*,*],1)
  a2[n:n+sz[0]-1,0:n-1]=reverse(a[*,0:n-1],2)
  a2[n:n+sz[0]-1,n+sz[1]:*]=reverse(a[*,sz[1]-n:*],2)
  a2[0:n-1,0:n-1]=rotate(a[0:n-1,0:n-1],2)
  a2[0:n-1,n+sz[1]:*]=rotate(a[0:n-1,sz[1]-n:*],2)
  a2[n+sz[0]:*,0:n-1]=rotate(a[sz[0]-n:*,0:n-1],2)
  a2[n+sz[0]:*,n+sz[1]:*]=rotate(a[sz[0]-n:*,sz[1]-n:*],2)
  if keyword_set(ramp) then begin
    r=findgen(n)/(n-1)
    r=rebin(reform(r,1,n),sz[0]+n*2,n)
    a2[*,0:n-1]=a2[*,0:n-1]*r
    a2[*,n+sz[1]:*]=a2[*,n+sz[1]:*]*rotate(r,7)
    a2[0:n-1,*]=a2[0:n-1,*]*rotate(r,4)
    a2[n+sz[0]:*,*]=a2[n+sz[0]:*,*]*rotate(r,1)
  endif
endif

return,a2

end