


;+
;
; NAME:
; unpad_3d
;
; PURPOSE:
; Remove, or unpad, a margin of n pixels from edges of 3D datacube. So if a is array of size [40,40,40], and n=5,
; b=unpad_3d(a,n) will create an array [30,30,30] with b extracted from a[5:34,5:34,5:34]. Useful reversal of pad_3d.
;
; CALLING SEQUENCE:
;   b=unpad_3d(a,n)
;
; INPUTS:
;   a = 3-dimensional array
;
;   n = required size of margin
;
; OUTPUTS:
;
;   b = unpadded array
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

function int2str_huw2,num,width=width
  return,string(num,format=n_elements(width) eq 0?'(I0)':$
    '(I' + int2str(width) + '.' + int2str(width) +')')
end

function int2str_huw,num,width
  return,string(num,format=n_params() lt 2?'(I0)':$
    '(I' + int2str_huw2(width) + '.' + int2str_huw2(width) +')')
end

function unpad_nd,a,n

  sz=size(a,/dim)-n*2
  ndim=(size(a))[0]
  arg=strarr(ndim)
  arg[*]='n:n+sz['+int2str_huw(indgen(ndim),1)+']-1'
  if ndim ge 2 then arg[0:ndim-2]=arg[0:ndim-2]+','
  arg='a=a['+totalstring(arg)+']'
  void=execute(arg)
  return,a

end