
;+
;
; NAME:
; totalstring
;
; PURPOSE:
; Given vector of strings e.g. sin=['feed','the','big cat'], sout=totalstring(sin) returns a single string sin='feedthebig cat' 
;
; CALLING SEQUENCE:
;   sout=totalstring(sin)
;
; INPUTS:
;   sin = array of strings
;
; OUTPUTS:
;
;   sout = single string
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

function totalstring,s

sm=s[0]
for i=1,n_elements(s)-1 do sm=sm+s[i]

return,sm

end