;+
;
; NAME:
; get_rev_ind
;
; PURPOSE:
; Given reverse indices as returned by histogram or hist_nd, and index
; of required histogram bin or bins, returns original points contributing to the histogram
; at those bins.
;
; CALLING SEQUENCE:
; index=get_rev_ind(revind,index_hist,n)
;
; INPUTS:
;   revind = reverse indices as returned by histogram function
;   index_hist = the histogram bin index for which the original contributing indices are required
;   n = the number of original contributing points that contribute to the current histogram bin(s).
;
;
; OUTPUTS:
;   Vector of index subscripts 
;
; OPTIONAL KEYWORD
;   none
;
; PROCEDURE:
;   For each member of index_hist, builds up index of contributing points
;   
; EXAMPLE:
;   array=randomn(seed,1001)
;   h=histogram(array,binsize=0.1,locations=xh,reverse_indices=revind)
;   void=min(abs(xh-0.1),indh);find histogram bin closest to 0.1
;   index=get_rev_ind(revind,indh,n)
;   print,'Number of array points that contribute to histogram bin = ',n
;   print,'Range of input values within this histogram bin = ',min(array[index]),' to ',max(array[index])
;
; USE & PERMISSIONS
;  If you reuse in your own code, please include acknowledgment to Huw Morgan (see below)
;
; ACKNOWLEDGMENTS:
;  This code was developed with the financial support of:
;  STFC Consolidated grant to Aberystwyth University (Morgan)
;
; MODIFICATION HISTORY:
; Created at Aberystwyth University 07/2019 - Huw Morgan hmorgan@aber.ac.uk
;
;
;-



function get_rev_ind,revind,index_hist,n

res=[-1ull]
for i=0l,n_elements(index_hist)-1 do begin
  res=[res,revind[revind[index_hist[i]]:revind[index_hist[i]+1]-1]]
endfor
res=res[1:*]
res=res[sort(res)]
n=n_elements(res)
return,res

end