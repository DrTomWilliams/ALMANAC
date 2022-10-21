function minmax_huw,a,nan=nan


mn=min(a,max=mx,nan=nan)

return,[mn,mx]

end

function gamma_transform,imin,gamma,a0,a1

if n_params() lt 2 then gamma=3.0
if n_params() lt 3 then a0=min(imin,/nan)
if n_params() lt 4 then a1=max(imin,/nan)

return,((float(imin)-a0)/(a1-a0))^(1/float(gamma))

end

;; $Id: //depot/Release/ENVI50_IDL82/idl/idldir/lib/gaussian_function.pro#1 $
;;
;; Copyright (c) 2010-2012, Exelis Visual Information Solutions, Inc. All
;       rights reserved. Unauthorized reproduction is prohibited.
;+
;; Gaussian_Function
;;
;; Purpose:
;;   This function returns a Gaussian kernel
;;
;-

;;---------------------------------------------------------------------------
;; Gaussian_Function
;;
;; Purpose:
;;   Create a guassian kernel
;;
;; Parameters:
;;   SIGMA - sigma value
;;
;; Keywords:
;;   DOUBLE - If set, use double precision
;;
;;   MAXIMUM - Set this keyword to the value to be used as the maximum
;;             value of the resulting array
;;
;;   NORMALIZE - If this keyword is set the peak height shall be
;;               determined such that the total of the Gaussian is 1.0
;;
;;   WIDTH - desired width of the array
;;
FUNCTION gaussian_function, sigmaIn, DOUBLE=doubleIn, MAXIMUM=maxIn, $
                                     NORMALIZE=normalIn, WIDTH=widthIn
  compile_opt hidden, idl2
on_error, 2                                     
    
  double = KEYWORD_SET(doubleIn)
  sigma = (N_ELEMENTS(sigmaIn) eq 0) ? 1.0 : $
    FIX(sigmaIn[*] > 0.00001d, TYPE=double+4)
  maximumSet = KEYWORD_SET(maxIn)
  maximum = (N_ELEMENTS(maxIn) eq 0) ? 1.0 : FIX(maxIn[0], TYPE=double+4)
  normal = KEYWORD_SET(normalIn)
  ;; fill in width if not provided
  if (N_ELEMENTS(widthIn) eq 0) then begin
    width = CEIL(sigma*3)
    ;; ensure width is odd
    width or= 1
    width *= 2
    width++
  endif else begin
    width = widthIn[*]
  endelse

  width = FIX(width) > 1
  nSigma = N_ELEMENTS(sigma)
  if (nSigma gt 8) then begin
    message, 'Sigma can have no more than 8 elements'
    return, 0
  endif
  nWidth = N_ELEMENTS(width)
  if (nWidth gt nSigma) then begin
    message, 'Incorrect width specification'
    return, 0
  endif
  if ((nWidth eq 1) && (nSigma gt 1)) then $
    width = REPLICATE(width, nSigma)

  kernel = replicate((keyword_set(double) ? 0.0d : 0.0), width)

  ;; Fill in all 8 dimensions
  temp = intarr(8)
  temp[0] = width
  width = temp  
    
  a = (b = (c = (d = (e = (f = (g = (h = 0)))))))  
  ;; create indices
  switch nSigma of
    8 : h = (keyword_set(double) ? $
      dindgen(width[7])-width[7]/2 : $
      findgen(width[7])-width[7]/2) + $
      (width[7] and 1 ? 0 : 0.5)
    7 : g = keyword_set(double) ? $
      dindgen(width[6])-width[6]/2 : $
      findgen(width[6])-width[6]/2 + $
      (width[6] and 1 ? 0 : 0.5)
    6 : f = keyword_set(double) ? $
      dindgen(width[5])-width[5]/2 : $
      findgen(width[5])-width[5]/2 + $
      (width[5] and 1 ? 0 : 0.5)
    5 : e = keyword_set(double) ? $
      dindgen(width[4])-width[4]/2 : $
      findgen(width[4])-width[4]/2 + $
      (width[4] and 1 ? 0 : 0.5)
    4 : d = keyword_set(double) ? $
      dindgen(width[3])-width[3]/2 : $
      findgen(width[3])-width[3]/2 + $
      (width[3] and 1 ? 0 : 0.5)
    3 : c = keyword_set(double) ? $
      dindgen(width[2])-width[2]/2 : $
      findgen(width[2])-width[2]/2 + $
      (width[2] and 1 ? 0 : 0.5)
    2 : b = keyword_set(double) ? $
      dindgen(width[1])-width[1]/2 : $
      findgen(width[1])-width[1]/2 + $
      (width[1] and 1 ? 0 : 0.5)
    1 : a = keyword_set(double) ? $
      dindgen(width[0])-width[0]/2 : $
      findgen(width[0])-width[0]/2 + $
      (width[0] and 1 ? 0 : 0.5)
  endswitch
  
  ;; create kernel
  for hh=0,width[7]-1>0 do $
    for gg=0,width[6]-1>0 do $
      for ff=0,width[5]-1>0 do $
        for ee=0,width[4]-1>0 do $
          for dd=0,width[3]-1>0 do $
            for cc=0,width[2]-1>0 do $
              for bb=0,width[1]-1>0 do $
                for aa=0,width[0]-1>0 do $
                  kernel[aa,bb,cc,dd,ee,ff,gg,hh] = $
                    exp(-((a[aa]^2)/(2*sigma[[0]]^2) + $
                          (b[bb]^2)/(2*sigma[[1]]^2) + $
                          (c[cc]^2)/(2*sigma[[2]]^2) + $
                          (d[dd]^2)/(2*sigma[[3]]^2) + $
                          (e[ee]^2)/(2*sigma[[4]]^2) + $
                          (f[ff]^2)/(2*sigma[[5]]^2) + $
                          (g[gg]^2)/(2*sigma[[6]]^2) + $
                          (h[hh]^2)/(2*sigma[[7]]^2)))
  
  if (KEYWORD_SET(maximumSet)) then begin
    kernel *= maximum
  endif else begin
    if (KEYWORD_SET(normal)) then begin
      kernel /= total(kernel, /PRESERVE_TYPE)
    endif
  endelse
  
  return, kernel

end


;hmorgan@aber.ac.uk 09/2013
;imin - input image
;a0 - minimum input value for gamma transform
;a1- maximum input value for gamma transform
;gamma - default is 3.5. Used for gamma transform
;h - optional weighting for combining of gamma-transformed image and MGN image. Default 0.9
;k - optional contrast stretching of MGN images. Default value 1


function mgn,imin,a0,a1,gamma=gamma,h=h,k=k,imp=imp
            

gamma=keyword_set(gamma)?gamma:3.5
h=keyword_set(h)?h:0.9
k=keyword_set(k)?k:1;0.6
imin2=float(imin);cumulative_clip(float(imin),[0.001,99.999],/preserve)
if n_elements(a0) eq 0 then a0=min(imin2,/nan)
if n_elements(a1) eq 0 then a1=max(imin2,/nan) 
im=gamma_transform(imin2,gamma,a0,a1)

;width of Gaussian kernels
w=[2.5,5,10,20,40,80]

nw=n_elements(w)
s=size(im)
imp=fltarr(s[1],s[2]);,nw)

;loop through Gaussian kernel widths (different spatial scales)
for iw=0,nw-1 do begin
  ;print,iw
  ker=gaussian_function(w[iw]/4.,/norm)
  ind=minmax_huw(where(ker ge max(ker)*0.05))
  ker=ker[ind[0]:ind[1]]
  if n_elements(ker) ge s[1] or n_elements(ker) ge s[2] then begin
    imp[*,*,iw]=!values.f_nan
    continue
  endif
  norm=total(ker)
  nk=n_elements(ker)
  
  ;calculate local mean
  m=rotate(convol(rotate(convol(imin2,ker,norm,/edge_trunc,/nan),4),ker,norm,/edge_trunc,/nan),4)
  
  ;calculate local standard deviation
  md=rotate(sqrt(convol(rotate(convol((imin2-m)^2,ker,norm,/edge_trunc,/nan),4), $
          ker,norm,/edge_trunc,/nan)),4)
  
  const=mean(md,/nan)/10;this constant factor included in denominator below
                          ;helps reduce noise in final image
  
  ;add current result to processed image. Subtract mean, divide by standard deviation
  imp=imp+(imin2-m)/(md+const)

endfor

imp=temporary(imp)/nw
imp=atan(k*temporary(imp));apply arctangent transform
imout=h*im+(1-h)*imp;final image combination of gamma-transformed and normalized


return,imout

end
