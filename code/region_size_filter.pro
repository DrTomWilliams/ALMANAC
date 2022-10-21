function region_size_filter,image,threshold,npix,label=label,mask=mask, $
            max_region_size=max_region_size,all_neighbors=all_neighbors,return_mask=return_mask, $
            smooth=smooth,npixinitial=npixinitial

sizearr,image,nx,ny,nz,ndim=ndim
if ndim gt 3 then begin
  print,'Region size filter only implemented for 2 or 3 dimensions, returning -1'
  return,-1
endif

case ndim of
  1:stop
  2:mask=pad_2d(image,1) gt threshold
  3:mask=pad_3d(image,1) gt threshold
endcase

if keyword_set(npixinitial) then begin
  label=label_region(mask,/ulong,all_neighbors=all_neighbors)
  h=histogram(label,min=1,binsize=1,loc=x,rev=ri)
  ind=where(h lt npixinitial,cnt)
  if cnt gt 0 then begin
    ind=get_rev_ind(ri,ind)
    mask[ind]=0
  endif
endif

if keyword_set(smooth) then begin
  nsmo=n_elements(smooth) lt ndim?fltarr(ndim)+smooth[0]:smooth
  case ndim of
    1:stop
    2:begin
      kx=gaussian_function(nsmo[0],/norm)
      ky=gaussian_function(nsmo[1],/norm)
      ky=reform(ky,1,n_elements(ky))
    end
    3:begin
      kx=gaussian_function(nsmo[0],/norm)
      ky=gaussian_function(nsmo[1],/norm)
      ky=reform(ky,1,n_elements(ky),1)
      kz=gaussian_function(nsmo[2],/norm)
      kz=reform(kz,1,1,n_elements(kz))
    end
  endcase   
  masksmo=convol(float(mask),kx,/edge_trunc,/nan)
  if ndim gt 1 then masksmo=convol(masksmo,ky,/edge_trunc,/nan)
  if ndim gt 2 then masksmo=convol(masksmo,kz,/edge_trunc,/nan)
  masksmo=masksmo gt 1.e-2
endif else masksmo=mask;no smoothing

ind=where(masksmo,cntmaskinit)
if cntmaskinit eq 0 then return,image*0

label=label_region(masksmo,/ulong,all_neighbors=all_neighbors)
label=label*mask

h=histogram(label,min=1,binsize=1,loc=x,rev=ri)
ind=keyword_set(max_region_size)?where(h lt npix or h gt max_region_size,cnt):where(h lt npix,cnt)
if cnt gt 0 then begin
  ind=get_rev_ind(ri,ind)
  mask[ind]=0
endif 
label=label_region(mask,/ulong);redo label region, just to get ordered incremental set

case ndim of
  1:stop
  2:begin
    mask=unpad_2d(mask,1)
    label=unpad_2d(label,1)
  end
  3:begin
    mask=unpad_3d(mask,1)
    label=unpad_3d(label,1)
  end
endcase

return,keyword_set(return_mask)?mask:image*mask

end