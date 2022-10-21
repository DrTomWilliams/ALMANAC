pro read_sdo_tom, files, index, data, xllp, yllp , nxp, nyp, $
    _extra=_extra, nodata=nodata, fnames_uncomp=fnames_uncomp, $
    mixed_comp=mixed_comp, parent_out=parent_out,outdir=outdir, $
    comp_header=comp_header, time_tag=time_tag, verbose=verbose, $
    count=count

version=2.1
if not file_exist(files(0)) then begin
   box_message,'IDL> read_sdo_tom,<filelist>,index [,data,llpx,llpy,nx,ny] [,/noshell] [/use_shared]
   return
endif

noshell=keyword_set(noshell)

if keyword_set(mixed_comp) then begin 
   nfiles=n_elements(files)
   fsize=file_size(files)
   aac=where(fsize lt 33569280,ccnt)
   if ccnt gt 0 and ccnt ne nfiles then begin 
      box_message,'Mixed compression'
      ifiles=files
      read_sdo_tom,files(aac),iizz,ddzz,/only_uncompress,fnames_uncomp=unames, $
          parent_out=parent_out, outdir=outdir, noshell=noshell
      files(aac)=unames
   endif ; else box_message,'/MIXED_COMP set but already homogenous'
endif
next=get_fits_nextend(files(0))
if next gt 0 then begin
  fits_info, files[0], extname=extname, /silent
  if stregex(/boo,/fold, extname[1],'desat_info') then next = 0
endif

proc=(['mreadfits_shm','mreadfits_tilecomp'])(next gt 0)

nodata=keyword_set(nodata)

use_index=required_tags(_extra,/use_index)
if use_index and data_chk(index,/struct) then orig_index=index ; save for .HISTORY update

case 1 of 
   n_params() lt 2: box_message,'IDL> read_sdo_tom,files,index [,data [,xll,yll,nx,ny]]
   n_params() eq 2 or nodata: mreadfits_header,files,index,exten=next,_extra=_extra
   n_params() eq 3: begin
      if next eq 0 then mreadfits_shm,files,index,data,_extra=_extra else $
         mreadfits_tilecomp,files,index,data,_extra=_extra, fnames_uncomp=fnames_uncomp, $
            parent_out=parent_out, outdir=outdir, time_tag=time_tag, /silent
   endcase
   else: begin
      if next eq 0 then $
         mreadfits_shm,files,index,data,xllp,yllp,nxp,nyp,_extra=_extra else $
         mreadfits_tilecomp,files,index,data,xllp,yllp,nxp,nyp,_extra=_extra, /silent, $
            fnames_uncomp=fnames_uncomp, parent_out=parent_out, outdir=outdir
   endcase
endcase

comp2head=1-keyword_set(comp_header)

if next gt 0 and data_chk(index,/struct) and comp2head then begin 
   ftags=['BITPIX','NAXIS1','NAXIS2']
   ztags='Z'+ftags
   if required_tags(index,ftags) and required_tags(index,ztags) then begin 
      for i=0,n_elements(ftags)-1 do begin 
         index.(tag_index(index(0),ftags(i)))=gt_tagval(index,ztags(i))      
      end   
   endif
endif

if keyword_set(use_index) then begin
  if exist(xllp) then $
     index.crpix1 = (orig_index.crpix1 - xllp)*(orig_index.cdelt1/index.cdelt1)
  if exist(yllp) then $
     index.crpix2 = (orig_index.crpix2 - yllp)*(orig_index.cdelt2/index.cdelt2)
endif

if data_chk(index,/struct) and 1-tag_exist(index(0),'xcen') then begin 
   ; add xcen/ycen
   if required_tags(index(0),'crpix1,cdelt1') then begin 
      xcen=comp_fits_cen(index.crpix1,index.cdelt1,index.naxis1,index.crval1)
      ycen=comp_fits_cen(index.crpix2,index.cdelt2,index.naxis2,index.crval2)
      index=add_tag(index,xcen,'xcen')
      index=add_tag(index,ycen,'ycen')
   endif
endif

if data_chk(index,/struct) then begin ; add some history
   update_history,index,version=version,/caller
   if n_elements(xllp) gt 0 then begin
      if keyword_set(verbose) then box_message,'FOV history'
      if n_elements(orig_index) eq 0 then $
         read_sdo_tom,files,orig_index, $
            only_tags='crval1,crval2,naxis1,naxis2,crpix1,crpix2,cdelt1,cdelt2,crota2,date'
      pinf=replicate(arr2str(strtrim([xllp,yllp,nxp,nyp],2)),n_elements(index))
      update_history,index,/caller,'xll,xyy,nx,ny: ' + pinf,/mode
      update_history,index,/caller,'Orig FILE: ' + files,/mode
      update_history,index,/caller,'Orig DATE: ' + gt_tagval(orig_index,/date,missing=''),/mode
      if not required_tags(orig_index,'crpix1,crpix2,cdelt1,cdelt2') then begin 
         box_message,'Warning: FITS data contains no pointing information!'
         update_history,index,/caller,'NO POINTING INFO!'
      endif else begin 
         fovx=strtrim(gt_tagval(index,/naxis1,miss=4096)*gt_tagval(index,/cdelt1,miss=.5),2)
         fovy=strtrim(gt_tagval(index,/naxis2,miss=4096)*gt_tagval(index,/cdelt2,miss=.5),2)
         update_history,index,/caller,'fovx,fovy: ' + fovx+','+fovy,/mode
         update_history,index,/caller,'Orig CRPIX1,CRPIX2: ' + $
            strtrim(orig_index.crpix1,2)+ ',' + strtrim(orig_index.crpix2,2),/mode
         update_history,index,/caller,'Orig CDELT1,CDELT2: ' + $
            strtrim(orig_index.cdelt1,2)+ ',' + strtrim(orig_index.cdelt2,2),/mode
         update_history,index,/caller,'Orig CROTA2: ' + $
            strtrim(gt_tagval(orig_index,/crota2, $
            missing=gt_tagval(index,/crota1,missing=0.)),2),/mode
      endelse
   endif 
endif

if n_elements(ifiles) eq n_elements(files) then begin
   if required_tags(_extra,'uncomp_delete') then begin 
      box_message,'removing mixed_comp uncompressed'
      ssw_file_delete,files(aac)
   endif
   files=ifiles ; restore input
endif

return
end


