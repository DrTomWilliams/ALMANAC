PRO almanac, location=location, onset_time=onset_time, $
             start_time=start_time, end_time=end_time, wave=wave, $
             nrt=nrt, test=test, usrdir=usrdir, flist=flist, $
             savdir=savdir, show=show, cadence=cadence, $
             movies=movies, movdir=movdir
  
  Compile_Opt idl2 ; ensure integers are 32-bit by default (defint32) and 
                   ; square brackets must be used for array elements
                   ; (strictarr).
  
  ; location: 1D or 2D array containing the longitude and latitude of the
  ;           CME origin in all wavelengths specified by wave.
  ;
  ; onset_time: String array containing the date and time of CME onset in
  ;             all wavelengths specified by wave.
  ;
  ; start_time: (OPTIONAL) to be used on archival data in format - 
  ;             'yyyy/mm/dd hh:mm'.
  ;
  ; end_time: (OPTIONAL) to be used on archival data in format - 
  ;           'yyyy/mm/dd hh:mm'.
  ;
  ; wave: A single or array of SDO/AIA wavelengths.
  ;
  ; nrt: Enable to use near real time data. DEFAULT: synoptic map data.
  ;
  ; test: To be used with archival data that has been downloaded already.
  ;
  ; usrdir: Only to be used with test on pre-downloaded data.
  ;
  ; flist: list of files rather than a directory. To be used in place of
  ;           usrdir when files are spread across more than one directory.
  ;           Test keyword must be defined.
  ;
  ; savdir: Use to specify a top-level directory for the output files for
  ;         each event detected by the ALMANAC routine.
  ;
  ; show: Set to show frame by frame movies and centre of mass of CME. 
  ;       NB: only works for the first element of wave if wave is an array.
  ;
  ; cadence: Set to desired cadence between images when downloading data
  ;          given in seconds.
  ;
  ; movies: Set to produce MP4s of the CME
  
  dirset=0
  fileset=0
  flag_wave = 0
  flag_synoptic = 1
  flag_gettime = 0
  flag_film = 0
  flag_test = 0
  flag_err1 = 0
  flag_err2 = 0
  flag_err3 = 0
  flag_err4 = 0
  flag_top = 0
  csys = 'stonyhurst'
  flag_show = 0
  frames=0
  flag_no_files = 0
  r_thresh = 2 ; 1.80   ; threshold value for ratio in Boolean mask
  px_thresh = 15e2  ; min px number for a region in Boolean mask
  min_time = 18d0
  smo_time = 30d0
  
  
  
  if keyword_set(usrdir) then begin
    datadir = usrdir
    dirset=1
  endif
  if keyword_set(flist) then begin
    filelist = flist
    fileset = 1
  endif
  if keyword_set(test) then flag_test = 1 else flag_test = 0
  if keyword_set(wave) then begin
    wvlnth = wave
    flag_wave = 1
  endif
  if keyword_set(nrt) then flag_synoptic=0
  if keyword_set(start_time) and keyword_set(end_time) then begin
    time0 = start_time
    time1 = end_time
    flag_gettime = 1
  endif
  if keyword_set(test) then flag_test = 1
  if keyword_set(start_time) and not keyword_set(end_time) then begin
    flag_err3 = 1
    goto, time_error
  endif
  if keyword_set(savdir) then begin
    topdir = savdir
    tmp = strlen(topdir)
    if strmid(topdir,tmp-1,1) ne '/' then topdir+='/'
    flag_top = 1
  endif
  if keyword_set(carrington) then csys = 'carrington'
  if keyword_set(show) then flag_show = 1
  if keyword_set(frames) then min_frames = frames
  if keyword_set(movies) then flag_film = 1
  if keyword_set(cadence) then step = cadence else step = 360
  
  if keyword_set(movdir) then flag_film_dir = 1 else flag_film_dir = 0
  
  if flag_film then restore, file='limb.dat' 
  ; contains the limb pixel coordinates
  
  ; this block of code determines the time and date of the event
  ; to be tracked and generates an observation window to download
  ; the quicklook AIA data (1024x1024 pixels).
  ;===============================================================
  if not flag_test then begin
    if not flag_synoptic then begin
      current = systime(/utc)
      hh = strmid(current,11,2)   ; hour
      mm = strmid(current,14,2)   ; minute
      ss = strmid(current,17,2)   ; second
      yr = strmid(current,20,4)   ; year
      mn = strmid(current,4,3)    ; month
      dy = strmid(current,8,2)    ; day
      check = fix(dy)
    endif else begin
      if typename(end_time) eq 'UNDEFINED' then begin
        flag_err3 = 1
        goto, time_error
      endif
      current = end_time
      hh = strmid(current,11,2)   ; hour
      mm = strmid(current,14,2)   ; minute
      ss = strmid(current,17,2)   ; second
      yr = strmid(current,0,4)    ; year
      mn = strmid(current,5,2)    ; month
      dy = strmid(current,8,2)    ; day
      check = fix(dy)
    endelse
    if check lt 10 then dy = '0'+strmid(current,9,1)
    if not dirset then begin
      datadir = yr+mn+dy
    endif ; dirset = 0
    
    if mn eq 'Jan' then mn = '01'
    if mn eq 'Feb' then mn = '02'
    if mn eq 'Mar' then mn = '03'
    if mn eq 'Apr' then mn = '04'
    if mn eq 'May' then mn = '05'
    if mn eq 'Jun' then mn = '06'
    if mn eq 'Jul' then mn = '07'
    if mn eq 'Aug' then mn = '08'
    if mn eq 'Sep' then mn = '09'
    if mn eq 'Oct' then mn = '10'
    if mn eq 'Nov' then mn = '11'
    if mn eq 'Dec' then mn = '12'
  
    min_chk = fix(round(fix(mm) + fix(ss)/60.))
    mm = 2*round(min_chk/2.)
    if long(mm)/2 ne mm/2. then begin
      if mm+1 gt 59 then mm = mm-1 else mm = mm+1
    endif
    
    year = fix(yr)  &  month = fix(mn)  &  day = fix(dy)
    hour = fix(hh)  &  minute = fix(mm)  &  second = 0
    
    time0 = mytime(YEAR=year, MONTH=month, DAY=day, HOUR=hour, $
                  MINUTE=minute,SECOND=second, OFFSET=8, /UTC)
    time0 = strmid(time0,0,10)+' '+strmid(time0,11,5)
    time1 = mytime(YEAR=year, MONTH=month, DAY=day, HOUR=hour, $
                  MINUTE=minute,SECOND=second, /UTC)
    time1 = strmid(time1,0,10)+' '+strmid(time1,11,5)
    
    ; download NRT data using system time
    if not flag_synoptic then begin
    get_aia_synoptic_data, time0, time1, wvlnth, $
                           increment=step, $
                           topsavedir=datadir, /nrt
  
    ; download data using synoptic data for a specified time range
    endif else begin
      yr = strmid(end_time,0,4)
      mn = strmid(end_time,5,2)
      dy = strmid(end_time,8,2)
      datadir = yr+mn+dy
      get_aia_synoptic_data, time0, time1, wvlnth, $
                             increment=step, $
                             topsavedir=datadir
    endelse
  endif ; not flag_test

  
  
  ; The following segment reads in the SDO quicklook data and then
  ; creates a 3D datacube. The differential rotation of the Sun is
  ; removed, and MGN is then applied.
  ;===============================================================
  if flag_test then begin
    wave = strcompress(wvlnth,/remove_all)
    nwave = wave.length
    if dirset then dir = datadir+'/'+wave+'/'
  endif else begin
    wave = strcompress(wvlnth,/remove_all)
    if not flag_synoptic then $
      dir = datadir+'/'+wave+'/' $
    else $
      dir = datadir+'/'+wave+'/'
    nwave = wave.length
    files = findfile(dir[0])
  endelse
  
  if keyword_set(flist) then files=flist; HUW

  onset_time = strarr(nwave)
  location = dblarr(nwave,2)
  
  for iwave=0,nwave-1 do begin
    if not fileset then begin
      fnames = findfile(dir[iwave])
      read_files = dir[iwave]+fnames
    endif else begin
      read_files = filelist
      a = read_files[0]
      len = strlen(a)
      fnames = strmid(read_files,len-26,26)
    endelse

    tmp = fnames.length-1
    print, '% Processing files from '+strmid(fnames[tmp],3,4)+ $
           '/'+strmid(fnames[tmp],7,2)+'/'+strmid(fnames[tmp],9,2)+ $
           ' in the '+wave[iwave]+' channel. . .'



    ; adjust the minimum number of frames and the smoothing
    ; range based on cadence such that minimum event and 
    ; smoothed time are always constant.
    f1 = read_files[0]
    f2 = read_files[1]
    f3 = read_files[2]
    f1 = fix(strmid(f1, strlen(f1)-12,2))
    f2 = fix(strmid(f2, strlen(f2)-12,2))
    f3 = fix(strmid(f3, strlen(f3)-12,2))
    if f2 le f1 then f2+= 60
    if f3 le f2 then f3+= 60
    td1 = f2-f1
    td2 = f3-f2
    td = min([td1,td2]) ; cadence in mins
    min_frames = round(min_time/td)
    if min_time/td lt 1 then begin
      flag_err4 = 1
      goto, error_frm
    endif
    range = round(smo_time/td)
    

    
    
    ; read in idl files with conditional statement to catch bad files
    ; pointers are used for hdr as these may not always be consistent
    im = !null
    index_good = intarr(read_files.length)  &  index_good[*] = 1
    for i=0, read_files.length-1 do begin

      if read_files.length lt 2*range + 1 then begin
        flag_no_files = 1
        goto, no_files
      endif

      
      ; create a catch for any errors when reading in SDO files.
      catch, error_status
      if error_status then begin
        catch, /cancel
        print, 'error with this file: ', strcompress(i,/remove_all)
        print, 'Error message: ', !ERROR_STATE.MSG
        index_good[i] = 0
        error_status = 0
        continue
      endif ; error_status = true   

      read_sdo_tom, read_files[i], tmp1, tmp2, /noshell, /uncomp_delete
      ; parent_out tells where to put tmp files.
      if im eq !null then begin
        hdr = ptrarr(read_files.length)
        sizearr, tmp2,nx,ny
        im = fltarr(nx,ny,read_files.length)
      endif
      hdr[i] = ptr_new(tmp1)
      im[*,*,i] = tmp2
    endfor ; i 0, read_files.length-1
    
    indgood=where(index_good,cntgood)
    if cntgood eq 0 then begin
      print,'CNTGOOD is 0!'
      continue
    endif
    hdr = hdr[indgood]
    im = im[*,*,indgood]
    
    catch, /cancel

    
    

    im=float(im>0);set negatives to zero
    sizearr,im,nx,ny,nt ;query size of datacube as nx,ny,nt
    for i=0,nt-1 do begin
      ihdr = *hdr[i]
      im[*,*,i]=temporary(im[*,*,i])/ihdr.exptime
    endfor
   

    ; rebin the data if not 1k x 1k
    if nx gt 1024 then begin
       sizearr,im,nx,ny,nt
       im = rebin(im,1024,1024,nt)
       ;prefer to use rebin, since uses local average to shrink images
       zoom=nx/1024.
       sizearr,im,nx,ny,nt
       for it=0,nt-1 do begin 
         h=*hdr[it]
         h.naxis1=nx
         h.naxis2=ny
         h.crpix1=h.crpix1/zoom
         h.crpix2=h.crpix2/zoom
         h.cdelt1=h.cdelt1*zoom
         h.cdelt2=h.cdelt2*zoom
         *hdr[it]=h
       endfor
    endif

    ; eliminate over/under-saturated data and re-define nt
    med = dblarr(nt)
    for i=0,nt-1 do med[i] = median(im[*,*,i])
    medval = median(med)
    good = where(med lt 1.3*medval and med gt 0.7*medval)
    im = im[*,*,good]
    hdr = hdr[good]
    nt = good.length
    if nt lt 2*range+1 then begin
      print, '    Insufficient acceptable files: median check error.'
      goto, opt_exit
    endif

    ; eliminate frame to frame jitter and store the interpolation
    ; shift in x and y in array, xyshift. Re-assign the crpix1 and
    ; crpix2 variables with the shift performed.
    get_shift, im
    
    ndx = (nx-nx*0.75)/2.
    ndy = round((nx-nx*0.65)/2.)
    im = im[ndx:nx-ndx,ndy:ny-ndy,*] ; shrink image array size
    sizearr,im,nx,ny,nt
    ;adjust header
    for it=0,nt-1 do begin
      h=*hdr[it]
      h.naxis1=nx
      h.naxis2=ny
      h.crpix1=h.crpix1-ndx
      h.crpix2=h.crpix2-ndy
      *hdr[it]=h
    endfor
    
    ; standardise image intensity
    val = median(im)
    immax = 2e3
    im /= val
    im *= 150
    im = (im>0)<immax  ; threshold input


    ; For the time differencing, some form of consideration is needed
    ; for the cadence between images in case there are drop-outs that
    ; occur during the observation window.
    ;==============================================================
    nker=2*range+1
    ker=fltarr(nker)+1/float(nker)
    ker=reform(ker,1,1,nker)
    avg=convol(im,ker,/edge_trunc,/nan)
    diff=abs(im-avg)

    ; determine the median absolute changes over time for each pixel
    ; within the SDO NRT data.
    medim = median(abs(diff),dim=3)
  
    ; Determine the ratio of each image to the median value. Median
    ; image is duplicated so that a fast subtraction can be made.
    ratio = diff/temporary(rebin(medim,nx,ny,nt))
    ratio = (ratio>0)<30  ; make sure values are sensible

    ; create mask using region size filter to find locations of
    ; significant variation to the median.
    mask = bytarr(nx,ny,nt);HUW see below couple of lines down
    ; actually create mask now
    for it=0,nt-1 do mask[*,*,it]= region_size_filter(ratio[*,*,it], $
          +                         r_thresh,0.7*px_thresh,/return_mask)
    ;let's smear out the mask to join regions
    kermask=gaussian_function(1.,/norm);create narrow gaussian kernel
    nkermask=kermask.length
    ;apply convolution/smoothing across first x-dimension then y
    masksmo=convol(float(mask),kermask,/edge_trunc,/nan);x-dimension convol
    masksmo=transpose(convol(transpose(masksmo,[1,0,2]),kermask,$
                             /edge_trunc,/nan),[1,0,2]);y-dimension convol
    ; I noticed there was some flickering where an obvious large region
    ; disappears for one frame then reappears in time, 
    ; happened after region_size_filter. This is caused by the hard
    ; threshold px_thresh in region_size_filter.
    ; Fix this by smoothing also in time (but very narrow kernel).
    kert=gaussian_function(0.7,/norm)
    masksmo=convol(masksmo,reform(kert,1,1,n_elements(kert)),$
                   /edge_trunc,/nan)
    ; this threshold not so important, arbitary small value.
    masksmo=masksmo gt median(masksmo)

    ;region size filter. Note applied to smoothed mask, but with result
    ; multiplied by mask to return original shape of regions
    mask2 = mask*0b
    for it=0,nt-1 do mask2[*,*,it] = region_size_filter(masksmo[*,*,it], $
                                                        0.8,px_thresh,$
                                                        /return_mask)

    labelmask = unpad_3d(label_region(pad_3d(mask2,1), /ulong),1)*mask
    ;unfortunately label_region ignores 1 pixel at margins, so use pad
  
    ; labelmask 0 values are where there are no regions of interest.
    ; determine the region that is most likely to be the CME by taking
    ; histogram of the labelled regions e.g. 1-8. Then we need to make
    ; a new array that only takes the values of mask where labelmask is
    ; equal to the most frequent value from the histogram.
    ;
    ; All regions of signficance are outputted to sav environments.
    ;==============================================================
    
    minvoxvol=9e3
    ; filter out region volumes < minvoxvol voxels
    
    if min(labelmask) eq max(labelmask) then begin
      print,'single labelmask region'
      goto, opt_exit
    endif
    hist = histogram(labelmask,location=xhist,min=1,reverse_ind=ri)
    indhistregion=where(hist gt minvoxvol,cntregion)
    
    ; HUW - note I'm using reverse indices here, very useful and quick.
    ; For each populated histogram bin, the reverse indices store the
    ; indices of the input array (labelmask) that contribute to that
    ; histogram bin. So we can loop through the histogram bins of interest
    ; and quickly extract the relevant indices of labelmask, makes life
    ; much easier and quicker. 
    if cntregion eq 0 then begin
      print,'No regions'
      goto, opt_exit
    endif
    
    tmp = nt  ; store nt of array as nt varies in loop.
    ; loop through substantial-sized regions
    for iregion=0,cntregion-1 do begin
      nt = tmp
      ; extract the cube locations of current region  
      indregionnow=get_rev_ind(ri,indhistregion[iregion],nnow)
      
      ;change cube locations to x, y, and t indices
      one2n,indregionnow,[nx,ny,nt],ix,iy,it,/dim
      trange=max(it)-min(it)+1
      
      ; skip if doesn't span min_frames time steps
      if trange lt min_frames then goto, skip_loop
      itime=min(it)
      
      hdrnow=*hdr[itime]
      wcs=fitshead2wcs(hdrnow)
      indtnow=where(it eq itime) ;consider only points at current time step
      ixnow=ix[indtnow]
      iynow=iy[indtnow]
      ixav=total(ixnow*ratio[indregionnow[indtnow]])/$
           total(ratio[indregionnow[indtnow]])
      iyav=total(iynow*ratio[indregionnow[indtnow]])/$
           total(ratio[indregionnow[indtnow]])
      coord=wcs_get_coord(wcs,[ixav,iyav])
      wcs_convert_from_coord, wcs, coord, 'hg', lon, lat
      wcs_convert_from_coord, wcs, coord, 'hg', clon, clat, /carrington
      

      ; create directory tree and save event details as a structure
      yr = strmid(hdrnow.t_obs,0,4)
      mn = strmid(hdrnow.t_obs,5,2)
      dy = strmid(hdrnow.t_obs,8,2)
      dirsav = yr+'/'+mn+'/'
      if flag_top then dirsav = topdir + dirsav
      
      file_mkdir, dirsav
      timestamp = yr+mn+dy+'_'+ $
                  strmid(hdrnow.t_obs,11,2)+ $
                  strmid(hdrnow.t_obs,14,2)

      ; segment region to separate variables for ouput
      cmehdr = hdr[min(it):max(it)]
      cmemask = mask
      notcme = where(labelmask ne indhistregion[iregion]+1)
      cmemask[notcme] = 0
      if flag_film then begin
        tom = max(it)+round(60./td)
        if tom ge nt-1 then tom = nt-1
        cmemask2 = cmemask[*,*,min(it):tom]
        cmeim2 = im[*,*,min(it):tom]
      endif ; flag_film
      cmemask=cmemask[*,*,min(it):max(it)]
      cmeim = im[*,*,min(it):max(it)]
      cmeind = where(cmemask)
      cmerat = ratio[*,*,min(it):max(it)]

      nt=trange
      one2n, cmeind,cmemask,ix,iy,it
      ituniq=it[uniq(it,sort(it))]
      n=ituniq.length
      ixra=lonarr(2,n)
      iyra=ixra
      for ixy=0,n-1 do begin
        ind=where(it eq ituniq[ixy])
        ixra[*,ixy]=minmax(ix[ind])
        iyra[*,ixy]=minmax(iy[ind])
      endfor

      com = dblarr(3,2)
      com[0,*] = [lon,lat]
      com[1,*] = [clon,clat]
      com[2,*] = coord
      
      
      info = {nx:"number of x data points", ny:"number of y datapoints", $
              nt:"number of frames", hdr:"header files for each frame", $
              index:"indices of non-zero elements", $
              values:"intensity values at CME location", $
              ratio:"Difference image(s) values at CME location;"+$
                    " used to define CoM", $
              xtrim:"number of pixels trimmed from left and right of"+ $
                    "each image; e.g. xtrim = 100, then total of 200"+ $
                    " px trimmed.", $
              ytrim:"number of pixels trimmed from top and bottom of"+ $
                    "each image; e.g. ytrim = 100, then total of 200"+ $
                    " px trimmed.", $
              wave:"string indicating wavelength in Angstroms", $
              CoM:"Center of mass for the CME origin given in longitude"+ $
                  ", latitude for Stonyhurst [0,*], Carrington [1,*], "+ $
                  "and Heliographic [3,*] coordinates", $
              totals:"Total intensities of ratio [0], SDO image [1], "+ $
                     "and no. of indices [2]"}

      save, file=dirsav+'almanac_'+timestamp+'_wvlnth_'+wave+'_'+$
            strcompress(iregion,/remove_all)+'.dat', info, nx, ny, nt, $
            cmehdr, cmeim, cmerat, cmemask, cmeind, wave, com, $
            ixra, iyra, /compress

      if flag_film then begin
        tobs = hdrnow.t_obs
        tobs = strmid(tobs,0,10)+'_'+strmid(tobs,11,8)
        tobs = tobs.replace('-','_')
        tobs = tobs.replace(':','_')
        
        ;aia_lct, wavelnth=193, /load
        a = cmeim2/max(cmeim2)*255
        b = cmemask2
        b[where(b gt 0)] = 255
        b[where(b eq 0)] = 1
        b[where(b eq 255)] = 0
        b *= 255
        c = byte([a,b])
        if flag_film_dir then mov_dir = movdir else mov_dir = 'movies/'
        moviemaker, images=c, dir=mov_dir, $
                    /use_mgn, movie_name='almanac_'+$
                    timestamp+'_wvlnth_'+wave+'_'+$
                    strcompress(iregion,/remove_all)+ $
                    '.mp4', xsize=2*nx, ysize=ny, $
                    limb=limbpix, /double_image, tstamp=nhdr, dtype='Z'
        loadct, 0,/silent
      endif ; flag_film
      skip_loop:
    endfor
    no_files: if flag_no_files then begin
      print, 'Files are missing. . .'
    endif ; flag_no_files
  endfor ; iwave
  
  ; exit statements when incorrect call procedure is detected.
  time_error: if flag_err3 then print, $
                'End-time not supplied. . .'
  error_onset: if flag_err2 then print, $
                'Onset time variable not supplied. . .'
  error_loc: if flag_err1 then print, $
                'Location variable is undefined. . .'
  error_frm: if flag_err4 then print, $
                'Insufficient cadence, try smaller values. . .'
  opt_exit: ; conditional exit when no regions.
END ; almanac
