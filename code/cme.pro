  
  ; create a new program that calls almanac to download and process a
  ; bunch of cmes in 171 and 193 and stores the data to a sav environment.
  
  ; Use WCS to get coordinates and set off-limb ratio values to 0. Might
  ; help with noise like in examples 8 and 20.

PRO CME
  aia_lct, wave=193, /load
  nt = 20
  shc = dblarr(nt,2)
  time = strarr(nt)
  cmelist = !NULL  &  randlist = !NULL
  
  ; used for 6 minute cadence
  date1e = '2010/08/14 10:30'  &  date1s = '2010/08/14 02:30'
  date2e = '2011/02/15 03:00'  &  date2s = '2011/02/14 19:00'
  date3e = '2011/09/06 03:00'  &  date3s = '2011/09/05 19:00'
  date4e = '2011/11/26 08:00'  &  date4s = '2011/11/26 00:00'
  date5e = '2012/04/07 21:45'  &  date5s = '2012/04/07 13:45'
  date6e = '2012/11/23 14:15'  &  date6s = '2012/11/23 06:15'
  date7e = '2012/11/27 03:00'  &  date7s = '2012/11/26 19:00'
  date8e = '2013/03/15 07:50'  &  date8s = '2013/03/14 23:50'
  date9e = '2013/04/11 07:55'  &  date9s = '2013/04/10 23:55'
  date10e = '2013/05/15 02:35'  &  date10s = '2013/05/14 18:35'
  date11e = '2013/07/09 15:45'  &  date11s = '2013/07/09 07:45'
  date12e = '2013/08/20 08:45'  &  date12s = '2013/08/20 00:45'
  date13e = '2013/09/29 22:45'  &  date13s = '2013/09/29 14:45'
  date14e = '2014/02/18 02:00'  &  date14s = '2014/02/17 18:00'
  date15e = '2014/03/23 04:10'  &  date15s = '2014/03/22 20:10'
  date16e = '2014/04/01 17:20'  &  date16s = '2014/04/01 09:20'
  date17e = '2014/04/29 23:59'  &  date17s = '2014/04/29 16:00'
  date18e = '2014/08/15 18:20'  &  date18s = '2014/08/15 10:20'
  date19e = '2014/08/22 11:45'  &  date19s = '2014/08/22 03:45'
  date20e = '2014/12/21 12:45'  &  date20s = '2014/12/21 04:45'
  
  ; used for 10 minute cadence
  date1e = '2010/08/14 10:30'  &  date1s = '2010/08/14 02:30'
  date2e = '2011/02/15 03:00'  &  date2s = '2011/02/14 19:00'
  date3e = '2011/09/06 03:00'  &  date3s = '2011/09/05 19:00'
  date4e = '2011/11/26 08:00'  &  date4s = '2011/11/26 00:00'
  date5e = '2012/04/07 21:40'  &  date5s = '2012/04/07 13:40'
  date6e = '2012/11/23 14:10'  &  date6s = '2012/11/23 06:10'
  date7e = '2012/11/27 03:00'  &  date7s = '2012/11/26 19:00'
  date8e = '2013/03/15 07:50'  &  date8s = '2013/03/14 23:50'
  date9e = '2013/04/11 07:50'  &  date9s = '2013/04/10 23:50'
  date10e = '2013/05/15 02:30'  &  date10s = '2013/05/14 18:30'
  date11e = '2013/07/09 15:40'  &  date11s = '2013/07/09 07:40'
  date12e = '2013/08/20 08:40'  &  date12s = '2013/08/20 00:40'
  date13e = '2013/09/29 22:40'  &  date13s = '2013/09/29 14:40'
  date14e = '2014/02/18 02:00'  &  date14s = '2014/02/17 18:00'
  date15e = '2014/03/23 04:10'  &  date15s = '2014/03/22 20:10'
  date16e = '2014/04/01 17:20'  &  date16s = '2014/04/01 09:20'
  date17e = '2014/04/29 23:50'  &  date17s = '2014/04/29 16:00'
  date18e = '2014/08/15 18:20'  &  date18s = '2014/08/15 10:20'
  date19e = '2014/08/22 11:40'  &  date19s = '2014/08/22 03:40'
  date20e = '2014/12/21 12:40'  &  date20s = '2014/12/21 04:40'

  tic
  for it =0,nt-1 do begin
  ;for it =8,8 do begin ; for de-bugging
  ;for it =12,12 do begin ; use for figures
    if it eq 0 then begin
      d1 = date1s  &  d2 = date1e
    endif  
    if it eq 1 then begin
      d1 = date2s  &  d2 = date2e
    endif  
    if it eq 2 then begin
      d1 = date3s  &  d2 = date3e
    endif  
    if it eq 3 then begin
      d1 = date4s  &  d2 = date4e
    endif  
    if it eq 4 then begin
      d1 = date5s  &  d2 = date5e
    endif  
    if it eq 5 then begin
      d1 = date6s  &  d2 = date6e
    endif  
    if it eq 6 then begin
      d1 = date7s  &  d2 = date7e
    endif  
    if it eq 7 then begin
      d1 = date8s  &  d2 = date8e
    endif
    if it eq 8 then begin
      d1 = date9s  &  d2 = date9e
    endif  
    if it eq 9 then begin
      d1 = date10s  &  d2 = date10e
    endif  
    if it eq 10 then begin
      d1 = date11s  &  d2 = date11e
    endif  
    if it eq 11 then begin
      d1 = date12s  &  d2 = date12e
    endif  
    if it eq 12 then begin
      d1 = date13s  &  d2 = date13e
    endif
    if it eq 13 then begin
      d1 = date14s  &  d2 = date14e
    endif  
    if it eq 14 then begin
      d1 = date15s  &  d2 = date15e
    endif  
    if it eq 15 then begin
      d1 = date16s  &  d2 = date16e
    endif  
    if it eq 16 then begin
      d1 = date17s  &  d2 = date17e
    endif
    if it eq 17 then begin
      d1 = date18s  &  d2 = date18e
    endif  
    if it eq 18 then begin
      d1 = date19s  &  d2 = date19e
    endif  
    if it eq 19 then begin
      d1 = date20s  &  d2 = date20e
    endif  
    
    ; call ALMANAC and process the example CMEs.
    dir = strmid(d2,0,4)+strmid(d2,5,2)+strmid(d2,8,2)
    almanac, location=tmp1, onset_time=tmp2, wave=[193], $
             usrdir=dir, /test, /movies;, $
             ;start_time=d1, end_time=d2, cadence=600
    
    ; query the cdaw catalogue for the CME info
    c = ssw_getcme_list(d1,d2,/cdaw)
    if typename(cmelist) eq 'UNDEFINED' then cmelist = ptr_new(c) else $
                   cmelist = [cmelist, ptr_new(c)]
  endfor ; it
  toc
  save, file='cmelist.sav', cmelist
  print,''
  print,''
  print,''
  
  ; used for 6 minute cadence
  date1e = '2017/10/25 14:45'  &  date1s = '2017/10/25 06:45'
  date2e = '2018/01/31 23:00'  &  date2s = '2018/01/31 15:00'
  date3e = '2018/06/27 01:00'  &  date3s = '2018/06/26 17:00'
  date4e = '2018/08/28 08:00'  &  date4s = '2018/08/28 00:00'
  date5e = '2019/04/11 21:45'  &  date5s = '2019/04/11 13:45'
  date6e = '2019/05/25 14:55'  &  date6s = '2019/05/25 06:55'
  date7e = '2019/12/26 23:30'  &  date7s = '2019/12/26 15:30'
  date8e = '2020/03/13 04:50'  &  date8s = '2020/03/12 20:50'
  date9e = '2020/04/19 12:05'  &  date9s = '2020/04/19 04:05'
  date10e = '2020/04/23 05:25'  &  date10s = '2020/04/22 21:25'
  date11e = '2020/06/07 03:00'  &  date11s = '2020/06/06 19:00'
  date12e = '2020/07/03 18:05'  &  date12s = '2020/07/03 12:05'
  date13e = '2020/09/08 13:45'  &  date13s = '2020/09/08 05:45'
  date14e = '2021/01/10 02:00'  &  date14s = '2021/01/09 18:05'
  date15e = '2021/03/02 06:10'  &  date15s = '2021/03/01 22:10'
  date16e = '2021/05/13 17:20'  &  date16s = '2021/05/13 09:20'
  date17e = '2021/06/28 23:35'  &  date17s = '2021/06/28 15:35'
  date18e = '2021/07/07 10:15'  &  date18s = '2021/07/07 02:15'
  date19e = '2021/10/31 19:30'  &  date19s = '2021/10/31 11:30'
  date20e = '2021/11/17 08:00'  &  date20s = '2021/11/17 00:00'

  ; used for 10 minute cadence
  date1e = '2017/10/25 14:40'  &  date1s = '2017/10/25 06:40'
  date2e = '2018/01/31 23:00'  &  date2s = '2018/01/31 15:00'
  date3e = '2018/06/27 01:00'  &  date3s = '2018/06/26 17:00'
  date4e = '2018/08/28 08:00'  &  date4s = '2018/08/28 00:00'
  date5e = '2019/04/11 21:40'  &  date5s = '2019/04/11 13:40'
  date6e = '2019/05/25 14:50'  &  date6s = '2019/05/25 06:50'
  date7e = '2019/12/26 23:30'  &  date7s = '2019/12/26 15:30'
  date8e = '2020/03/13 04:50'  &  date8s = '2020/03/12 20:50'
  date9e = '2020/04/19 12:00'  &  date9s = '2020/04/19 04:00'
  date10e = '2020/04/23 05:20'  &  date10s = '2020/04/22 21:20'
  date11e = '2020/06/07 03:00'  &  date11s = '2020/06/06 19:00'
  date12e = '2020/07/03 18:00'  &  date12s = '2020/07/03 12:00'
  date13e = '2020/09/08 13:40'  &  date13s = '2020/09/08 05:40'
  date14e = '2021/01/10 02:00'  &  date14s = '2021/01/09 18:00'
  date15e = '2021/03/02 06:10'  &  date15s = '2021/03/01 22:10'
  date16e = '2021/05/13 17:20'  &  date16s = '2021/05/13 09:20'
  date17e = '2021/06/28 23:30'  &  date17s = '2021/06/28 15:30'
  date18e = '2021/07/07 10:10'  &  date18s = '2021/07/07 02:10'
  date19e = '2021/10/31 19:30'  &  date19s = '2021/10/31 11:30'
  date20e = '2021/11/17 08:00'  &  date20s = '2021/11/17 00:00'

  tic
  for it =0,nt-1 do begin
  ;for it =1,1 do begin
    if it eq 0 then begin
      d1 = date1s  &  d2 = date1e
    endif  
    if it eq 1 then begin
      d1 = date2s  &  d2 = date2e
    endif  
    if it eq 2 then begin
      d1 = date3s  &  d2 = date3e
    endif  
    if it eq 3 then begin
      d1 = date4s  &  d2 = date4e
    endif  
    if it eq 4 then begin
      d1 = date5s  &  d2 = date5e
    endif  
    if it eq 5 then begin
      d1 = date6s  &  d2 = date6e
    endif  
    if it eq 6 then begin
      d1 = date7s  &  d2 = date7e
    endif  
    if it eq 7 then begin
      d1 = date8s  &  d2 = date8e
    endif
    if it eq 8 then begin
      d1 = date9s  &  d2 = date9e
    endif  
    if it eq 9 then begin
      d1 = date10s  &  d2 = date10e
    endif  
    if it eq 10 then begin
      d1 = date11s  &  d2 = date11e
    endif  
    if it eq 11 then begin
      d1 = date12s  &  d2 = date12e
    endif  
    if it eq 12 then begin
      d1 = date13s  &  d2 = date13e
    endif
    if it eq 13 then begin
      d1 = date14s  &  d2 = date14e
    endif  
    if it eq 14 then begin
      d1 = date15s  &  d2 = date15e
    endif  
    if it eq 15 then begin
      d1 = date16s  &  d2 = date16e
    endif  
    if it eq 16 then begin
      d1 = date17s  &  d2 = date17e
    endif
    if it eq 17 then begin
      d1 = date18s  &  d2 = date18e
    endif  
    if it eq 18 then begin
      d1 = date19s  &  d2 = date19e
    endif  
    if it eq 19 then begin
      d1 = date20s  &  d2 = date20e
    endif  
    
    ; call ALMANAC and process the example CMEs.
    dir = strmid(d2,0,4)+strmid(d2,5,2)+strmid(d2,8,2)
    almanac, location=tmp1, onset_time=tmp2, wave=[193], $
             usrdir=dir, /test, /movies
             ;start_time=d1, end_time=d2, cadence=600
    ; query the cdaw catalogue for the CME info
    c = ssw_getcme_list(d1,d2,/cdaw)
    if typename(randlist) eq 'UNDEFINED'  then randlist = ptr_new(c) else $
                   randlist = [randlist, ptr_new(c)]
  endfor ; it
  toc
  save, file='randlist.sav', randlist
END ; CME
