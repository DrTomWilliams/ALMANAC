function int2str_huw2,num,width=width
  return,string(num,format=n_elements(width) eq 0?'(I0)':$
    '(I' + int2str(width) + '.' + int2str(width) +')')
end

function int2str_huw,num,width
  return,string(num,format=n_params() lt 2?'(I0)':$
    '(I' + int2str_huw2(width) + '.' + int2str_huw2(width) +')')
end


pro wrap,test=test

;dates=['2011/01/01 00:00','2011/01/01 01:00','2011/01/01 02:00','2011/01/01 03:00','2011/01/01 04:00']
;dates=[,]; 17:22','2017/08/21 17:30','2017/08/21 17:50']
;sd='2013/11/15 00:00'
;ed='2013/11/22 00:00'
;sd='2011/02/12 00:00'
;ed='2011/02/18 12:00'
sd='2011/03/04 00:00'
ed='2011/03/11 00:00'
wl=[094,131,171,193,211,335,304]
get_aia_synoptic_data,sd,ed,wl,increment=10*60.,test=test



end

;OPTIONAL KEYWORDS
;attempted_files: returns 1d array of all the filenames that the program has attempted to download, regardless of download success.
;status=status, returns 1d binary array of zeros and ones (False and True) indicating download success.
;Status = 1:either file was downloaded successfully, or file already exists locally.
;Status = 0: failure in download, and file does not exist locally.

pro get_aia_synoptic_data,tstartusr,tendusr,wlusr,test=test, $
      increment=increment,topsavedir=topsavedir,verbose=verbose, $
      nrt=nrt,status=status,attempted_files=attempted_files


if n_params() eq 0 then begin
  tstartusr='2018/10/01 00:00'
  tendusr='2020/05/01 00:00';+systime2yymmdd()
  wlusr=[094,131,193,335]
  increment=3600;120.
endif

tstart=tstartusr
tend=tendusr
wl=n_elements(wlusr) eq 0?[094,131,171,193,211,335,304]:wlusr

if strlen(tstart) lt 12 then tstart=tstart+' 00:00'
if strlen(tend) lt 12 then tend=tend+' 23:59'

;increment in seconds
increment=keyword_set(increment)?increment:60.

reqtime=anytim2cal(timegrid(tstart,tend,seconds=increment),form=11)
reqdate11=anytim2cal(reqtime,/da,form=11)
reqdate08=anytim2cal(reqtime,/da,form=8)
reqtime08=keyword_set(nrt)?strmid(anytim2cal(reqtime,/tim,form=8),0,6):strmid(anytim2cal(reqtime,/tim,form=8),0,4)
reqhour='H'+strmid(reqtime08,0,2)+'00'

savedir=keyword_set(topsavedir)?topsavedir:getenv('AIA_SYNOPTIC')
;if keyword_set(nrt) then savedir=savedir+'/nrt'

wlstr3=int2str_huw(wl,3)
wlstr4=int2str_huw(wl,4)

remotetop='http://jsoc2.stanford.edu/data/aia/synoptic'
if keyword_set(nrt) then remotetop+='/nrt'


oUrl = OBJ_NEW('IDLnetUrl')   ; create a new IDLnetURL object
oUrl->SetProperty, VERBOSE = keyword_set(verbose) ; Set verbose to 1 to see more info on the transacton

status=0b
attempted_files=''

for idate=0,n_elements(reqtime)-1 do begin
  for iwl=0,n_elements(wl)-1 do begin
    
    localdir=savedir+'/'+wlstr3[iwl]
    filenamenow='AIA'+reqdate08[idate]+'_'+reqtime08[idate]+'_'+wlstr4[iwl]+'.fits'
    
    attempted_files=[attempted_files,filenamenow]
    
    if file_exist(localdir+'/'+filenamenow) then begin
      print,'File exists ',localdir+'/'+filenamenow
      status=[status,1b]
      continue
    endif 
    
    file_mkdir,localdir
    
    url=remotetop+'/'+reqdate11[idate]+'/'+reqhour[idate]+'/'+filenamenow

    if keyword_set(test) then begin
      print,localdir+'/'+filenamenow
      print,url
      continue
    endif
    
    print,'Downloading ',localdir+'/'+filenamenow
    
    parts=parse_url(url)
    oUrl->SetProperty, url_scheme = parts.scheme ; Set the transfer protocol
    oUrl->SetProperty, URL_HOST = parts.host;
    oUrl->SetProperty, URL_PATH = parts.path
    sock_copy,url,out_dir=localdir,status=statusnow
    wait,0.2;courtesy to JSOC!
    status=[status,statusnow ne 0]
    
  endfor
endfor
    
obj_destroy,oURL
if n_elements(status) gt 1 then begin
  status=status[1:*]
  attempted_files=attempted_files[1:*]
endif

end

