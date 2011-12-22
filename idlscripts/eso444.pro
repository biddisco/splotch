pro setpath
   COMMON filepath, path, nx, ny, nxi, nyi, galid
   path='../ObsData/'
   nxi=180
   nyi=180
   nx=180
   ny=180
   galid='ESO444-G084'
   
   PRINT, path, nx, ny, nxi, nyi, galid
   RETURN
END

pro hi_colors
   COMMON filepath, path, nx, ny, nxi, nyi, galid

   mt=[[0.0,0.5,1.],[0.,1.,1.],[0.,1.,0.5]]
   mt=[[0.1,0.25,1.],[0.25,1.,0.1],[1.,0.1,0.25]]
   namein=['med']

   window,xsize=nxi,ysize=nyi

   dd=fltarr(3,nx,ny)
   ddi=fltarr(3,nxi,nyi)
   dd(*,*,*)=0
   ddi(*,*,*)=0

   files=galid+'.HImom0-'+namein+'.regrid.fits'
   ofile=galid+'_hi_color.dat'

   min_use=[0.0,0.0]
   max_use=[2.0,0.7]

   for i=0,N_ELEMENTS(namein)-1 do begin
      aa=readfits(path+files(i),/silent)

      bb=aa
      jj=where(bb GT -1e10)
      print,'Min/Max Data (',i,'):',minmax(bb(jj))

; log scale
;     bb(*,*)=min(alog10(aa(jj)))
;     bb(jj)=alog10(aa(jj))
; sqrt scale
;     bb(*,*)=0
;     bb(jj)=sqrt(aa(jj))
; linear scale
     bb(*,*)=0
     bb(jj)=aa(jj)

      cc=congrid(bb,nxi,nyi)

      for k=0,2 do begin
         ddi(k,*,*) = ddi(k,*,*) + mt(k,i)*cc(*,*)*255
         dd(k,*,*) = dd(k,*,*) + mt(k,i)*bb(*,*)*255
      end
   end

   ddi(*,*,*) = ddi(*,*,*)/max(ddi)*255
   dd(*,*,*) = dd(*,*,*)/max(dd)*255

   eei=bytarr(3,nxi,nyi)
   eei(*,*,*)=ddi(*,*,*)

   ee=bytarr(3,nx,ny)
   ee(*,*,*)=dd(*,*,*)

   tv,eei,/true

   openw,1,ofile
   writeu,1,ee
   close,1


   stop
end

pro hi_thick
   COMMON filepath, path, nx, ny, nxi, nyi, galid

   window,xsize=nxi,ysize=nyi
   namein=['med']

   dd=fltarr(nx,ny)
   ddi=fltarr(nxi,nyi)
   dd = 0
   ddi = 0

   files=galid+'.HImom2-'+namein+'.regrid.fits'
   ofile=galid+'_hi_height.dat'
   min_use=[0.,1.0]
   max_use=[14.0,15.0]

   for i=0,N_ELEMENTS(namein)-1 do begin
      print, 'iteration', i
      aa=readfits(path+files(i),/silent)
      bb=aa
      jj=where(bb GT -1e10)
      print,'Min/Max Data (',i,'):',minmax(bb(jj))

      bb(*,*)=0
      bb(jj)=aa(jj)

      cc=congrid(bb,nxi,nyi)

      ddi = ddi + cc
      dd = dd + bb

   end

   eei=bytarr(nxi,nyi)
   eei(*,*)=sqrt(ddi(*,*)/max(ddi))*255

   ee=bytarr(nx,ny)
   ee(*,*)=sqrt(dd(*,*)/max(dd))*255

   tv,eei

   openw,1,ofile
   writeu,1,ee
   close,1


   stop
end

pro optical_colors
   COMMON filepath, path, nx, ny, nxi, nyi, galid

   mt=[[1.,0.25,0.25],[0.25,1.,1.],[0.,0.,1.]]
   mt=[[1.,0.25,0.25],[1.0,1.,0.5],[0.,0.,1.]]

   mt=[[1.0,0.60,0.25],[0.50,0.30,0.20],[0.,0.,0.]]

;   device,true_color=24,decomposed=0,retain=2
   window,xsize=nxi,ysize=nyi

   dd=fltarr(3,nx,ny)
   ddc=fltarr(3,nx,ny)
   ddi=fltarr(3,nxi,nyi)

   files=galid+'.dss2'+['red','blue']+'.mask.regrid.fits'

   min_use=[21000.,16000.]
   max_use=[62000.,82000.]

   for i=0,N_ELEMENTS(files)-1 do begin
      aa=readfits(path+files(i),/silent)

      bb=aa
      bbc=aa
      jj=where(bb GT -1e10)
      print,'Min/Max Data (',i,'):',minmax(bb(jj))

      ii=where(bb LT min_use(i),ndo)

      if(ndo GT 0) then bb(ii)=min_use(i)
      ii=where(bb GT max_use(i),ndo)
      if(ndo GT 0) then bb(ii)=max_use(i)
      print,'Min/Max used (',i,'):',minmax(bb(jj))

      bb(jj)=bb(jj)-min(bb(jj))
      bb(jj)=bb(jj)/(max(bb(jj))-min(bb(jj)))

      aa=bb

     bb(*,*)=0
     bbc(*,*)=0
; log scalin
;     bb(*,*)=min(alog10(aa(jj)))
;     bb(jj)=alog10(aa(jj))
; sqrt scaling
     bb(jj)=sqrt(aa(jj))
     bbc(jj)=sqrt(aa(jj))
; power scaling
;      bb(jj)=aa(jj)^2
; linear scaling
;     bb(jj)=aa(jj)

      cc=congrid(bb,nxi,nyi)
      print,'Map ',i,' color index =',mt(0,i),mt(1,i),mt(2,i)
      for k=0,2 do begin
         ddi(k,*,*) = ddi(k,*,*) + mt(k,i)*cc(*,*)*255
         dd(k,*,*) = dd(k,*,*) + mt(k,i)*bb(*,*)*255
         ddc(k,*,*) = ddc(k,*,*) + mt(k,i)*bbc(*,*)*255
      end
          
   end

   ddi = ddi(*,*,*)/max(ddi)*255
   dd = dd(*,*,*)/max(dd)*255
   ddc = ddc(*,*,*)/max(ddc)*255

   eei=bytarr(3,nxi,nyi)
   x0=(nx-nyi)/2
   y0=(ny-nyi)/2

   x00=nx/2
   y00=ny/2

   eei(*,*,*)=ddi(*,*,*)
;   eei(*,*,*)=dd(*,x0:x0+nxi-1,y0:y0+nyi-1)

   ee=bytarr(3,nx,ny)
   ee(*,*,*)=ddc(*,*,*)
   print, 'Min/Max Color (R):',minmax(ee(0,*,*))
   print, 'Min/Max Color (G):',minmax(ee(1,*,*))
   print, 'Min/Max Color (B):',minmax(ee(2,*,*))

;   tv,eei,/true,xsize=nxi,ysize=nyi
;   save_screen,path+galid+'_dss2_color'

   openw,1,galid+'_optical_color.dat'
   writeu,1,ee
   close,1


   ff=fltarr(nx,ny)
   ff(*,*)=dd(0,*,*)+dd(1,*,*)+dd(2,*,*)

   ;gg=filter_image(ff,fwhm=60)
   ;jj=where(gg LT 30)

   ;ff(*,*)=dd(0,*,*)+dd(1,*,*)+dd(2,*,*)
   ;ff(jj) = 0
   ;jj=where(ff LT 1)
   ;ff(jj)=0


   iff=congrid(ff,nxi,nyi)
   iff(*,*)=ff(x0:x0+nxi-1,y0:y0+nyi-1)
   tv,iff/max(iff)*255
;   save_screen,path+galid+'_optical_mask'

   ff=ff-min(ff)
   ff=ff/max(ff)*255
   mm=bytarr(nx,ny)
   mm(*,*)=ff(*,*)

   openw,1,galid+'_optical_mask.dat'
   writeu,1,mm
   close,1

  stop
end

pro uv_colors,log=log,wurzel=wurzel,quad=quad
   COMMON filepath, path, nx, ny, nxi, nyi, galid

   mt=[[0.25,0.25,1.],[1.0,0.25,0.5],[0.,0.,1.]]
   mt=[[0.5,0.25,1.],[0.5,1.0,1.0],[0.,0.,1.]]

   mt=[[0.2,0.2,1.],[0.3,0.3,1.0],[0.,0.,1.]]

   window,xsize=nxi,ysize=nyi

   dd=fltarr(3,nx,ny)
   ddi=fltarr(3,nxi,nyi)

   files=galid+['.FUV','.NUV']+'.convol5.mask.regrid.fits'

   min_use=[-2.,-2.]
   max_use=[0.,1.]

   for i=0,N_ELEMENTS(files)-1 do begin
      aa=readfits(path+files(i),/silent)

      bb=aa
      jj=where(bb GT 0)

; log scaling
      bb(*,*)=0
      bb(jj)=alog10(aa(jj))
; lin scaling
;      bb(*,*)=0
;      bb(jj)=aa(jj)

      print,'Min/Max Data (',i,'):',minmax(bb(jj))

      ii=where(bb LT min_use(i),ndo)
      if(ndo GT 0) then bb(ii)=min_use(i)
      ii=where(bb GT max_use(i),ndo)
      if(ndo GT 0) then bb(ii)=max_use(i)
      print,'Min/Max used (',i,'):',minmax(bb)

      bb(jj)=bb(jj)-min(bb(jj))
      bb(jj)=bb(jj)/(max(bb(jj))-min(bb(jj)))

      cc=congrid(bb,nxi,nyi)
      print,'Map ',i,' color index =',mt(0,i),mt(1,i),mt(2,i)
      for k=0,2 do begin
         ddi(k,*,*) = ddi(k,*,*) + mt(k,i)*cc(*,*)*255
         dd(k,*,*) = dd(k,*,*) + mt(k,i)*bb(*,*)*255
      end
   end

    ddi = ddi(*,*,*)/max(ddi)*255
    dd = dd(*,*,*)/max(dd)*255

   eei=bytarr(3,nxi,nyi)
   x0=(nx-nyi)/2
   y0=(ny-nyi)/2
   eei(*,*,*)=ddi(*,*,*)
;   eei(*,*,*)=dd(*,x0:x0+nxi-1,y0:y0+nyi-1)

   ee=bytarr(3,nx,ny)
   ee(*,*,*)=dd(*,*,*)

   tv,eei,/true
;   save_screen,path+galid+'_uv_color'

   openw,1,galid+'_uv_color.dat'
   writeu,1,ee
   close,1
   
   ff=fltarr(nx,ny)
   ff(*,*)=dd(0,*,*)+dd(1,*,*)+dd(2,*,*)

   iff=congrid(ff,nxi,nyi)
   iff(*,*)=ff(x0:x0+nxi-1,y0:y0+nyi-1)
   tv,iff/max(iff)*255

   ff=ff-min(ff)
   ff=ff/max(ff)*255
   mm=bytarr(nx,ny)
   mm(*,*)=ff(*,*)

   openw,1,galid+'_uv_mask.dat'
   writeu,1,mm
   close,1

  stop
end

