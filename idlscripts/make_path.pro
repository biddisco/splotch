pro make_path
   snr=00.0

   lookat_x=     0
   lookat_y=     0
   lookat_z=     0
   x=-3.0
   z=-2.0

   nframes_old=720
   nframes=400
   dx=6.0/nframes
   dz=4.0/nframes
   
   r0=2.0

  readnew,'ic4247_000',h,'HEAD'
  readnew,'ic4247_000',x1,'POS',part=0
  readnew,'ic4247_000',x2,'POS',part=4
  ii1=lindgen(h.npart(0)/500)*500
  ii2=lindgen(h.npart(4)/100)*100
  plot,x1(0,ii1),x1(1,ii1),psym=3,xr=[-2,2],yr=[-1,3],xst=1,yst=1
  oplot,x2(0,ii2),x2(1,ii2),psym=3,col=255*256LL
;   plot,[0],[0],psym=4,xr=[-2,2],yr=[-1,3],xst=1,yst=1,symsize=2,thick=2

   openw,1,'path.txt',width=1000
   printf,1,'camera_x camera_y camera_z brightness0 brightness1 brightness2 fidx'
   ;printf,1,'camera_x camera_y camera_z fidx'
   for i=0,nframes-1 do begin
      fac=float(i)/(nframes_old-1)
      phi=-fac*4*!Pi

      y=x^2-2.5

      ;rr=(1.2+cos(phi/2))*r0/2
      modulation =  1 - fac * 1.5
      ;x=lookat_x + 0.5*rr*sin(phi)
      ;y=lookat_y + r0*0.85 + rr*cos(phi) -0.75
      ;z=lookat_z + modulation

      bright = (x^2+y^2+z^2)/10
      brightstars = sqrt(bright)
      smooth = bright

      oplot,[x],[y],psym=1,col=255

      bfac=((cos(phi)+1)/2)^(1./4.)
      bbb=1.0
      printf,1,x,y,z,0.5*bbb,0.3*bbb,bbb,snr
      ;printf,1,x,y,z,bright*0.005,bright*0.02,brightstars*0.04,smooth*5,snr
      ;printf,1,x,y,z,snr
      x=x+dx
      z=z+dz
      ;printf,1,x,y,z,snr
   end
   close,1
end
