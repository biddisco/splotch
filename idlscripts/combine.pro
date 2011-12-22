pro combine

name_in=['../m83/Model/m83_000','../ic4247/Model/ic4247_000','../eso444/Model/eso444_000']
;name_in=['/DATA/SERPENS_2/dol038/m83/Model/snap_051','../ic4247/Model/ic4247_000','../eso444/Model/eso444_000']
center=[[0,0,0],[0.5,0.5,0.5],[-0.5,-0.5,-0.5]]
scale=[1.0,0.10,0.15]

name_out='combined_000'

readnew,name_in(0),h,'HEAD'
FOR l=1,N_ELEMENTS(name_in)-1 DO BEGIN
  readnew,name_in(l),hthis,'HEAD'
  h.npart(*) += hthis.npart(*)
  h.parttotal(*) += hthis.parttotal(*)
END

write_head,name_out,h

nall=h.npart(0)+h.npart(1)+h.npart(2)+h.npart(3)+h.npart(4)+h.npart(5)

x=fltarr(3,nall)
icount=0L

FOR i=0,5 DO BEGIN
   FOR l=0,N_ELEMENTS(name_in)-1 DO BEGIN
      readnew,name_in(l),hthis,'HEAD'
      nthis=hthis.npart(i)
      IF(nthis GT 0) THEN BEGIN
         readnew,name_in(l),xthis,'POS',part=i
         FOR k=0,2 DO BEGIN
	     x(k,icount:icount+nthis-1)=xthis(k,*)*scale(l)+center(k,l)
         END
         icount=icount+nthis
      END
   END
END

add_block,name_out,x,'POS'

BASE_BLOCKS=['HSM','COL','INT']
BASE_VECTOR=[0,1,0]

FOR i=0,5 DO BEGIN
   IF h.npart(i) GT 0 THEN BEGIN

      FOR j=0,N_ELEMENTS(BASE_BLOCKS)-1 DO BEGIN
         IF BASE_VECTOR(j) EQ 0 THEN x=fltarr(h.npart(i)) $
                                ELSE x=fltarr(3,h.npart(i))   
         icount=0L

         label=BASE_BLOCKS(j)+string(i,form='(i1)')
         IF BASE_VECTOR(j) EQ 0 THEN BEGIN
            type="FLOAT"
            ndim=1
         END ELSE BEGIN
            type="FLOATN"
            ndim=3
         END
         is_present=[0,0,0,0,0,0]
         is_present[i]=1

         FOR l=0,N_ELEMENTS(name_in)-1 DO BEGIN
            readnew,name_in(l),hthis,'HEAD'
            nthis=hthis.npart(i)
            IF(nthis GT 0) THEN BEGIN
               readnew,name_in(l),xthis,label,type=type,ndim=ndim,is_present=is_present
	       IF BASE_VECTOR(j) EQ 0 THEN x(icount:icount+nthis-1) = xthis(*) $
                                      ELSE x(*,icount:icount+nthis-1) = xthis(*,*)
               IF j EQ 0 THEN BEGIN
	          x(icount:icount+nthis-1) = x(icount:icount+nthis-1)*scale(l)
               END
               IF j EQ 2 THEN BEGIN
	          x(icount:icount+nthis-1) = x(icount:icount+nthis-1)*scale(l)*scale(l)
               END
               icount += nthis
            END
         END

         IF icount NE h.npart(i) THEN BEGIN
             print,'Something is wrong, i have to stop here ...'
             stop
         END

         add_block,name_out,x,label

      END
   END
END

end
