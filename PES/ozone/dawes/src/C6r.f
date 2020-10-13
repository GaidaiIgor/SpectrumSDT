      subroutine c6r(mla,mpla,Ma,ipp,Mb,c6lb)
      use splin
* a partir du fichier tabCn.res
      implicit real*8(a-h,o-z)
      parameter (Nip=1,lmax=2)
      dimension C6lb(0:lmax)
      dimension xi(nr),yi(nr),b(nr),c(nr),d(nr)
      common/droo/roo
      common/dcc/ic
      data l1,l1p/1,1/
      data l2,l2p/1,1/
      data z0,un/0.d0,1.d0/
      n=l1+l1p+l2+l2p+2
       xl1=dfloat(l1)
       xl1p=dfloat(l1p)
       xl2=dfloat(l2)
       xl2p=dfloat(l2p)
c      write(6,101)l1,l1p,l2,l2p
 101   format('# l1 l1p l2 l2p =',4i3)
       do lb=0,lmax
       c6lb(lb)=0.d0
       enddo
       do lb=iabs(Mb),lmax
       c6t=0.d0
       do la=iabs(Ma),lmax
       ml=min0(la,lb)
       do lambda=iabs(la-lb),la+lb
       enddo
       enddo
       ic=ic+1
       call initspline(ic,xi,yi,b,c,d)
       c6t=ispline(roo,xi,yi,b,c,d,nr)
!      read(10,*)imla,impla,jpp,iMa,iMb,iLb,c6t
       imla=ll(ic,1)
       impla=ll(ic,2)
       jpp=ll(ic,3)
       iMa=ll(ic,4)
       iMb=ll(ic,5)
       iLb=ll(ic,6)
c       write(11,106),imla,impla,jpp,iMa,iMb,iLb,c6t
       if(iabs(imla-mla).ne.0)then
       print *,'# erreur de lecture de Mla',imla,mla
       stop
       endif
       if(iabs(impla-mpla).ne.0)then
       print *,'# erreur de lecture de Mpla',impla,mpla
       stop
       endif
       if(iabs(ima-ma).ne.0)then
       print *,'# erreur de lecture de Ma',ima,ma
       stop
       endif
       if(iabs(imb-mb).ne.0)then
       print *,'# erreur de lecture de Mb',imb,mb
       stop
       endif
       if(iabs(ilb-lb).ne.0)then
       print *,'# erreur de lecture de Lb',ilb,lb
       stop
       endif
       if(iabs(ipp-jpp).ne.0)then
       print *,'# erreur de lecture de ipp',jpp,ipp
       stop
       endif
       c6lb(lb)=c6t
       enddo
  106  format(i3,'&',i3,'&',i3,'&',i3,'&',i3,'&',i3,'&',g16.8)
       return
       end 
