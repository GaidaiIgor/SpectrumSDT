      subroutine cnelr(mla,mpla,ipp,Ma,Mb,clb)
      use splin
* energie electrostatique de X(3P)+OH(X2Pi)
      implicit real*8(a-h,o-z)
      parameter (Nip=1,lmax=2)
      dimension clb(0:lmax)
      dimension xi(nr),yi(nr),b(nr),c(nr),d(nr)
      common/droo/roo
      common/dcc/ic
       if(Ma.ne.mla-mpla)then
       write(6,*)'Erreur de Ma',Ma
       stop
       endif
       do lb=1,lmax
       clb(lb)=0.d0
       enddo
       l1=2
       do l2=1,lmax
       Lb=l2
       if(iabs(Mb).le.Lb)then
       ic=ic+1
       call initspline(ic,xi,yi,b,c,d)
       cc=ispline(roo,xi,yi,b,c,d,nr)
c      read(10,108)imla,impla,jpp,iMa,iMb,iLb,cc
       imla=ll(ic,1)
       impla=ll(ic,2)
       jpp=ll(ic,3)
       iMa=ll(ic,4)
       iMb=ll(ic,5)
       iLb=ll(ic,6)
c       write(11,108) imla,impla,jpp,iMa,iMb,iLb,cc
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
  108  format(i3,'&',i3,'&',i3,'&',i3,'&',i3,'&',i3,'&',g16.8)
*      print *
       clb(l2)=cc
       else
       clb(l2)=0.d0
       endif
       enddo
       return
       end 
