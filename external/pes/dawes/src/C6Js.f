c xgw fix c6(:,-1) bug, search xgw below. Sep2013
      subroutine c6jmat(ipp,c6jj)
      implicit real*8(a-h,o-z)
      parameter (lmax=2,nfj=9)
      dimension c6lb(0:lmax),c6(0:lmax,0:2*lmax+1),C6mmp(4,0:3,0:3,0:2)
      dimension ci(nfj,3),ims(nfj,3),iml(nfj,3)
      dimension jnt(nfj),ijj(nfj),imj(nfj)
      dimension c6jj(nfj,nfj,0:lmax,0:2*lmax+1)
* element de matrice de O2 
C      write(9,*)'ipp=',ipp
      if(ipp.eq.1)then
      mmb=0
       else
        print *,'erreur de ipp',ipp
       stop
       endif
* Decomposition de la fonction d'onde |JMj> de l'atome  en |LSMlMs>
      l=1
      is=1
      xl=dfloat(l)
      xs=dfloat(is)
*     print *,'J Mj Ml Ms Clebsh'
      nf=0
      do jj=iabs(l-is),l+is
      xj=dfloat(jj)
      do mj=-jj,jj
      nf=nf+1
      ijj(nf)=jj
      imj(nf)=mj
      xmj=dfloat(mj)
      nt=0
      do ml=-l,l
      nt=nt+1
      xml=dfloat(ml)
      ms=mj-ml
      xms=dfloat(ms)
      call cleb(xl,xs,xj,xml,xms,xmj,1,xcleb)
*     print '(4f6.2,f16.8)',xj,xmj,xml,xms,xcleb
      ci(nf,nt)=xcleb
      ims(nf,nt)=ms
      iml(nf,nt)=ml
       enddo
      jnt(nf)=nt
*     print *
       enddo
       enddo
      nft=nf
*     print *,'nft=',nft
      if(nft.ne.nfj)then
      print *,'erreur de dimension de ci,ims et iml',nft
      stop
      endif
* Calcul des C6(mla,mpla,ipp)
      do mla=-l,l
      do mpla=-l,l
      mma=mla-mpla
      call C6r(mla,mpla,mma,ipp,mmb,c6lb)
      do lb=0,lmax
      m1=mla+l
      m2=mpla+l
      c6mmp(ipp,m1,m2,lb)=c6lb(lb) 
      enddo
      enddo
      enddo
* Calcul des C6 pour les elements de matrice <JMj| |>< | |J'Mj'>
      do i=1,nft
      jj=ijj(i)
      mj=imj(i)
      do j=1,nft
      jjp=ijj(j)
      mpj=imj(j)
C      write(9,100)jj,mj,jjp,mpj
      ni=nft-i+1
      nj=nft-j+1
      	do il=0,lmax
      	do mm=-il,il
      	imm=mm+il
      	if(imm.gt.5)print *,'erreur de dim de C6',imm
      	c6(il,imm)=0.d0
      	enddo
      	enddo
      		do it=1,jnt(i)
      ms=ims(i,it)
      mla=iml(i,it)
      m1=mla+l
      xci=ci(i,it)
      if(dabs(xci).gt.1.e-5)then 
      		do jt=1,jnt(j)
      mps=ims(j,jt)
      mpla=iml(j,jt)
      m2=mpla+l
      xcj=ci(j,jt)
      if(dabs(xcj).gt.1.e-5)then
      if(iabs(ms-mps).lt.1.e-4)then
C      write(9,99),ms,mla,mps,mpla
      	do lb=0,lmax
      ima=mla-mpla+lb
      xc6=c6mmp(ipp,m1,m2,lb)
      if(ima.lt.0.and.dabs(xc6).gt.1.e-5)then
      print *,'|Ma| >  Lb : STOP',mla,mpla,lb,xc6
      stop
      endif
*     print *,ipp,m1,m2,lb,xc6
c     xgw,  fix c6(:,-1) bug        
        if ( ima >= 0 ) then 
      c6(lb,ima)=c6(lb,ima)+xc6*xci*xcj
        end if
      	enddo
      endif
      endif
      		enddo
      endif
      		enddo
      	do il=0,lmax
      	do mm=-il,il
        imm=mm+il
        if(dabs(c6(il,imm)).gt.1.e-5)then
C      	write(9,101)il,mm,mmb,c6(il,imm)
        endif
      c6jj(ni,nj,il,imm)=c6(il,imm)
      	enddo
      	enddo
      enddo
      enddo
  99  format('Ms=',i3,' Ml=',i3,' Mps=',i3,' Mpl=',i3)
 100  format('<J=',i3,' Mj=',i3,' Jp=',i3,' Mjp=',i3,'>')
 101  format('C6(Lb=',i3,'Ma=',i3,'Mb=',i3,')=',f16.8)
      return
      end
