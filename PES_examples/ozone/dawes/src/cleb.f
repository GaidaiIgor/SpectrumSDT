      subroutine cleb(a,b,c,d,e,f,kselc,rac)
c#######################################################################
c     taro tamura, c.p.c. 1(1970)337-342.  (correcte4 version 19/1171@ #
c     calculates clebsch-gordon coefficients.                          #
c     see g.racah, phys.rev. 62(1942)438.                              #
c --- this subroutine has been modified according to the comments made #
c     by j. g. wills in c.p.c. 2(1971)381-382.                         #
c----------------------------------------------------------------------#
c     rac = clebsh (j1=a, j2=b, j3=c; m1=d, m2=e, m3=f)                #
c     kselc # 0  : test for selection rule before calculation          #
c#######################################################################
      implicit double precision (a-h,o-z)
      parameter (nfctmx=250)
      common/logfac/facl(nfctmx)
      data ntimes/0/
      if (ntimes .eq.0) call faclog
      ntimes = ntimes+1
c
      ia = ifix(sngl(2.*a))
      ib = ifix(sngl(2.*b))
      ic = ifix(sngl(2.*c))
      id = ifix(sngl(2.*d))
      ie = ifix(sngl(2.*e))
      if = ifix(sngl(2.*f))
c
      rac=0.0
      ig4=ia+ib+ic
      if(kselc.eq.0) go to 170
      if(id+ie.ne.if) go to 1000
      if(mod(ig4,2).ne.0) go to 1000
      if(ia+ib-ic.lt.0.or.ic-iabs(ia-ib).lt.0) go to 1000
      if(min0(ia-iabs(id),ib-iabs(ie),ic-iabs(if)).lt.0) go to 1000
      if(mod(ib+ie,2).ne.0.or.mod(ic+if,2).ne.0) go to 1000
  170 if(ia.eq.0.or.ib.eq.0) go to 175
      if(ic) 1000,180,185
  175 rac=1.0
      go to 1000
  180 fb=ib+1
      rac=((-1.0)**((ia-id)/2))/dsqrt(fb)
      go to 1000
  185 if(id.ne.0.or.ie.ne.0) go to 250
      ig2=ig4/2
      if(mod(ig2,2).ne.0) go to 1000
      ig=ig2/2
      i1=(ia+ib-ic)/2+1
      i2=(ic+ia-ib)/2+1
      i3=(ic+ib-ia)/2+1
      i4=ig2+2
      i5=ig+1
      i6=(ig2-ia)/2+1
      i7=(ig2-ib)/2+1
      i8=(ig2-ic)/2+1
      f1=dexp(.5*(facl(i1)+facl(i2)+facl(i3)-facl(i4))
     1         +(facl(i5)-facl(i6)-facl(i7)-facl(i8)))
      f2=ic+1
      f2=dsqrt(f2)
      s1=1-2*mod((ig2+ic)/2,2)
      rac=s1*f1*f2
      go to 1000
  250 fc2=ic+1
      iabcp=ig4/2+1
      iabc=iabcp-ic
      icab=iabcp-ib
      ibca=iabcp-ia
      iapd=(ia+id)/2+1
      iamd=iapd-id
      ibpe=(ib+ie)/2+1
      ibme=ibpe-ie
      icpf=(ic+if)/2+1
      icmf=icpf-if
      sqfclg=0.5*(dlog(fc2)-facl(iabcp+1)
     1      +facl(iabc)+facl(icab)+facl(ibca)
     2      +facl(iapd)+facl(iamd)+facl(ibpe)
     3      +facl(ibme)+facl(icpf)+facl(icmf))
      nzmic2=(ib-ic-id)/2
      nzmic3=(ia-ic+ie)/2
      nzmi=max0(0,nzmic2,nzmic3)+1
      nzmx=min0(iabc,iamd,ibpe)
      s1=(-1.0)**(nzmi-1)
      nz = nzmi
      nzm1 = nzmi - 1
      nzt1=iabc-nzm1
      nzt2=iamd-nzm1
      nzt3=ibpe-nzm1
      nzt4=nz-nzmic2
      nzt5=nz-nzmic3
      termlg=sqfclg-facl(nz)-facl(nzt1)-facl(nzt2)
     1           -facl(nzt3)-facl(nzt4)-facl(nzt5)
      rac = s1*dexp(termlg)
c
      nterm = nzmx - nzmi + 1
      a1 = 1.
      i = nterm
      do 400 j=1,nterm
      x = float((iabc-nzm1-i)*(iamd-nzm1-i)*(ibpe-nzm1-i))
      y = float((nzm1+i)*(-nzmic2+nzm1+i)*(-nzmic3+nzm1+i))
      a1 = 1. - a1*x/y
      i = i - 1
  400 continue
      rac = rac*a1
      if(dabs(rac).lt.0.000001) rac=0.0
c
 1000 return
      end
