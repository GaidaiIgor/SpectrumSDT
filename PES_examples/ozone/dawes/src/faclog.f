      subroutine faclog
c#######################################################################
c#    initialisation of logarithms of factorials array                 #
c#######################################################################
      implicit double precision (a-h,o-z)
      parameter (nfctmx=1001)
      common /logfac/ fct(nfctmx)
      data ntimes /0/
c
      ntimes = ntimes+1
      if (ntimes .gt. 1) return
      fct(1) = 0.d0
      do 10 i = 1,nfctmx-1
         ai = i
         fct(i+1) = fct(i)+dlog(ai)
 10   continue
c
      return
      end
