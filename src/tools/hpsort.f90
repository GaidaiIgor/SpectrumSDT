!-----------------------------------------------------------------------
!  Sorting procedure.
!-----------------------------------------------------------------------
      subroutine hpsort(n,ra,ja)
      integer n,ja(n)
      real*8 ra(n)
c       Sorts an array ra(1:n) into ascending numerical order using
c       the Heapsort algorithm.
      integer i,ir,j,l,kja
      real*8 rra
      do i=1,n
        ja(i)=i
      enddo
      if ( n<2 ) return
      l=n/2+1
      ir=n
   10 continue
        if ( l>1 ) then
      l=l-1
      rra=ra(l)
      kja=ja(l)
        else
      rra=ra(ir)
      kja=ja(ir)
      ra(ir)=ra(1)
      ja(ir)=ja(1)
      ir=ir-1
      if ( ir==1 ) then
        ra(1)=rra
        ja(1)=kja
        return
      endif
        endif
      i=l
      j=l+l
   20   if ( j<=ir ) then
          if ( j<ir ) then
        if ( ra(j)<ra(j+1) ) j=j+1
      endif
      if ( rra<ra(j) ) then
        ra(i)=ra(j)
        ja(i)=ja(j)
        i=j
        j=j+j
      else
        j=ir+1
      endif
      goto 20
      endif
      ra(i)=rra
      ja(i)=kja
      goto 10
      end subroutine


