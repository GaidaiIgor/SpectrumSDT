module index_conversion_mod

contains

!-----------------------------------------------------------------------
!  Convert global index to local index in block-cyclic distribution.
!-----------------------------------------------------------------------
  subroutine g2l(i,n,np,nb,p,il)
    implicit none
    integer :: i    ! global array index, input
    integer :: n    ! global array dimension, input
    integer :: np   ! processor array dimension, input
    integer :: nb   ! block size, input
    integer :: p    ! processor array index, output
    integer :: il   ! local array index, output
    integer :: im1
    
    im1 = i-1
    p   = mod((im1/nb),np)
    il  = (im1/(np*nb))*nb + mod(im1,nb) + 1
  end

!-----------------------------------------------------------------------
!  Convert local index to global index in block-cyclic distribution.
!-----------------------------------------------------------------------
  subroutine l2g(il,p,n,np,nb,i)
    implicit none
    integer :: il   ! local array index, input
    integer :: p    ! processor array index, input
    integer :: n    ! global array dimension, input
    integer :: np   ! processor array dimension, input
    integer :: nb   ! block size, input
    integer :: i    ! global array index, output
    integer :: ilm1

    ilm1 = il-1
    i    = (((ilm1/nb) * np) + p)*nb + mod(ilm1,nb) + 1
  end
end module
