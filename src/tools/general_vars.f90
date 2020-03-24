module general_vars
  implicit none
  real*8,allocatable::freq1(:),freq2(:),freq3(:)
  real*8,allocatable::der1(:,:),der2(:,:),der3(:,:)
  real*8,allocatable::grho2(:),sintet2(:)
  complex*16,allocatable::freq1z(:),der1z(:,:)
end module