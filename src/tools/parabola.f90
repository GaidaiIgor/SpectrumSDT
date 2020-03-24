module parabola
contains
!-----------------------------------------------------------------------
!  Finds parabola using three points.
!-----------------------------------------------------------------------
  subroutine find_parabola(r_l,pot_l,r_m,pot_m,r_r,pot_r,r,pot,a)
    implicit none
    real*8 mat(3,3),r_l,pot_l,r_m,pot_m,r_r,pot_r,r,pot
    real*8 det0,deta,detb,detc,a,b,c

    mat(1,1) = r_l**2
    mat(1,2) = r_l
    mat(1,3) = 1
    mat(2,1) = r_m**2
    mat(2,2) = r_m
    mat(2,3) = 1
    mat(3,1) = r_r**2
    mat(3,2) = r_r
    mat(3,3) = 1

    det0 = det(mat)       

    mat(1,1) = pot_l
    mat(1,2) = r_l
    mat(1,3) = 1
    mat(2,1) = pot_m
    mat(2,2) = r_m
    mat(2,3) = 1
    mat(3,1) = pot_r
    mat(3,2) = r_r
    mat(3,3) = 1

    deta = det(mat)

    mat(1,1) = r_l**2
    mat(1,2) = pot_l
    mat(1,3) = 1
    mat(2,1) = r_m**2
    mat(2,2) = pot_m
    mat(2,3) = 1
    mat(3,1) = r_r**2
    mat(3,2) = pot_r
    mat(3,3) = 1

    detb = det(mat)      

    mat(1,1) = r_l**2
    mat(1,2) = r_l
    mat(1,3) = pot_l
    mat(2,1) = r_m**2
    mat(2,2) = r_m
    mat(2,3) = pot_m
    mat(3,1) = r_r**2
    mat(3,2) = r_r
    mat(3,3) = pot_r       

    detc = det(mat)    

    a = deta / det0
    b = detb / det0
    c = detc / det0    

    r  = -b/(2*a)
    pot= a*r**2 + b*r + c

  contains
    real*8 function det(a)
      real*8 a(3,3)
      det = a(1,1)*(a(2,2)*a(3,3) - a(3,2)*a(2,3)) + a(1,2)*(a(3,1)*a(2,3) - a(2,1)*a(3,3)) + a(1,3)*(a(2,1)*a(3,2) - a(3,1)*a(2,2))
    end function
  end subroutine
end module