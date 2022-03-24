!#######################################################################
      SUBROUTINE DIAG(N,A,NA,D,E)
!#######################################################################
!
! diagonalize the N by N symmetric matrix A, returning ordered
! eigenvalues in D and eigenvectors as columns of A; E is a workspace
      IMPLICIT none
      integer :: N,NA
      REAL*8 :: A(N,N),D(N),E(N)
      integer :: ierr=0

! TRED2 to tridiagonalise A, TQL2 to obtain eigenvalues and vectors,
! SORT2 to put them both in increasing order of eigenvalue
      CALL TRED2(A,N,NA,D,E)
      CALL TQL2(N,NA,D,E,A,ierr)
      if(ierr.ne.0) stop 'probleme: base d elongation'
      CALL SORT2(N,D,A,NA)

      RETURN
      END
!-----------------------------------------------------------------------
!
! ######################################################################
      SUBROUTINE TRED2(A,N,NP,D,E)
!#######################################################################
!
! this routine is copied direct from NR, except for the IMPLICIT
! statement and corresponding replacements of 0. by 0D0
!
! Householder reduction of a real, symmetric, N by N matrix A, stored in
! an NP by NP physical array.  On output, A is replaced by the
! orthogonal matrix Q effecting the transformation.  D returns the
! diagonal elements of the tridiagonal matrix, and E the off-diagonal
! elements, with E(1)=0.  Several statements, as noted in comments, can
! be omitted if only eigenvalues are to be found, in which case A
! contains no useful information on output.  Otherwise they are to be
! included.
      IMPLICIT none
      integer :: I,J,K,L,N,NP
      real*8 :: SCALE,F,G,H,HH
      REAL*8 :: A(N,N),D(N),E(N)

      DO 18 I=N,2,-1
       L=I-1
       H=0D0
       SCALE=0D0
       IF(L.GT.1) THEN
        DO 11 K=1,L
         SCALE=SCALE+ABS(A(I,K))
11      CONTINUE
        IF(SCALE.EQ.0D0) THEN
         E(I)=A(I,L)
        ELSE
         DO 12 K=1,L
          A(I,K)=A(I,K)/SCALE
          H=H+A(I,K)**2
12       CONTINUE
         F=A(I,L)
         G=-SIGN(SQRT(H),F)
         E(I)=SCALE*G
         H=H-F*G
         A(I,L)=F-G
         F=0D0
         DO 15 J=1,L
! Omit following line if finding only eigenvalues
          A(J,I)=A(I,J)/H
          G=0D0
          DO 13 K=1,J
           G=G+A(J,K)*A(I,K)
13        CONTINUE
          DO 14 K=J+1,L
           G=G+A(K,J)*A(I,K)
14        CONTINUE
          E(J)=G/H
          F=F+E(J)*A(I,J)
15       CONTINUE
         HH=F/(H+H)
         DO 17 J=1,L
          F=A(I,J)
          G=E(J)-HH*F
          E(J)=G
          DO 16 K=1,J
           A(J,K)=A(J,K)-F*E(K)-G*A(I,K)
16        CONTINUE
17       CONTINUE
        ENDIF
       ELSE
        E(I)=A(I,L)
       ENDIF
       D(I)=H
18    CONTINUE
! Omit following line if finding only eigenvalues
      D(1)=0D0
      E(1)=0D0
      DO 23 I=1,N
! Delete lines from here ...
       L=I-1
       IF(D(I).NE.0D0) THEN
        DO 21 J=1,L
         G=0D0
         DO 19 K=1,L
          G=G+A(I,K)*A(K,J)
19       CONTINUE
         DO 20 K=1,L
          A(K,J)=A(K,J)-G*A(K,I)
20       CONTINUE
21      CONTINUE
       ENDIF
! ... to here when finding only eigenvalues
       D(I)=A(I,I)
! Also delete lines from here ...
       A(I,I)=1D0
       DO 22 J=1,L
        A(I,J)=0D0
        A(J,I)=0D0
22     CONTINUE
! ... to here when finding only eigenvalues
23    CONTINUE
      RETURN
      END
!-----------------------------------------------------------------------
!
      SUBROUTINE SORT2(N,RA,V,NV)
! ######################################################################
!
! sorts an array RA of length N into ascending numerical order using the
! Heapsort algorithm, and reorders the vectors y in V(x,y)
! correspondingly; WRK is a workspace
!
! this routine is copied direct from the NR routine SORT, except for the
! IMPLICIT statement and the STOP statement, and additional code to
! reorder the vectors
!
      IMPLICIT none
      integer :: N,NV,L,I,J,K,IR
      real*8 :: RRA
      REAL*8 :: RA(N),V(N,N),WRK(N)

      IF(N.LT.1) STOP'unnatural length in SORTE'
      L=N/2+1
      IR=N
10    CONTINUE
       IF(L.GT.1) THEN
        L=L-1
        RRA=RA(L)
        DO 1 K=1,N
         WRK(K)=V(K,L)
1       CONTINUE
       ELSE
        RRA=RA(IR)
        DO 2 K=1,N
         WRK(K)=V(K,IR)
2       CONTINUE
        RA(IR)=RA(1)
        DO 3 K=1,N
         V(K,IR)=V(K,1)
3       CONTINUE
        IR=IR-1
        IF(IR.EQ.1) THEN
         RA(1)=RRA
         DO 4 K=1,N
          V(K,1)=WRK(K)
4        CONTINUE
         RETURN
        ENDIF
       ENDIF
       I=L
       J=L+L
20     IF(J.LE.IR) THEN
        IF(J.LT.IR) THEN
         IF(RA(J).LT.RA(J+1))J=J+1
        ENDIF
        IF(RRA.LT.RA(J)) THEN
         RA(I)=RA(J)
         DO 5 K=1,N
          V(K,I)=V(K,J)
5        CONTINUE
         I=J
         J=J+J
        ELSE
         J=IR+1
        ENDIF
       GO TO 20
       ENDIF
       RA(I)=RRA
       DO 6 K=1,N
        V(K,I)=WRK(K)
6      CONTINUE
      GO TO 10
      END



!#######################################################################
!#######################################################################
!
      subroutine tql2(nm,n,d,e,z,ierr)

      implicit none
      integer i,j,k,l,m,n,ii,l1,l2,nm,mml,ierr
      real*8 d(n),e(n),z(nm,n)
      real*8 c,c2,c3,dl1,el1,f,g,h,p,r,s,s2,tst1,tst2,pythag
!
!     this subroutine is a translation of the algol procedure tql2,
!     num. math. 11, 293-306(1968) by bowdler, martin, reinsch, and
!     wilkinson.
!     handbook for auto. comp., vol.ii-linear algebra, 227-240(1971).
!
!     this subroutine finds the eigenvalues and eigenvectors
!     of a symmetric tridiagonal matrix by the ql method.
!     the eigenvectors of a full symmetric matrix can also
!     be found if  tred2  has been used to reduce this
!     full matrix to tridiagonal form.
!
!     on input
!
!        nm must be set to the row dimension of two-dimensional
!!         array parameters as declared in the calling program
!          dimension statement.
!
!        n is the order of the matrix.
!
!        d contains the diagonal elements of the input matrix.
!
!        e contains the subdiagonal elements of the input matrix
!          in its last n-1 positions.  e(1) is arbitrary.
!
!        z contains the transformation matrix produced in the
!          reduction by  tred2, if performed.  if the eigenvectors
!          of the tridiagonal matrix are desired, z must contain
!          the identity matrix.
!
!      on output
!
!        d contains the eigenvalues in ascending order.  if an
!          error exit is made, the eigenvalues are correct but
!          unordered for indices 1,2,...,ierr-1.
!
!        e has been destroyed.
!
!        z contains orthonormal eigenvectors of the symmetric
!          tridiagonal (or full) matrix.  if an error exit is made,
!          z contains the eigenvectors associated with the stored
!          eigenvalues.
!
!        ierr is set to
!          zero       for normal return,
!          j          if the j-th eigenvalue has not been
!                     determined after 30 iterations.
!
!
!     questions and comments should be directed to burton s. garbow,
!     mathematics and computer science div, argonne national laboratory
!
!     this version dated august 1983.
!
!     ------------------------------------------------------------------
!
      ierr = 0
      if (n .eq. 1) go to 1001
!
      do 100 i = 2, n
  100 e(i-1) = e(i)
!
      f = 0.0d0
      tst1 = 0.0d0
      e(n) = 0.0d0
!
      do 240 l = 1, n
         j = 0
         h = abs(d(l)) + abs(e(l))
         if (tst1 .lt. h) tst1 = h
!     .......... look for small sub-diagonal element ..........
         do 110 m = l, n
            tst2 = tst1 + abs(e(m))
            if (tst2 .eq. tst1) go to 120
!     .......... e(n) is always zero, so there is no exit
!                through the bottom of the loop ..........
  110    continue
!
  120    if (m .eq. l) go to 220
  130    if (j .eq. 30) go to 1000
         j = j + 1
!     .......... form shift ..........
         l1 = l + 1
         l2 = l1 + 1
         g = d(l)
         p = (d(l1) - g) / (2.0d0 * e(l))
         r=sqrt(p*p+1.0d0)
         d(l) = e(l) / (p + sign(r,p))
         d(l1) = e(l) * (p + sign(r,p))
         dl1 = d(l1)
         h = g - d(l)
         if (l2 .gt. n) go to 145
!
         do 140 i = l2, n
  140    d(i) = d(i) - h
!
  145    f = f + h
!     .......... ql transformation ..........
         p = d(m)
         c = 1.0d0
         c2 = c
         el1 = e(l1)
         s = 0.0d0
         mml = m - l
!     .......... for i=m-1 step -1 until l do -- ..........
         do 200 ii = 1, mml
            c3 = c2
            c2 = c
            s2 = s
            i = m - ii
            g = c * e(i)
            h = c * p
            r=sqrt(p*p+e(i)*e(i))
            e(i+1) = s * r
            s = e(i) / r
            c = p / r
            p = c * d(i) - s * g
            d(i+1) = h + s * (c * g + s * d(i))
!     .......... form vector ..........
            do 180 k = 1, n
               h = z(k,i+1)
               z(k,i+1) = s * z(k,i) + c * h
               z(k,i) = c * z(k,i) - s * h
  180       continue
!
  200    continue
!
         p = -s * s2 * c3 * el1 * e(l) / dl1
         e(l) = s * p
         d(l) = c * p
         tst2 = tst1 + abs(e(l))
         if (tst2 .gt. tst1) go to 130
  220    d(l) = d(l) + f
  240 continue
!     .......... order eigenvalues and eigenvectors ..........
      do 300 ii = 2, n
         i = ii - 1
         k = i
         p = d(i)
!
         do 260 j = ii, n
            if (d(j) .ge. p) go to 260
            k = j
            p = d(j)
  260    continue
!
         if (k .eq. i) go to 300
         d(k) = d(i)
         d(i) = p
!
         do 280 j = 1, n
            p = z(j,i)
            z(j,i) = z(j,k)
            z(j,k) = p
  280    continue
!
  300 continue
!
      go to 1001
!     .......... set error -- no convergence to an
!                eigenvalue after 30 iterations ..........
 1000 ierr = l
 1001 return
      end

!#######################################################################
C
C David J. Heisterberg
C The Ohio Supercomputer Center
C 1224 Kinnear Rd.
C Columbus, OH  43212-1163
C (614)292-6036
C djh@ccl.net    djh@ohstpy.bitnet    ohstpy::djh
C
C ANALYZE
C analyze the x to y fit and the fitting matrix
C
      SUBROUTINE analyz (n, x, y, w, u)
      IMPLICIT CHARACTER (A-Z)
C
      INTEGER n
      DOUBLEPRECISION x (3, n)
      DOUBLEPRECISION y (3, n)
      DOUBLEPRECISION w (n)
      DOUBLEPRECISION u (3, 3)
C
      DOUBLEPRECISION err
      DOUBLEPRECISION wnorm
      DOUBLEPRECISION urnorm (3)
      DOUBLEPRECISION ucnorm (3)
      INTEGER i, j
C
C find the mean sum of squares error
C
      err = 0.0D0
      wnorm = 0.0D0
      DO 10000 i = 1, n
        err = err + w (i) * ((y (1, i) - x (1, i)) ** 2 +
     &                       (y (2, i) - x (2, i)) ** 2 +
     &                       (y (3, i) - x (3, i)) ** 2)
        wnorm = wnorm + w (i)
10000   CONTINUE
      err = SQRT (err / wnorm)
C
C find the row and column norms of u
C
      ucnorm (1) = u (1, 1) ** 2 + u (2, 1) ** 2 + u (3, 1) ** 2
      ucnorm (2) = u (1, 2) ** 2 + u (2, 2) ** 2 + u (3, 2) ** 2
      ucnorm (3) = u (1, 3) ** 2 + u (2, 3) ** 2 + u (3, 3) ** 2
      urnorm (1) = u (1, 1) ** 2 + u (1, 2) ** 2 + u (1, 3) ** 2
      urnorm (2) = u (2, 1) ** 2 + u (2, 2) ** 2 + u (2, 3) ** 2
      urnorm (3) = u (3, 1) ** 2 + u (3, 2) ** 2 + u (3, 3) ** 2
C
C write the error and u norms
C
      WRITE (*, *)
      DO 11000 i = 1, 3
        WRITE (*, 1000) (u (i, j), j = 1, 3), urnorm (i)
11000   CONTINUE
      WRITE (*, *)
      WRITE (*, 1000) (ucnorm (j), j = 1, 3), err
C
      RETURN
C
 1000 FORMAT (1X, 3E18.10, 4X, E18.10)
C
      END
C
C CENTER
C center a molecule about its weighted centroid or other origin
C
      SUBROUTINE center (n, x, w, io, o, t)
      IMPLICIT CHARACTER (A-Z)
C
      INTEGER n
      DOUBLEPRECISION x (3, n)
      DOUBLEPRECISION w (n)
      LOGICAL io
      DOUBLEPRECISION o (3)
      DOUBLEPRECISION t (3)
C
      DOUBLEPRECISION wnorm
      INTEGER i
C
      IF (io) THEN
        t (1) = o (1)
        t (2) = o (2)
        t (3) = o (3)
      ELSE
        t (1) = 0.0D0
        t (2) = 0.0D0
        t (3) = 0.0D0
        wnorm = 0.0D0
        DO 10000 i = 1, n
          t (1) = t (1) + x (1, i) * SQRT (w (i))
          t (2) = t (2) + x (2, i) * SQRT (w (i))
          t (3) = t (3) + x (3, i) * SQRT (w (i))
          wnorm = wnorm + SQRT (w (i))
10000     CONTINUE
        t (1) = t (1) / wnorm
        t (2) = t (2) / wnorm
        t (3) = t (3) / wnorm
      ENDIF
C
      DO 10100 i = 1, n
        x (1, i) = x (1, i) - t (1)
        x (2, i) = x (2, i) - t (2)
        x (3, i) = x (3, i) - t (3)
10100   CONTINUE
C
      RETURN
      END

C GENROT
C Generate a left rotation matrix from a normalized rotation axis
C and cosine and sine of the rotation angle.
C
C INPUT
C   a      - rotation axis
C   cost   - cosine of rotation angle
C   sint   - sine of rotation angle
C
C OUTPUT
C   u      - the rotation matrix
C
      SUBROUTINE genrot (a, cost, sint, u)
      IMPLICIT CHARACTER (A-Z)
C
      DOUBLEPRECISION a (3)
      DOUBLEPRECISION cost, sint
      DOUBLEPRECISION u (3, 3)
C
      u (1, 1) = a (1) ** 2 + (1.0D0 - a (1) ** 2) * cost
      u (2, 1) = a (1) * a (2) * (1.0D0 - cost) + a (3) * sint
      u (3, 1) = a (1) * a (3) * (1.0D0 - cost) - a (2) * sint
C
      u (1, 2) = a (1) * a (2) * (1.0D0 - cost) - a (3) * sint
      u (2, 2) = a (2) ** 2 + (1.0D0 - a (2) ** 2) * cost
      u (3, 2) = a (2) * a (3) * (1.0D0 - cost) + a (1) * sint
C
      u (1, 3) = a (1) * a (3) * (1.0D0 - cost) + a (2) * sint
      u (2, 3) = a (2) * a (3) * (1.0D0 - cost) - a (1) * sint
      u (3, 3) = a (3) ** 2 + (1.0D0 - a (3) ** 2) * cost
C
      RETURN
      END
C
C GENXYW
C generate random reference and test molecules and weights
C
      SUBROUTINE genxyw (angle, pert, n, x, y, w)
      IMPLICIT CHARACTER (A-Z)
C
      DOUBLEPRECISION angle
      DOUBLEPRECISION pert
      INTEGER n
      DOUBLEPRECISION x (3, n)
      DOUBLEPRECISION y (3, n)
      DOUBLEPRECISION w (n)
C
      DOUBLEPRECISION u (3, 3)            
      INTEGER i
C
      DOUBLEPRECISION random
      EXTERNAL random
C
      CALL ranmol (n, y)
      CALL ranrot (angle, u)
      CALL rotmol (n, y, x, u)
      CALL pertur (pert, n, x)
      CALL ranrot (angle, u)
      CALL rotmol (n, x, x, u)
C
      DO 10000 i = 1, n
        w (i) = random ()
10000   CONTINUE
C
      RETURN
      END
C
C IDENTM
C generate an identity matrix
C
      SUBROUTINE identm (n, np, u)
      IMPLICIT CHARACTER (A-Z)
C
      INTEGER n, np
      DOUBLEPRECISION u (np, n)
C
      INTEGER i, j
C
      DO 10000 j = 1, n
        DO 10010 i = 1, n
          u (i, j) = 0.0D0
10010     CONTINUE
        u (j, j) = 1.0D0
10000   CONTINUE
C
      RETURN
      END
C
C JACOBI
C Jacobi diagonalizer with sorted output.  Same calling sequence as
C EISPACK routine, but must specify nrot!
C
      SUBROUTINE jacobi (a, n, np, d, v, nrot)
      IMPLICIT CHARACTER (A-Z)
C
      INTEGER n, np, nrot
      DOUBLEPRECISION a (np, n)
      DOUBLEPRECISION d (n)
      DOUBLEPRECISION v (np, n)
C
      DOUBLEPRECISION onorm, dnorm
      DOUBLEPRECISION b, dma, q, t, c, s
      DOUBLEPRECISION atemp, vtemp, dtemp
      INTEGER i, j, k, l
C
      DO 10000 j = 1, n
        DO 10010 i = 1, n
          v (i, j) = 0.0D0
10010     CONTINUE
        v (j, j) = 1.0D0
        d (j) = a (j, j)
10000   CONTINUE
C
      DO 20000 l = 1, nrot
        dnorm = 0.0D0
        onorm = 0.0D0
        DO 20100 j = 1, n
          dnorm = dnorm + ABS (d (j))
          DO 20110 i = 1, j - 1
            onorm = onorm + ABS (a (i, j))
20110       CONTINUE
20100     CONTINUE
        IF (onorm / dnorm .LE. 0.0D0) GOTO 19999
        DO 21000 j = 2, n
          DO 21010 i = 1, j - 1
            b = a (i, j)
            IF (ABS (b) .GT. 0.0D0) THEN
              dma = d (j) - d (i)
              IF (ABS (dma) + ABS (b) .LE. ABS (dma)) THEN
                t = b / dma
              ELSE
                q = 0.5D0 * dma / b
                t = SIGN (1.0D0 / (ABS (q) + SQRT (1.0D0 + q * q)), q)
              ENDIF
              c = 1.0D0 / SQRT (t * t + 1.0D0)
              s = t * c
              a (i, j) = 0.0D0
              DO 21110 k = 1, i - 1
                atemp    = c * a (k, i) - s * a (k, j)
                a (k, j) = s * a (k, i) + c * a (k, j)
                a (k, i) = atemp
21110           CONTINUE
              DO 21120 k = i + 1, j - 1
                atemp    = c * a (i, k) - s * a (k, j)
                a (k, j) = s * a (i, k) + c * a (k, j)
                a (i, k) = atemp
21120           CONTINUE
              DO 21130 k = j + 1, n
                atemp    = c * a (i, k) - s * a (j, k)
                a (j, k) = s * a (i, k) + c * a (j, k)
                a (i, k) = atemp
21130           CONTINUE
              DO 21140 k = 1, n
                vtemp    = c * v (k, i) - s * v (k, j)
                v (k, j) = s * v (k, i) + c * v (k, j)
                v (k, i) = vtemp
21140           CONTINUE
              dtemp = c * c * d (i) + s * s * d (j) -
     &                2.0D0 * c * s * b
              d (j) = s * s * d (i) + c * c * d (j) +
     &                2.0D0 * c * s * b
              d (i) = dtemp
              ENDIF
21010       CONTINUE
21000     CONTINUE
20000   CONTINUE
19999 CONTINUE
      nrot = l
C
      DO 30000 j = 1, n - 1
        k = j
        dtemp = d (k)
        DO 30100 i = j + 1, n
          IF (d (i) .LT. dtemp) THEN
            k = i
            dtemp = d (k)
            ENDIF
30100     CONTINUE
        IF (k .GT. j) THEN
          d (k) = d (j)
          d (j) = dtemp
          DO 30200 i = 1, n
            dtemp    = v (i, k)
            v (i, k) = v (i, j)
            v (i, j) = dtemp
30200       CONTINUE
          ENDIF
30000   CONTINUE
C
      RETURN
      END
C
C PERTUR
C apply a random perturbation
C
      SUBROUTINE pertur (pert, n, x)
      IMPLICIT CHARACTER (A-Z)
C
      DOUBLEPRECISION pert
      INTEGER n
      DOUBLEPRECISION x (3, n)
C
      INTEGER i, j
C
      DOUBLEPRECISION random
      EXTERNAL random
C
      DO 10000 j = 1, 3
        DO 10010 i = 1, n
          x (j, i) = x (j, i) *
     &               (1.0D0 + 2.0D0 * pert * (0.5D0 - random ()))
10010     CONTINUE
10000   CONTINUE
C
      RETURN
      END
C
C Q2MAT
C Generate a left rotation matrix from a normalized quaternion
C
C INPUT
C   q      - normalized quaternion
C
C OUTPUT
C   u      - the rotation matrix
C
      SUBROUTINE q2mat (q, u)
      IMPLICIT NONE
C
      DOUBLEPRECISION q (0 : 3)
      DOUBLEPRECISION u (3, 3)
C
      u (1, 1) = q (0) ** 2 + q (1) ** 2 - q (2) ** 2 - q (3) ** 2
      u (2, 1) = 2.0D0 * (q (1) * q (2) - q (0) * q (3))
      u (3, 1) = 2.0D0 * (q (1) * q (3) + q (0) * q (2))
C
      u (1, 2) = 2.0D0 * (q (2) * q (1) + q (0) * q (3))
      u (2, 2) = q (0) ** 2 - q (1) ** 2 + q (2) ** 2 - q (3) ** 2
      u (3, 2) = 2.0D0 * (q (2) * q (3) - q (0) * q (1))
C
      u (1, 3) = 2.0D0 * (q (3) * q (1) - q (0) * q (2))
      u (2, 3) = 2.0D0 * (q (3) * q (2) + q (0) * q (1))
      u (3, 3) = q (0) ** 2 - q (1) ** 2 - q (2) ** 2 + q (3) ** 2
C
      RETURN
      END
C
C QTRFIT
C Find the quaternion, q, [and left rotation matrix, u] that minimizes
C
C   |qTXq - Y| ^ 2    [|uX - Y| ^ 2]
C
C This is equivalent to maximizing Re (qTXTqY).
C
C This is equivalent to finding the largest eigenvalue and corresponding
C eigenvector of the matrix
C
C   [A2   AUx  AUy  AUz ]
C   [AUx  Ux2  UxUy UzUx]
C   [AUy  UxUy Uy2  UyUz]
C   [AUz  UzUx UyUz Uz2 ]
C
C where
C
C   A2   = Xx Yx + Xy Yy + Xz Yz
C   Ux2  = Xx Yx - Xy Yy - Xz Yz
C   Uy2  = Xy Yy - Xz Yz - Xx Yx
C   Uz2  = Xz Yz - Xx Yx - Xy Yy
C   AUx  = Xz Yy - Xy Yz
C   AUy  = Xx Yz - Xz Yx
C   AUz  = Xy Yx - Xx Yy
C   UxUy = Xx Yy + Xy Yx
C   UyUz = Xy Yz + Xz Yy
C   UzUx = Xz Yx + Xx Yz
C
C The left rotation matrix, u, is obtained from q by
C
C   u = qT1q
C
C INPUT
C   n      - number of points
C   x      - test vector
C   y      - reference vector
C   w      - weight vector
C
C OUTPUT
C   q      - the best-fit quaternion
C   u      - the best-fit left rotation matrix
C   nr     - number of jacobi sweeps required
C
      SUBROUTINE qtrfit (n, x, y, w, q, u, nr)
      IMPLICIT CHARACTER (A-Z)
C
      INTEGER n
      DOUBLEPRECISION x (3, n)
      DOUBLEPRECISION y (3, n)
      DOUBLEPRECISION w (n)
      DOUBLEPRECISION q (0 : 3)
      DOUBLEPRECISION u (3, 3)
      INTEGER nr
C
      DOUBLEPRECISION xxyx, xxyy, xxyz
      DOUBLEPRECISION xyyx, xyyy, xyyz
      DOUBLEPRECISION xzyx, xzyy, xzyz
      DOUBLEPRECISION c (0 : 3, 0 : 3), v (0 : 3, 0 : 3)
      DOUBLEPRECISION d (0 : 3)
      INTEGER i
C
C generate the upper triangle of the quadratic form matrix
C
      xxyx = 0.0D0
      xxyy = 0.0D0
      xxyz = 0.0D0
      xyyx = 0.0D0
      xyyy = 0.0D0
      xyyz = 0.0D0
      xzyx = 0.0D0
      xzyy = 0.0D0
      xzyz = 0.0D0
      DO 11000 i = 1, n
        xxyx = xxyx + x (1, i) * y (1, i) * w (i)
        xxyy = xxyy + x (1, i) * y (2, i) * w (i)
        xxyz = xxyz + x (1, i) * y (3, i) * w (i)
        xyyx = xyyx + x (2, i) * y (1, i) * w (i)
        xyyy = xyyy + x (2, i) * y (2, i) * w (i)
        xyyz = xyyz + x (2, i) * y (3, i) * w (i)
        xzyx = xzyx + x (3, i) * y (1, i) * w (i)
        xzyy = xzyy + x (3, i) * y (2, i) * w (i)
        xzyz = xzyz + x (3, i) * y (3, i) * w (i)
11000   CONTINUE
C
      c (0, 0) = xxyx + xyyy + xzyz
C
      c (0, 1) = xzyy - xyyz
      c (1, 1) = xxyx - xyyy - xzyz
C
      c (0, 2) = xxyz - xzyx
      c (1, 2) = xxyy + xyyx
      c (2, 2) = xyyy - xzyz - xxyx
C
      c (0, 3) = xyyx - xxyy
      c (1, 3) = xzyx + xxyz
      c (2, 3) = xyyz + xzyy
      c (3, 3) = xzyz - xxyx - xyyy
C
C diagonalize c
C
      nr = 16
      CALL jacobi (c, 4, 4, d, v, nr)
C
C extract the desired quaternion
C
      q (0) = v (0, 3)
      q (1) = v (1, 3)
      q (2) = v (2, 3)
      q (3) = v (3, 3)
C
C generate the rotation matrix
C
  
      CALL q2mat (q, u)
C
      RETURN
      END
C
C RANDOM
C random number generator after Knuth
C
      DOUBLEPRECISION FUNCTION random ()
      IMPLICIT CHARACTER (A-Z)
C
      INTEGER ntable
      INTEGER im1, im2, im3
      INTEGER ia1, ia2, ia3
      INTEGER ic1, ic2, ic3
      DOUBLEPRECISION rm1, rm2, rm3
      PARAMETER (ntable = 97)
      PARAMETER (im1 = 714025, im2 = 214326, im3 = 139968)
      PARAMETER (ia1 = 1366,   ia2 = 3613,   ia3 = 3877)
      PARAMETER (ic1 = 150889, ic2 = 45289,  ic3 = 29573)
      PARAMETER (rm1 = 1.0D0 / im1)
      PARAMETER (rm2 = 1.0D0 / im2)
      PARAMETER (rm3 = 1.0D0 / im3)
C
      INTEGER iseed0, iseed1, iseed2, iseed3
      DOUBLEPRECISION r (ntable)
      INTEGER i
C
      COMMON /rtable/ iseed0, iseed1, iseed2, iseed3, r
      SAVE /rtable/
C
      iseed1 = MOD (ia1 * iseed1 + ic1, im1)
      iseed2 = MOD (ia2 * iseed2 + ic2, im2)
      iseed3 = MOD (ia3 * iseed3 + ic3, im3)
C
      i = 1 + (ntable * iseed3) * INT(rm3)
      random = r (i)
      r (i) = (iseed1 + iseed2 * rm2) * rm1
C
      RETURN
      END

      INTEGER FUNCTION ranget ()
      IMPLICIT CHARACTER (A-Z)
C
      INTEGER ntable
      INTEGER im1, im2, im3
      INTEGER ia1, ia2, ia3
      INTEGER ic1, ic2, ic3
      DOUBLEPRECISION rm1, rm2, rm3
      PARAMETER (ntable = 97)
      PARAMETER (im1 = 714025, im2 = 214326, im3 = 139968)
      PARAMETER (ia1 = 1366,   ia2 = 3613,   ia3 = 3877)
      PARAMETER (ic1 = 150889, ic2 = 45289,  ic3 = 29573)
      PARAMETER (rm1 = 1.0D0 / im1)
      PARAMETER (rm2 = 1.0D0 / im2)
      PARAMETER (rm3 = 1.0D0 / im3)
C
      INTEGER iseed0, iseed1, iseed2, iseed3
      DOUBLEPRECISION r (ntable)
C
      COMMON /rtable/ iseed0, iseed1, iseed2, iseed3, r
      SAVE /rtable/
C
      ranget = iseed0
C
      RETURN
      END
C
C RANINI
C seed the random number generator
C
*     SUBROUTINE ranini ()
*     IMPLICIT CHARACTER (A-Z)
C
*     INTEGER bintim (2)
C
*     CALL sys$gettim (bintim)
*     CALL ranset (IAND (bintim (1), '7FFFFFFF'X))
C
*     RETURN
*     END 
C
C RANMOL
C generate a random molecule
C
      SUBROUTINE ranmol (n, x)
      IMPLICIT CHARACTER (A-Z)
C
      INTEGER n
      DOUBLEPRECISION x (3, n)
C
      INTEGER i, j
C
      DOUBLEPRECISION random
      EXTERNAL random
C
C use nested loops to get same result for scalar/vector
C
      DO 10000 j = 1, 3
        DO 10010 i = 1, n
          x (j, i) = random ()
10010     CONTINUE
10000   CONTINUE
C
      RETURN
      END
C
C RANROT
C generate a random (almost) rotation matrix
C
      SUBROUTINE ranrot (angle, u)
      IMPLICIT CHARACTER (A-Z)
C
      DOUBLEPRECISION angle
      DOUBLEPRECISION u (3, 3)
C
      DOUBLEPRECISION a (3)
      DOUBLEPRECISION anorm
      DOUBLEPRECISION theta
      DOUBLEPRECISION pi
C
      DOUBLEPRECISION random
      EXTERNAL random
C
      pi = 4.0D0 * ATAN (1.0D0)
C
      a (1) = 0.5D0 - random ()
      a (2) = 0.5D0 - random ()
      a (3) = 0.5D0 - random ()
C
      anorm = SQRT (a (1) ** 2 + a (2) ** 2 + a (3) ** 2)
      a (1) = a (1) / anorm
      a (2) = a (2) / anorm
      a (3) = a (3) / anorM
C
      theta = angle * pi * random ()
C
      CALL genrot (a, COS (theta), SIN (theta), u)
C
      RETURN
      END
      SUBROUTINE ranset (iseed)
      IMPLICIT CHARACTER (A-Z)
C
      INTEGER ntable
      INTEGER im1, im2, im3
      INTEGER ia1, ia2, ia3
      INTEGER ic1, ic2, ic3
      DOUBLEPRECISION rm1, rm2, rm3
      PARAMETER (ntable = 97)
      PARAMETER (im1 = 714025, im2 = 214326, im3 = 139968)
      PARAMETER (ia1 = 1366,   ia2 = 3613,   ia3 = 3877)
      PARAMETER (ic1 = 150889, ic2 = 45289,  ic3 = 29573)
      PARAMETER (rm1 = 1.0D0 / im1)
      PARAMETER (rm2 = 1.0D0 / im2)
      PARAMETER (rm3 = 1.0D0 / im3)
C
      INTEGER iseed
C
      INTEGER iseed0, iseed1, iseed2, iseed3
      DOUBLEPRECISION r (ntable)
      INTEGER i
C
      COMMON /rtable/ iseed0, iseed1, iseed2, iseed3, r
      SAVE /rtable/
C
      iseed0 = iseed
      iseed1 = MOD (iseed0, im1)
      iseed2 = MOD (iseed0, im2)
      iseed3 = MOD (iseed0, im3)
C
      DO 10000 i = 1, ntable
        iseed1 = MOD (ia1 * iseed1 + ic1, im1)
        iseed2 = MOD (ia2 * iseed2 + ic2, im2)
        iseed3 = MOD (ia3 * iseed3 + ic3, im3)
        r (i) = (iseed1 + iseed2 * rm2) * rm1
10000   CONTINUE
C
      RETURN
      END
C
C ROTMOL
C rotate a molecule
C
      SUBROUTINE rotmol (n, x, y, u)
      IMPLICIT CHARACTER (A-Z)
C
      INTEGER n
      DOUBLEPRECISION x (3, n)
      DOUBLEPRECISION y (3, n)
      DOUBLEPRECISION u (3, 3)
C
      DOUBLEPRECISION yx, yy, yz
      INTEGER i
C
      DO 10000 i = 1, n
        yx = u(1, 1) * x(1, i) + u(1, 2) * x(2, i) + u(1, 3) * x(3, i)
        yy = u(2, 1) * x(1, i) + u(2, 2) * x(2, i) + u(2, 3) * x(3, i)
        yz = u(3, 1) * x(1, i) + u(3, 2) * x(2, i) + u(3, 3) * x(3, i)
C
        y (1, i) = yx
        y (2, i) = yy
        y (3, i) = yz
10000   CONTINUE
C
      RETURN
      END

        SUBROUTINE LPMN(MM,M,N,X,PM,PD)
C
C       =====================================================
C       Purpose: Compute the associated Legendre functions 
C                Pmn(x) and their derivatives Pmn'(x)
C       Input :  x  --- Argument of Pmn(x)
C                m  --- Order of Pmn(x),  m = 0,1,2,...,n
C                n  --- Degree of Pmn(x), n = 0,1,2,...,N
C                mm --- Physical dimension of PM and PD
C       Output:  PM(m,n) --- Pmn(x)
C                PD(m,n) --- Pmn'(x)
C       =====================================================
C
        !IMPLICIT DOUBLE PRECISION (P,X)
        INTEGER I,J,N,M,MM
        !DIMENSION PM(0:MM,0:MM),PD(0:MM,0:MM)
        DOUBLE PRECISION PM(0:MM,0:MM),PD(0:MM,0:MM),X,LS,XQ,XS

        DO 10 I=0,N
        DO 10 J=0,M
           PM(J,I)=0.0D0
10         PD(J,I)=0.0D0
        PM(0,0)=1.0D0
        IF (DABS(X).EQ.1.0D0) THEN
           DO 15 I=1,N
              PM(0,I)=X**I
15            PD(0,I)=0.5D0*I*(I+1.0D0)*X**(I+1)
           DO 20 J=1,N
           DO 20 I=1,M
              IF (I.EQ.1) THEN
                 PD(I,J)=1.0D+300
              ELSE IF (I.EQ.2) THEN
                 PD(I,J)=-0.25D0*(J+2)*(J+1)*J*(J-1)*X**(J+1)
              ENDIF
20         CONTINUE
           RETURN
        ENDIF
        LS=1
        IF (DABS(X).GT.1.0D0) LS=-1
        XQ=DSQRT(LS*(1.0D0-X*X))
        XS=LS*(1.0D0-X*X)
        DO 30 I=1,M
30         PM(I,I)=-LS*(2.0D0*I-1.0D0)*XQ*PM(I-1,I-1)
        DO 35 I=0,M
35         PM(I,I+1)=(2.0D0*I+1.0D0)*X*PM(I,I)
        DO 40 I=0,M
        DO 40 J=I+2,N
           PM(I,J)=((2.0D0*J-1.0D0)*X*PM(I,J-1)-
     &             (I+J-1.0D0)*PM(I,J-2))/(J-I)
40      CONTINUE
        PD(0,0)=0.0D0
        DO 45 J=1,N
45         PD(0,J)=LS*J*(PM(0,J-1)-X*PM(0,J))/XS
        DO 50 I=1,M
        DO 50 J=I,N
           PD(I,J)=LS*I*X*PM(I,J)/XS+(J+I)
     &             *(J-I+1.0D0)/XQ*PM(I-1,J)
50      CONTINUE
        RETURN
        END


