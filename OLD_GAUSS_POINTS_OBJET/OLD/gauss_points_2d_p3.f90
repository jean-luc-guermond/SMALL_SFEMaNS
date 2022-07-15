!===
!Author: Jean-Luc Guermond, Copyright December 24th 2018
!===
MODULE mod_gauss_points_2d_p3
  PRIVATE
  PUBLIC gauss_points_2d_p3
CONTAINS
  SUBROUTINE gauss_points_2d_p3(mesh)
    USE def_type_mesh
    IMPLICIT NONE
    INTEGER, PARAMETER :: k_d=2, n_w=10, l_G=7, n_ws=4, l_Gs=4
    TYPE(mesh_type), TARGET :: mesh 
    INTEGER,                             POINTER :: me, mes
    INTEGER,      DIMENSION(:, :),       POINTER :: jj
    INTEGER,      DIMENSION(:, :),       POINTER :: js
    REAL(KIND=8), DIMENSION(:, :),       POINTER :: rr
    REAL(KIND=8), DIMENSION(:, :),       POINTER :: ww
    REAL(KIND=8), DIMENSION(:, :),       POINTER :: wws
    REAL(KIND=8), DIMENSION(:, :, :, :), POINTER :: dw
    REAL(KIND=8), DIMENSION(:,    :, :), POINTER :: rnorms
    REAL(KIND=8), DIMENSION(:, :),       POINTER :: rj
    REAL(KIND=8), DIMENSION(:, :),       POINTER :: rjs
    REAL(KIND=8), DIMENSION(:, :, :, :), POINTER :: dw_s
    REAL(KIND=8), DIMENSION(:, :, :, :), POINTER :: dwps !special!
    REAL(KIND=8), DIMENSION(:, :, :, :), POINTER :: dws  !SPECIAL!
    REAL(KIND=8), DIMENSION(k_d, n_w,  l_G)      :: dd
    REAL(KIND=8), DIMENSION(1,   n_ws, l_Gs)     :: dds
    REAL(KIND=8), DIMENSION(l_G)                 :: pp
    REAL(KIND=8), DIMENSION(l_Gs)                :: pps
    REAL(KIND=8), DIMENSION(n_w)                 :: r
    REAL(KIND=8), DIMENSION(n_ws)                :: rs
    REAL(KIND=8), DIMENSION(k_d, k_d)            :: dr
    REAL(KIND=8), DIMENSION(1,   k_d)            :: drs
    REAL(KIND=8) :: rjac, rjacs, x
    INTEGER :: m, l, k, k1, n,  ms, ls, face, cote
    REAL(KIND=8), DIMENSION(2) :: rnor, rsd ! JLG April 2009

    !TEST
    REAL(KIND=8), DIMENSION(n_w,  l_G)       :: ww_s
    !TEST
    me => mesh%me
    mes => mesh%mes
    jj => mesh%jj
    js => mesh%jjs
    rr => mesh%rr

    ! JLG (April 2009)
    IF (mesh%edge_stab) THEN
       ALLOCATE(mesh%gauss%dwni(n_w, l_Gs, 2, mesh%mi))
       ALLOCATE(mesh%gauss%rji(l_Gs, mesh%mi))
       ALLOCATE(mesh%gauss%rnormsi(k_d, l_Gs, 2, mesh%mi)) !JLG June 4 2012
       ALLOCATE(mesh%gauss%wwsi(n_ws, l_gs)) !JLG June 4 2012
    END IF
    ! JLG (April 2009)

    ALLOCATE(mesh%gauss%ww(n_w, l_g))
    ALLOCATE(mesh%gauss%wws(n_ws, l_gs))
    ALLOCATE(mesh%gauss%dw(k_d, n_w, l_G,  me ))
    ALLOCATE(mesh%gauss%rnorms(k_d,  l_Gs, mes))
    ALLOCATE(mesh%gauss%rj(l_G, me ))
    ALLOCATE(mesh%gauss%rjs(l_Gs, mes))
    !NEW Juin 8, 2004
    ALLOCATE(mesh%gauss%dw_s(k_d, n_w,  l_Gs,  mesh%mes))
    !NEW Juin 8, 2004
    ALLOCATE(mesh%gauss%dwps(1, n_ws, l_Gs, mes))
    ALLOCATE(mesh%gauss%dws(1, n_ws, l_Gs, mes))

    ww => mesh%gauss%ww
    wws => mesh%gauss%wws
    dw => mesh%gauss%dw
    rnorms => mesh%gauss%rnorms
    rj => mesh%gauss%rj
    rjs => mesh%gauss%rjs
    !NEW Mai 26, 2004
    dw_s => mesh%gauss%dw_s
    !NEW Mai 26, 2004
    dwps => mesh%gauss%dwps
    dws => mesh%gauss%dws

    mesh%gauss%k_d  = k_d
    mesh%gauss%n_w  = n_w
    mesh%gauss%l_G  = l_G
    mesh%gauss%n_ws = n_ws
    mesh%gauss%l_Gs = l_Gs

    !===Evaluate and store the values of derivatives and of the
    !===Jacobian determinant at Gauss points of all volume elements
    CALL element_2d_p3(ww, dd, pp)
    DO m = 1, me
       DO l = 1, l_G
          DO k = 1, k_d
             r = rr(k, jj(:,m))
             DO k1 = 1, k_d
                dr(k, k1) = SUM(r * dd(k1,:,l))
             ENDDO
          ENDDO
          rjac = dr(1,1)*dr(2,2) - dr(1,2)*dr(2,1)
          DO n = 1, n_w
             dw(1, n, l, m)   &
                  = (+ dd(1,n,l)*dr(2,2) - dd(2,n,l)*dr(2,1))/rjac
             dw(2, n, l, m)   &
                  = (- dd(1,n,l)*dr(1,2) + dd(2,n,l)*dr(1,1))/rjac
          ENDDO
          rj(l, m) = ABS(rjac)*pp(l) !===Sign of Jacobian determinant unknown
       ENDDO
    ENDDO

    CALL element_1d_p3(wws, dds, pps)
    !===surface elements   
    DO ms = 1, mes
       m = mesh%neighs(ms)
       !===Determine which face of the reference element is associated with ms
       DO n = 1, 3
          IF (MINVAL(ABS(js(:,ms)-jj(n,m)))==0) CYCLE
          face = n 
       END DO
       CALL element_2d_p3_boundary (face, dd, ww_s)
       DO ls = 1, l_Gs
          DO k = 1, k_d
             r = rr(k, jj(:,m))
             DO k1 = 1, k_d
                dr(k, k1) = SUM(r * dd(k1,:,ls))
             ENDDO
          ENDDO
          rjac = (dr(1,1)*dr(2,2) - dr(1,2)*dr(2,1))
          DO n = 1, n_w
             dw_s(1,n,ls,ms) = (+dd(1,n,ls)*dr(2,2) - dd(2,n,ls)*dr(2,1))/rjac
             dw_s(2,n,ls,ms) = (-dd(1,n,ls)*dr(1,2) + dd(2,n,ls)*dr(1,1))/rjac
          ENDDO
       ENDDO
    ENDDO

   DO ms = 1, mes
       DO ls = 1, l_Gs
          DO k = 1, k_d
             rs = rr(k, js(:,ms))
             drs(1, k) = SUM(rs * dds(1,:,ls))
          ENDDO
          rjacs = SQRT( drs(1,1)**2 + drs(1,2)**2 )
          rnorms(1, ls, ms) = - drs(1,2)/rjacs
          rnorms(2, ls, ms) = + drs(1,1)/rjacs
          rjs(ls, ms) = rjacs * pps(ls)
          m = mesh%neighs(ms)
          DO n = 1, n_w
             IF (MINVAL(ABS(js(:,ms)-jj(n,m)))==0) CYCLE
             face = n 
          END DO
          rs(1:k_d) = rr(:,jj(face,m)) - (rr(:,js(1,ms))+rr(:,js(2,ms)))/2
          x = SUM(rnorms(:,ls,ms)*rs(1:k_d))
          IF (x>0) THEN
             rnorms(:,ls,ms) = -rnorms(:,ls,ms)
          END IF
       ENDDO
    ENDDO

    DO ms = 1, mes  !===Evaluation of tangent gradient. Obsolete
       DO ls = 1, l_Gs
          dwps(1,:,ls,ms) = dds(1,:,ls)*pps(ls)
          dws(1,:,ls,ms)  = dds(1,:,ls)
       ENDDO
    ENDDO

    !===Cell interface (JLG, April 2009)
    IF (mesh%edge_stab) THEN
       CALL element_1d_p3(mesh%gauss%wwsi, dds, pps) !JLG June 4 2012
       DO ms = 1, mesh%mi
          DO cote = 1, 2
             m = mesh%neighi(cote,ms)
             DO n = 1, 3
                IF (MINVAL(ABS(mesh%jjsi(1:2,ms)-jj(n,m)))==0) CYCLE
                face = n 
             END DO
             CALL element_2d_p3_boundary (face, dd, ww_s)
             DO ls = 1, l_Gs
                DO k = 1, k_d
                   rs = rr(k, mesh%jjsi(:,ms))
                   drs(1, k) = SUM(rs * dds(1,:,ls))
                ENDDO
                rjacs = SQRT(drs(1,1)**2 + drs(1,2)**2)
                mesh%gauss%rji(ls, ms) = rjacs * pps(ls)
                rnor(1) =  drs(1,2)/rjacs
                rnor(2) = - drs(1,1)/rjacs
                mesh%gauss%rnormsi(1, ls, cote, ms) = rnor(1)
                mesh%gauss%rnormsi(2, ls, cote, ms) = rnor(2)
                rsd = rr(:,jj(face,m)) - (rr(:,mesh%jjsi(1,ms))+rr(:,mesh%jjsi(2,ms)))/2
                x = SUM(rnor(:)*rsd)
                IF (x>0) THEN !===Verify the normal is outward
                   rnor = -rnor
                END IF
                DO k = 1, k_d
                   r = rr(k, jj(:,m))
                   DO k1 = 1, k_d
                      dr(k, k1) = SUM(r * dd(k1,:,ls))
                   ENDDO
                ENDDO
                rjac = (dr(1,1)*dr(2,2) - dr(1,2)*dr(2,1))
                DO n = 1, n_w
                   mesh%gauss%dwni(n,ls, cote, ms)   &
                        = rnor(1)*(+ dd(1,n,ls)*dr(2,2) - dd(2,n,ls)*dr(2,1))/rjac &
                        + rnor(2)*(- dd(1,n,ls)*dr(1,2) + dd(2,n,ls)*dr(1,1))/rjac
                ENDDO
             ENDDO
          END DO
       ENDDO

       !===Some verifications
       DO ms = 1, mesh%mi
          DO cote = 1, 2
             m = mesh%neighi(cote,ms)
             DO ls = 1, l_Gs
                x = ABS(SQRT(SUM(mesh%gauss%rnormsi(:, ls, cote, ms)**2))-1.d0)
                IF (x.GE.1.d-13) THEN
                   write(*,*) 'Bug1 in normal '
                   STOP
                END IF
                rjac = (rr(1, mesh%jjsi(2,ms)) -  rr(1, mesh%jjsi(1,ms)))*mesh%gauss%rnormsi(1, ls, cote, ms) &
                     + (rr(2, mesh%jjsi(2,ms)) -  rr(2, mesh%jjsi(1,ms)))*mesh%gauss%rnormsi(2, ls, cote, ms)
                x = SQRT((rr(1, mesh%jjsi(2,ms)) -  rr(1, mesh%jjsi(1,ms)))**2 +&
                     (rr(2, mesh%jjsi(2,ms)) -  rr(2, mesh%jjsi(1,ms)))**2)
                IF (ABS(rjac/x) .GE. 1.d-10)  THEN
                   write(*,*) 'Bug2 in normal ', ms, cote, ls, ABS(rjac/x)
                   STOP
                END IF
             END DO
          END DO
       END DO
    END IF

    !===Some verifications
    DO ms = 1, mesh%mes
       m = mesh%neighs(ms)
       DO ls = 1, l_Gs
          x = ABS(SQRT(SUM(mesh%gauss%rnorms(:, ls, ms)**2))-1.d0)
          IF (x.GE.1.d-13) THEN
             write(*,*) 'Bug1 in normal '
             STOP
          END IF
          rjac = (rr(1, mesh%jjs(2,ms)) -  rr(1, mesh%jjs(1,ms)))*mesh%gauss%rnorms(1, ls, ms) &
               + (rr(2, mesh%jjs(2,ms)) -  rr(2, mesh%jjs(1,ms)))*mesh%gauss%rnorms(2, ls, ms)
          x = SQRT((rr(1, mesh%jjs(2,ms)) -  rr(1, mesh%jjs(1,ms)))**2 +&
               (rr(2, mesh%jjs(2,ms)) -  rr(2, mesh%jjs(1,ms)))**2)
          IF (ABS(rjac/x) .GE. 1.d-10)  THEN
             write(*,*) 'Bug2 in normal ', ms, cote, ls, ABS(rjac/x)
             STOP
          END IF
       END DO
    END DO
 PRINT*, 'end of gen_Gauss'

CONTAINS

 SUBROUTINE element_2d_p3 (w, d, p)
   !===triangular element with cubic interpolation
   !===and seven Gauss integration points
   !===w(n_w, l_G) : values of shape functions at Gauss points
   !===d(2, n_w, l_G) : derivatives values of shape functions at Gauss points
   !===p(l_G) : weight for Gaussian quadrature at Gauss points
   !=== 3
   !=== 7  5
   !=== 6  10 4
   !=== 1  8  9  2
   IMPLICIT NONE
   REAL(KIND=8), DIMENSION(   n_w, l_G), INTENT(OUT) :: w
   REAL(KIND=8), DIMENSION(2, n_w, l_G), INTENT(OUT) :: d
   REAL(KIND=8), DIMENSION(l_G),         INTENT(OUT) :: p
   REAL(KIND=8), DIMENSION(l_G)                      :: xx, yy
   INTEGER :: j
   REAL(KIND=8) :: one=1.d0, two=2.d0, three=3.d0, five=5.d0, nine=9.d0
   REAL(KIND=8) ::  f1, f2, f3, f4, f5, f6, f7, f8, f9, f10,   &
        df1x, df2x, df3x, df4x, df5x, df6x, df7x, df8x, df9x, df10x, &
        df1y, df2y, df3y, df4y, df5y, df6y, df7y, df8y, df9y, df10y, &
        x, y, r, a, s1, s2, t1, t2, b1, b2, area, sq

   f1(x,y) = -0.9d1/0.2d1*x**3 - 0.27d2/0.2d1*x**2*y - 0.27d2/0.2d1*x*y**2 &
        - 0.9d1/0.2d1*y**3 + 0.9d1*x**2 + 0.18d2*x*y + 0.9d1*y**2 &
        - 0.11d2/0.2d1*x - 0.11d2/0.2d1*y + 0.1d1
   f2(x,y) = 0.9d1/0.2d1*x**3 - 0.9d1/0.2d1*x**2 + x
   f3(x,y) = 0.9d1/0.2d1*y**3 - 0.9d1/0.2d1*y**2 + y
   f4(x,y) = 0.27d2/0.2d1*x**2*y - 0.9d1/0.2d1*x*y
   f5(x,y) = 0.27d2/0.2d1*x*y**2 - 0.9d1/0.2d1*x*y
   f6(x,y) = 0.27d2/0.2d1*x**2*y + 0.27d2*x*y**2 + 0.27d2 / 0.2d1*y**3 &
        - 0.45d2/0.2d1*x*y - 0.45d2/0.2d1*y**2 + 0.9d1*y
   f7(x,y) = -0.27d2/0.2d1*x*y**2 - 0.27d2/0.2d1*y**3 + 0.9d1/0.2d1*x*y &
        + 0.18d2*y**2 - 0.9d1/0.2d1*y
   f8(x,y) = 0.27d2/0.2d1*x**3 + 0.27d2*x**2*y + 0.27d2/0.2d1*x*y**2 &
        - 0.45d2/0.2d1*x**2 - 0.45d2/0.2d1*x*y + 0.9d1*x
   f9(x,y) = -0.27d2/0.2d1*x**3 - 0.27d2/0.2d1*x**2*y + 0.18d2*x**2 &
        + 0.9d1/0.2d1*x*y - 0.9d1/0.2d1*x
   f10(x,y) = -27*x**2*y - 27*x*y**2 + 27*x*y

   df1x(x,y) = -0.27d2/0.2d1*x**2 - 0.27d2*x*y - 0.27d2/0.2d1*y**2 &
        + 0.18d2*x + 0.18d2*y - 0.11d2/0.2d1
   df2x(x,y) = 0.27d2/0.2d1*x**2 - 0.9d1*x + 0.1d1
   df3x(x,y) = 0
   df4x(x,y) = 27*x*y - 0.9d1/0.2d1 *y
   df5x(x,y) = 0.27d2/0.2d1*y**2 - 0.9d1/0.2d1*y
   df6x(x,y) = 27*x*y + 27*y**2 - 0.45d2/0.2d1*y
   df7x(x,y) = -0.27d2/0.2d1*y**2 + 0.9d1/0.2d1*y
   df8x(x,y) = 0.81d2/0.2d1*x**2 + 0.54d2*x*y - 0.45d2*x &
        + 0.27d2/0.2d1*y**2 - 0.45d2/0.2d1*y + 0.9d1
   df9x(x,y) = -0.81d2/0.2d1*x**2 + 0.36d2*x - 0.27d2*x*y &
        + 0.9d1/0.2d1*y - 0.9d1/0.2d1
   df10x(x,y) = -54*x*y - 27*y**2 + 27*y

   df1y(x,y) = -0.27d2/0.2d1*x**2 - 0.27d2*x*y - 0.27d2/0.2d1* y**2 &
        + 0.18d2*x + 0.18d2*y - 0.11d2/0.2d1
   df2y(x,y) = 0.d0
   df3y(x,y) = 0.27d2/0.2d1*y**2 - 0.9d1*y + 0.1d1
   df4y(x,y) = 0.27d2/0.2d1*x**2 - 0.9d1/0.2d1*x
   df5y(x,y) = 27*x*y - 0.9d1/0.2d1*x
   df6y(x,y) = 54*x*y + 0.81d2/0.2d1*y**2 &
        - 45*y + 0.27d2/0.2d1*x**2 - 0.45d2/0.2d1*x + 0.9d1
   df7y(x,y) = -0.81d2/0.2d1*y**2 + 0.36d2*y - 0.27d2*x*y + 0.9d1/0.2d1*x - 0.9d1/0.2d1
   df8y(x,y) = 27*x**2 + 27*x*y - 0.45d2/0.2d1*x
   df9y(x,y) = -0.27d2/0.2d1*x**2 + 0.9d1/0.2d1*x
   df10y(x,y) = -27*x**2 - 54*x*y + 27*x

   !===Degree 5; 7 Points;  Stroud: p. 314, Approximate calculation of
   !===Multiple integrals (Prentice--Hall), 1971.
   area = one/two
   sq = SQRT(three*five)
   r  = one/three;                          a = area * nine/40
   s1 = (6 - sq)/21;  t1 = (9 + 2*sq)/21;  b1 = area * (155 - sq)/1200
   s2 = (6 + sq)/21;  t2 = (9 - 2*sq)/21;  b2 = area * (155 + sq)/1200
   xx(1) = r;    yy(1) = r;    p(1) = a
   xx(2) = s1;   yy(2) = s1;   p(2) = b1
   xx(3) = s1;   yy(3) = t1;   p(3) = b1
   xx(4) = t1;   yy(4) = s1;   p(4) = b1
   xx(5) = s2;   yy(5) = s2;   p(5) = b2
   xx(6) = s2;   yy(6) = t2;   p(6) = b2
   xx(7) = t2;   yy(7) = s2;   p(7) = b2

   DO j = 1, l_G
      w(1, j)    =   f1(xx(j), yy(j))
      d(1, 1, j) = df1x(xx(j), yy(j))
      d(2, 1, j) = df1y(xx(j), yy(j))
      w(2, j)    =   f2(xx(j), yy(j))
      d(1, 2, j) = df2x(xx(j), yy(j))
      d(2, 2, j) = df2y(xx(j), yy(j))
      w(3, j)    =   f3(xx(j), yy(j))
      d(1, 3, j) = df3x(xx(j), yy(j))
      d(2, 3, j) = df3y(xx(j), yy(j))
      w(4, j)    =   f4(xx(j), yy(j))
      d(1, 4, j) = df4x(xx(j), yy(j))
      d(2, 4, j) = df4y(xx(j), yy(j))
      w(5, j)    =   f5(xx(j), yy(j))
      d(1, 5, j) = df5x(xx(j), yy(j))
      d(2, 5, j) = df5y(xx(j), yy(j))
      w(6, j)    =   f6(xx(j), yy(j))
      d(1, 6, j) = df6x(xx(j), yy(j))
      d(2, 6, j) = df6y(xx(j), yy(j))
      w(7, j)    =   f7(xx(j), yy(j))
      d(1, 7, j) = df7x(xx(j), yy(j))
      d(2, 7, j) = df7y(xx(j), yy(j))
      w(8, j)    =   f8(xx(j), yy(j))
      d(1, 8, j) = df8x(xx(j), yy(j))
      d(2, 8, j) = df8y(xx(j), yy(j))
      w(9, j)    =   f9(xx(j), yy(j))
      d(1, 9, j) = df9x(xx(j), yy(j))
      d(2, 9, j) = df9y(xx(j), yy(j))
      w(10, j)   =  f10(xx(j), yy(j))
      d(1, 10, j)=df10x(xx(j), yy(j))
      d(2, 10, j)=df10y(xx(j), yy(j))
   ENDDO

 END SUBROUTINE element_2d_p3

 SUBROUTINE element_2d_p3_boundary (face, d, w)
   !===triangular element with cubic interpolation
   !===and seven Gauss integration points
   !===d(2, n_w, l_G) : derivatives values of shape functions at Gauss points
   IMPLICIT NONE
   INTEGER,                                INTENT(IN)  :: face
   REAL(KIND=8), DIMENSION(2, n_w,  l_Gs), INTENT(OUT) :: d
   REAL(KIND=8), DIMENSION(   n_w, l_Gs),  INTENT(OUT) :: w
   REAL(KIND=8), DIMENSION(l_Gs)                       :: xx, yy, gp
   INTEGER :: j
   REAL(KIND=8) :: half=0.5d0
   REAL(KIND=8) ::  f1, f2, f3, f4, f5, f6, f7, f8, f9, f10,   &
        df1x, df2x, df3x, df4x, df5x, df6x, df7x, df8x, df9x, df10x, &
        df1y, df2y, df3y, df4y, df5y, df6y, df7y, df8y, df9y, df10y, &
        x, y

   f1(x,y) = -0.9d1/0.2d1*x**3 - 0.27d2/0.2d1*x**2*y - 0.27d2/0.2d1*x*y**2 &
        - 0.9d1/0.2d1*y**3 + 0.9d1*x**2 + 0.18d2*x*y + 0.9d1*y**2 &
        - 0.11d2/0.2d1*x - 0.11d2/0.2d1*y + 0.1d1
   f2(x,y) = 0.9d1/0.2d1*x**3 - 0.9d1/0.2d1*x**2 + x
   f3(x,y) = 0.9d1/0.2d1*y**3 - 0.9d1/0.2d1*y**2 + y
   f4(x,y) = 0.27d2/0.2d1*x**2*y - 0.9d1/0.2d1*x*y
   f5(x,y) = 0.27d2/0.2d1*x*y**2 - 0.9d1/0.2d1*x*y
   f6(x,y) = 0.27d2/0.2d1*x**2*y + 0.27d2*x*y**2 + 0.27d2/0.2d1*y**3 &
        - 0.45d2/0.2d1*x*y - 0.45d2/0.2d1*y**2 + 0.9d1*y
   f7(x,y) = -0.27d2/0.2d1*x*y**2 - 0.27d2/0.2d1*y**3 + 0.9d1/0.2d1*x*y &
        + 0.18d2*y**2 - 0.9d1/0.2d1*y
   f8(x,y) = 0.27d2/0.2d1*x**3 + 0.27d2*x**2*y + 0.27d2/0.2d1*x*y**2 &
        - 0.45d2/0.2d1*x**2 - 0.45d2/0.2d1*x*y + 0.9d1*x
   f9(x,y) = -0.27d2/0.2d1*x**3 - 0.27d2/0.2d1*x**2*y + 0.18d2*x**2 &
        + 0.9d1/0.2d1*x*y - 0.9d1/0.2d1*x
   f10(x,y) = -27*x**2*y - 27*x*y**2 + 27*x*y

   df1x(x,y) = -0.27d2/0.2d1*x**2 - 0.27d2*x*y - 0.27d2/0.2d1*y**2 &
        + 0.18d2*x + 0.18d2*y - 0.11d2/0.2d1
   df2x(x,y) = 0.27d2/0.2d1*x**2 - 0.9d1*x + 0.1d1
   df3x(x,y) = 0
   df4x(x,y) = 27*x*y - 0.9d1/0.2d1 *y
   df5x(x,y) = 0.27d2/0.2d1*y**2 - 0.9d1/0.2d1*y
   df6x(x,y) = 27*x*y + 27*y**2 - 0.45d2/0.2d1*y
   df7x(x,y) = -0.27d2/0.2d1*y**2 + 0.9d1/0.2d1*y
   df8x(x,y) = 0.81d2/0.2d1*x**2 + 0.54d2*x*y - 0.45d2*x &
        + 0.27d2/0.2d1*y**2 - 0.45d2/0.2d1*y + 0.9d1
   df9x(x,y) = -0.81d2/0.2d1*x**2 + 0.36d2*x - 0.27d2*x*y &
        + 0.9d1/0.2d1*y - 0.9d1/0.2d1
   df10x(x,y) = -54*x*y - 27*y**2 + 27*y

   df1y(x,y) = -0.27d2/0.2d1*x**2 - 0.27d2*x*y - 0.27d2/0.2d1*y**2 &
        + 0.18d2*x + 0.18d2*y - 0.11d2/0.2d1
   df2y(x,y) = 0.d0
   df3y(x,y) = 0.27d2/0.2d1*y**2 - 0.9d1*y + 0.1d1
   df4y(x,y) = 0.27d2/0.2d1*x**2 - 0.9d1/0.2d1*x
   df5y(x,y) = 27*x*y - 0.9d1/0.2d1*x
   df6y(x,y) = 54*x*y + 0.81d2/0.2d1*y**2 &
        - 45*y + 0.27d2/0.2d1*x**2 - 0.45d2/0.2d1*x + 0.9d1
   df7y(x,y) = -0.81d2/0.2d1*y**2 + 0.36d2*y - 0.27d2*x*y + 0.9d1/0.2d1*x - 0.9d1/0.2d1
   df8y(x,y) = 27*x**2 + 27*x*y - 0.45d2/0.2d1*x
   df9y(x,y) = -0.27d2/0.2d1*x**2 + 0.9d1/0.2d1*x
   df10y(x,y) = -27*x**2 - 54*x*y + 27*x

   gp(1) = -SQRT((15.d0+2*SQRT(30.d0))/35.d0)
   gp(2) = -SQRT((15.d0-2*SQRT(30.d0))/35.d0)
   gp(3) =  SQRT((15.d0-2*SQRT(30.d0))/35.d0)
   gp(4) =  SQRT((15.d0+2*SQRT(30.d0))/35.d0)

   IF (face==1) THEN
      DO l = 1, l_Gs
         xx(l) = 1.d0-(gp(l)+1.d0)*half
         yy(l) = 0.d0+(gp(l)+1.d0)*half
      END DO
   ELSE IF (face==2) THEN
      DO l = 1, l_Gs
         xx(l) = 0.d0
         yy(l) = (gp(l)+1.d0)*half
      END DO
   ELSE
      DO l = 1, l_Gs
         xx(l) = (gp(l)+1.d0)*half
         yy(l) = 0.d0
      END DO
   END IF

   DO j = 1, l_Gs
      w(1, j)    =   f1(xx(j), yy(j))
      d(1, 1, j) = df1x(xx(j), yy(j))
      d(2, 1, j) = df1y(xx(j), yy(j))
      w(2, j)    =   f2(xx(j), yy(j))
      d(1, 2, j) = df2x(xx(j), yy(j))
      d(2, 2, j) = df2y(xx(j), yy(j))
      w(3, j)    =   f3(xx(j), yy(j))
      d(1, 3, j) = df3x(xx(j), yy(j))
      d(2, 3, j) = df3y(xx(j), yy(j))
      w(4, j)    =   f4(xx(j), yy(j))
      d(1, 4, j) = df4x(xx(j), yy(j))
      d(2, 4, j) = df4y(xx(j), yy(j))
      w(5, j)    =   f5(xx(j), yy(j))
      d(1, 5, j) = df5x(xx(j), yy(j))
      d(2, 5, j) = df5y(xx(j), yy(j))
      w(6, j)    =   f6(xx(j), yy(j))
      d(1, 6, j) = df6x(xx(j), yy(j))
      d(2, 6, j) = df6y(xx(j), yy(j))
      w(7, j)    =   f7(xx(j), yy(j))
      d(1, 7, j) = df7x(xx(j), yy(j))
      d(2, 7, j) = df7y(xx(j), yy(j))
      w(8, j)    =   f8(xx(j), yy(j))
      d(1, 8, j) = df8x(xx(j), yy(j))
      d(2, 8, j) = df8y(xx(j), yy(j))
      w(9, j)    =   f9(xx(j), yy(j))
      d(1, 9, j) = df9x(xx(j), yy(j))
      d(2, 9, j) = df9y(xx(j), yy(j))
      w(10, j)   =  f10(xx(j), yy(j))
      d(1, 10, j)=df10x(xx(j), yy(j))
      d(2, 10, j)=df10y(xx(j), yy(j))
   ENDDO

 END SUBROUTINE element_2d_p3_boundary

 !------------------------------------------------------------------------------

 SUBROUTINE element_1d_p3(w, d, p)
   !===one-dimensional element with linear interpolation
   !===and 3 Gauss integration points
   !===w(n_w, l_G)    : values of shape functions at Gauss points
   !===d(1, n_w, l_G) : derivatives values of shape functions at Gauss points
   !===p(l_G)         : weight for Gaussian quadrature at Gauss points
   !===Enumeration: 1  3  4  2
   IMPLICIT NONE
   REAL(KIND=8), DIMENSION(   n_ws, l_Gs), INTENT(OUT) :: w
   REAL(KIND=8), DIMENSION(1, n_ws, l_Gs), INTENT(OUT) :: d
   REAL(KIND=8), DIMENSION(l_Gs),          INTENT(OUT) :: p
   REAL(KIND=8), DIMENSION(l_Gs) :: xx
   INTEGER :: j
   REAL(KIND=8) :: f1, f2, f3, f4, df1, df2, df3, df4, x

   f1(x) = -0.1D1/0.16D2 + x/0.16D2 + 0.9D1/0.16D2*x**2 - 0.9D1/0.16D2*x**3
   f2(x) = -x/0.16D2 + 0.9D1/0.16D2*x**2 + 0.9D1/0.16D2*x**3 - 0.1D1/0.16D2
   f3(x) = -0.9D1/0.16D2*x**2 + 0.27D2/0.16D2*x**3 + 0.9D1/0.16D2 - 0.27D2/0.16D2*x
   f4(x) = -0.9D1/0.16D2*x**2 - 0.27D2/0.16D2*x**3 + 0.9D1/0.16D2 + 0.27D2/0.16D2*x

   df1(x) = 0.1D1/0.16D2 - 0.27D2/0.16D2*x**2 + 0.9D1/0.8D1*x
   df2(x)  = -0.1D1/0.16D2 + 0.27D2/0.16D2*x**2 + 0.9D1/0.8D1*x
   df3(x) = -0.27D2/0.16D2 - 0.9D1/0.8D1*x + 0.81D2/0.16D2*x**2
   df4(x) = -0.9D1/0.8D1*x - 0.81D2/0.16D2*x**2 + 0.27D2/0.16D2

   xx(1) = -SQRT((15.d0+2*SQRT(30.d0))/35.d0)
   xx(2) = -SQRT((15.d0-2*SQRT(30.d0))/35.d0)
   xx(3) =  SQRT((15.d0-2*SQRT(30.d0))/35.d0)
   xx(4) =  SQRT((15.d0+2*SQRT(30.d0))/35.d0)
   p(1) = (18.d0-SQRT(30.d0))/36.d0
   p(2) = (18.d0+SQRT(30.d0))/36.d0
   p(3) = (18.d0+SQRT(30.d0))/36.d0
   p(4) = (18.d0-SQRT(30.d0))/36.d0

   DO j = 1, l_Gs
      w(1, j) =  f1(xx(j))
      d(1, 1, j) = df1(xx(j))
      w(2, j) =  f2(xx(j))
      d(1, 2, j) = df2(xx(j))
      w(3, j) =  f3(xx(j))
      d(1, 3, j) = df3(xx(j))
      w(4, j) =  f4(xx(j))
      d(1, 4, j) = df4(xx(j))
   ENDDO
 END SUBROUTINE element_1d_p3

END SUBROUTINE gauss_points_2d_p3

END MODULE mod_gauss_points_2d_p3

