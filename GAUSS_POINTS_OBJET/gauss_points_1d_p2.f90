!===
!Authors: Jean-Luc Guermond, Lugi Quartapelle, Copyright 1994
!Author:  Jean-Luc Guermond, Copyright December 24th 2018
!===
MODULE mod_gauss_points_1d_p2 
  PRIVATE
  PUBLIC gauss_points_1d_p2
CONTAINS
  SUBROUTINE GAUSS_POINTS_1d_p2(mesh)
    !===one-dimensional element with quadratic interpolation
    !===and 3 or 4 Gauss integration points
    !===Enumeration: 1  3  2
    USE def_type_mesh
    IMPLICIT NONE
    TYPE(mesh_type), TARGET :: mesh
    INTEGER,      PARAMETER       :: l_G = 4, k_dim=1
    REAL(KIND=8), DIMENSION(l_G) :: xx, p
    REAL(KIND=8) :: zero = 0,  one = 1,  two = 2,  three = 3,  &
         five = 5, eight= 8,   nine = 9
    REAL(KIND=8) :: f1, f2, f3, df1, df2, df3, x, dhatxdx, dx, xhat, eps=1.d-12
    INTEGER :: l, m, ms, cote
    f1(x) = (x - one)*x/two
    f2(x) = (x + one)*x/two
    f3(x) = (x + one)*(one - x)
    df1(x) = (two*x - one)/two
    df2(x) = (two*x + one)/two
    df3(x) = -two*x

    mesh%gauss%k_d = 1
    mesh%gauss%n_w = 3
    mesh%gauss%l_G = l_G
    mesh%gauss%n_ws = 1
    mesh%gauss%l_Gs = 1
    ALLOCATE(mesh%gauss%ww(mesh%gauss%n_w,mesh%gauss%l_G))
    ALLOCATE(mesh%gauss%dw(k_dim,mesh%gauss%n_w,mesh%gauss%l_G,mesh%me))
    ALLOCATE(mesh%gauss%rj(mesh%gauss%l_G,mesh%me))
    
!!$    xx(1) = -SQRT(three/five)
!!$    xx(2) =  zero
!!$    xx(3) =  SQRT(three/five) 
!!$    p(1)  =  five/nine
!!$    p(2)  =  eight/nine 
!!$    p(3)  =  five/nine
!!$    
    xx(1) = -SQRT((15.d0+2*SQRT(30.d0))/35.d0)
    xx(2) = -SQRT((15.d0-2*SQRT(30.d0))/35.d0)
    xx(3) =  SQRT((15.d0-2*SQRT(30.d0))/35.d0)
    xx(4) =  SQRT((15.d0+2*SQRT(30.d0))/35.d0)
    p(1) = (18.d0-SQRT(30.d0))/36.d0
    p(2) = (18.d0+SQRT(30.d0))/36.d0
    p(3) = (18.d0+SQRT(30.d0))/36.d0
    p(4) = (18.d0-SQRT(30.d0))/36.d0
    
!!$    xx(1) = -SQRT(5+2*SQRT(10.d0/7))/3
!!$    p(1)  = (322-13*SQRT(70.d0))/(900)
!!$    xx(2) = -SQRT(5-2*SQRT(10.d0/7))/3
!!$    p(2)  = (322+13*SQRT(70.d0))/(900)
!!$    xx(3) = 0.d0
!!$    p(3) = 128.d0/225
!!$    xx(4) = SQRT(5-2*SQRT(10.d0/7))/3
!!$    p(4)  = (322+13*SQRT(70.d0))/(900)
!!$    xx(5) = SQRT(5+2*SQRT(10.d0/7))/3
!!$    p(5)  = (322-13*SQRT(70.d0))/(900)

    DO l = 1, mesh%gauss%l_G
       mesh%gauss%ww(1,l) = f1(xx(l))
       mesh%gauss%ww(2,l) = f2(xx(l))
       mesh%gauss%ww(3,l) = f3(xx(l))
    ENDDO

    DO m = 1, mesh%me
       dhatxdx = 2/ABS(mesh%rr(1,mesh%jj(1,m)) - mesh%rr(1,mesh%jj(2,m)))
       DO l = 1, mesh%gauss%l_G
          mesh%gauss%dw(1,1,l,m) = df1(xx(l))*dhatxdx
          mesh%gauss%dw(1,2,l,m) = df2(xx(l))*dhatxdx
          mesh%gauss%dw(1,3,l,m) = df3(xx(l))*dhatxdx
          mesh%gauss%rj(l,m) = p(l)/dhatxdx          
       END DO
    END DO

    IF (mesh%edge_stab) THEN
       ALLOCATE(mesh%gauss%wwsi(mesh%gauss%n_ws, mesh%gauss%l_Gs))
       mesh%gauss%wwsi = 1.d0
       ALLOCATE(mesh%gauss%rnormsi(1, mesh%gauss%l_Gs, 2, mesh%mi))
       ALLOCATE(mesh%gauss%rji(mesh%gauss%l_Gs ,mesh%mi))
       mesh%gauss%rji = 1.d0
       ALLOCATE(mesh%gauss%dwni(mesh%gauss%n_w,mesh%gauss%l_Gs,2,mesh%mi)) 
       DO ms = 1, mesh%mi
          DO cote = 1, 2
             m =  mesh%neighi(cote,ms)
             dx = mesh%rr(1,mesh%jj(1,m)) - mesh%rr(1,mesh%jj(2,m))
             IF (abs(mesh%rr(1,mesh%jjsi(1, ms)) - mesh%rr(1,mesh%jj(1,m))).LE.eps*abs(dx)) THEN
                mesh%gauss%rnormsi(1, :, cote, ms) = dx/abs(dx)
                IF (cote==1) WRITE(*,*) ' BUG cote 1', cote 
             ELSE
                mesh%gauss%rnormsi(1, :, cote, ms) = -dx/abs(dx)
                IF (cote==2) WRITE(*,*) ' BUG cote 2', cote 
             END IF
             dx = mesh%gauss%rnormsi(1, 1, cote, ms)
             dhatxdx = 2/ABS(mesh%rr(1,mesh%jj(1,m)) - mesh%rr(1,mesh%jj(2,m)))
             IF (mesh%jjsi(1,ms)==mesh%jj(1,m)) THEN
                xhat = -1
             ELSE
                xhat = 1
             END IF
             mesh%gauss%dwni(1,1,cote,ms) = df1(xhat) *dhatxdx*dx
             mesh%gauss%dwni(2,1,cote,ms) = df2(xhat) *dhatxdx*dx
             mesh%gauss%dwni(3,1,cote,ms) = df3(xhat) *dhatxdx*dx
          END DO
       END DO
    END IF

  END SUBROUTINE GAUSS_POINTS_1d_p2
END MODULE mod_gauss_points_1d_p2
