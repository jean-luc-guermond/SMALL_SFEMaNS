!===
!Authors: Jean-Luc Guermond, Lugi Quartapelle, Copyright 1994
!Author:  Jean-Luc Guermond, Copyright December 24th 2018
!===
MODULE mod_gauss_points_1d_p1 
  PRIVATE
  PUBLIC gauss_points_1d_p1
CONTAINS
  SUBROUTINE GAUSS_POINTS_1d_p1(mesh)
    USE def_type_mesh
    IMPLICIT NONE
    TYPE(mesh_type), TARGET :: mesh 
    REAL(KIND=8) ::  one = 1.d0,  two = 2.d0,  three = 3.d0
    REAL(KIND=8) :: f1, f2, x, dhatxdx, dx, eps=1.d-12
    REAL(KIND=8), DIMENSION(2) :: xx
    INTEGER :: l, m, ms, cote, k_dim=1
    f1(x) = (one - x)/two
    f2(x) = (x + one)/two
    mesh%gauss%k_d = 1
    mesh%gauss%n_w = 2
    mesh%gauss%n_ws = 1
    mesh%gauss%l_G = 2
    mesh%gauss%l_Gs = 1
    ALLOCATE(mesh%gauss%ww(mesh%gauss%n_w,mesh%gauss%l_G))
    ALLOCATE(mesh%gauss%dw(k_dim,mesh%gauss%n_w,mesh%gauss%l_G,mesh%me))
    ALLOCATE(mesh%gauss%rj(mesh%gauss%l_G,mesh%me))
    xx(1) = - one/SQRT(three)
    xx(2) = + one/SQRT(three)
    DO l = 1, mesh%gauss%l_G
       mesh%gauss%ww(1,l) = f1(xx(l))
       mesh%gauss%ww(2,l) = f2(xx(l))
    ENDDO
    DO m = 1, mesh%me
       dhatxdx = 2/ABS(mesh%rr(1,mesh%jj(1,m)) - mesh%rr(1,mesh%jj(2,m)))
       DO l = 1, mesh%gauss%l_G
          mesh%gauss%dw(1,1,l,m) = - one/two*dhatxdx
          mesh%gauss%dw(1,2,l,m) =   one/two*dhatxdx
          mesh%gauss%rj(l,m) = 1/dhatxdx
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
             dhatxdx = 2/ABS(mesh%rr(1,mesh%jj(1,m)) - mesh%rr(1,mesh%jj(2,m)))
             mesh%gauss%dwni(1,1,cote,ms) = - one/two*dhatxdx
             mesh%gauss%dwni(2,1,cote,ms) =   one/two*dhatxdx

             mesh%gauss%dwni(1,1,cote,ms) = (-1)**cote    *one/two*dhatxdx
             mesh%gauss%dwni(2,1,cote,ms) = (-1)**(cote+1)*one/two*dhatxdx
             dx = mesh%rr(1,mesh%jj(1,m)) - mesh%rr(1,mesh%jj(2,m))
             IF (abs(mesh%rr(1,mesh%jjsi(1, ms)) - mesh%rr(1,mesh%jj(1,m))).LE.eps*abs(dx)) THEN
                mesh%gauss%rnormsi(1, :, cote, ms) = dx/abs(dx)
                IF (cote==1) WRITE(*,*) ' BUG cote 1', cote 
             ELSE
                mesh%gauss%rnormsi(1, :, cote, ms) = -dx/abs(dx)
                IF (cote==2) WRITE(*,*) ' BUG cote 2', cote 
             END IF
          END DO
       END DO
    END IF

  END SUBROUTINE GAUSS_POINTS_1d_p1
END MODULE mod_gauss_points_1d_p1
