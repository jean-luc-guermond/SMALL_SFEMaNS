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
    !===and 3 Gauss integration points
    !===Enumeration: 1  3  2
    USE def_type_mesh
    IMPLICIT NONE
    TYPE(mesh_type), TARGET :: mesh
    INTEGER,      PARAMETER       :: l_Gs = 3
    REAL(KIND=8), DIMENSION(l_Gs) :: xx, p
    INTEGER :: l, m
    REAL(KIND=8) :: zero = 0,  one = 1,  two = 2,  three = 3,  &
         five = 5, eight= 8,   nine = 9
    REAL(KIND=8) :: f1, f2, f3, df1, df2, df3, x, dhatxdx 
    f1(x) = (x - one)*x/two
    f2(x) = (x + one)*x/two
    f3(x) = (x + one)*(one - x)
    df1(x) = (two*x - one)/two
    df2(x) = (two*x + one)/two
    df3(x) = -two*x
    mesh%gauss%k_d = 1
    mesh%gauss%n_w = 3
    mesh%gauss%l_G = l_Gs
    mesh%gauss%l_Gs = 1
    xx(1) = -SQRT(three/five)
    xx(2) =  zero
    xx(3) =  SQRT(three/five) 
    p(1)  =  five/nine
    p(2)  =  eight/nine 
    p(3)  =  five/nine
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
          mesh%gauss%rj(l,m) = 1/dhatxdx
       END DO
    END DO
  END SUBROUTINE GAUSS_POINTS_1d_p2
END MODULE mod_gauss_points_1d_p2
