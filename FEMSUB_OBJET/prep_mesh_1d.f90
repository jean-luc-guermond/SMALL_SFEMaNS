MODULE prep_mesh_1d
  IMPLICIT NONE
CONTAINS
    SUBROUTINE load_mesh_1d(dir, fil, type_fe, mesh)
    USE input_data
    USE mod_gauss_points_1d
    USE def_type_mesh
    IMPLICIT NONE
    CHARACTER(len=64),     INTENT(IN) :: dir, fil
    INTEGER,               INTENT(IN) :: type_fe
    TYPE(mesh_type)                   :: mesh
    INTEGER, PARAMETER :: in_unit=30 
    INTEGER :: i, n, m, nn
    REAL(KIND=8) :: x0, x1, dx
    OPEN(in_unit,FILE=TRIM(ADJUSTL(dir))//&
         '/'//TRIM(ADJUSTL(fil)),FORM='formatted')
    READ(in_unit,*) mesh%me
    READ(in_unit,*) x0, x1
    ALLOCATE(mesh%jj(type_fe+1,mesh%me))
    ALLOCATE(mesh%neigh(2,mesh%me))
    DO m = 1, mesh%me
       mesh%jj(1,m) = type_fe*(m-1) + 1
       mesh%jj(2,m) = type_fe*(m-1) + type_fe+1
       nn = 2
       DO n = 2, type_fe
          nn = nn + 1
          mesh%jj(n,m) = type_fe*(m-1) + n 
       END DO
       mesh%neigh(1,m) = m+1
       mesh%neigh(2,m) = m-1
    END DO
    mesh%np=type_fe*mesh%me+1

    mesh%neigh(2,1) = 0
    mesh%neigh(1,mesh%me)= 0
    ALLOCATE(mesh%i_d(mesh%me))
    mesh%i_d = 1

    ALLOCATE(mesh%rr(1,mesh%np))
    dx = (x1-x0)/(mesh%np-1)
    DO i = 1, mesh%np
       mesh%rr(1,i) = x0+(i-1)*dx
    END DO

    mesh%mes=2
    ALLOCATE(mesh%jjs(1,mesh%mes))
    mesh%jjs(1,1) = 1
    mesh%jjs(1,mesh%mes) = mesh%np 
    ALLOCATE(mesh%sides(mesh%mes))
    READ(in_unit,*) mesh%sides

    SELECT CASE (type_fe)
    CASE(1)
       CALL GAUSS_POINTS_1d_p1(mesh)
    CASE(2)
       CALL GAUSS_POINTS_1d_p2(mesh)
    CASE(3)
       CALL GAUSS_POINTS_1d_p3(mesh)
    CASE DEFAULT
       WRITE(*,*) ' BUG load_mesh_1d: FE not programmed yet'
       STOP
    END SELECT

  END SUBROUTINE load_mesh_1d
END MODULE prep_mesh_1d
