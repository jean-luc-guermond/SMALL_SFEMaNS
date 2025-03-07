INCLUDE 'mkl_pardiso.f90'
MODULE pardiso_solve
  USE mkl_pardiso
  TYPE pardiso_parameters
     INTEGER, POINTER :: parm(:)
     INTEGER          :: mtype
     INTEGER          :: phase
  END TYPE pardiso_parameters
  TYPE(pardiso_parameters), ALLOCATABLE, PUBLIC :: pardiso_param(:)
  PUBLIC :: allocate_pardiso_parameters, solve_pardiso
  PRIVATE
  TYPE pardiso_pt
     TYPE(MKL_PARDISO_HANDLE), ALLOCATABLE :: pt(:)
  END TYPE pardiso_pt
  TYPE(pardiso_pt), ALLOCATABLE :: my_pardiso_pt(:)
  INTEGER :: max_nb_of_matrices_in_pardiso
CONTAINS
  SUBROUTINE allocate_pardiso_parameters(maxfct)
    INTEGER :: maxfct, n
    max_nb_of_matrices_in_pardiso = maxfct
    ALLOCATE(pardiso_param(maxfct))
    ALLOCATE(my_pardiso_pt(maxfct))
    DO n = 1, maxfct
       ALLOCATE(pardiso_param(n)%parm(64))
       pardiso_param(n)%parm = 0
       !===Inputs
       pardiso_param(n)%parm(1)  = 1 !===no solver default
       pardiso_param(n)%parm(2)  = 2 !===fill-in reordering from METIS
       pardiso_param(n)%parm(4)  = 0 !===no iterative-direct algorithm
       pardiso_param(n)%parm(5)  = 0 !===no user fill-in reducing permutation
       pardiso_param(n)%parm(6)  = 0 !=== 0 solution on the first n compoments of x
       pardiso_param(n)%parm(8)  = 8 !==numbers of iterative refinement steps
       pardiso_param(n)%parm(10) = 7 !==perturb the pivot elements with 1E-7
       pardiso_param(n)%parm(11) = 1 !===use nonsymmetric permutation and scaling MPS
       pardiso_param(n)%parm(13) = 0 !===maximum weighted matching switched-off (default for symmetric)
       !===Outputs
       pardiso_param(n)%parm(14) = 0 !===number of perturbed pivots
       pardiso_param(n)%parm(18) =-1 !===number of nonzeros in the factor LU
       pardiso_param(n)%parm(19) =-1 !===Mflops for LU factorization
       pardiso_param(n)%parm(20) = 0 !===Numbers of CG Iterations
       !===Default
       pardiso_param(n)%mtype    = 1  !===Real and symmetric positive definite
       pardiso_param(n)%phase    = 33 !===Direct solve
    END DO
 
  END SUBROUTINE allocate_pardiso_parameters
  
  SUBROUTINE solve_pardiso(a,ia,ja,b,x,mnum_in,opt_memory_release)
    USE mkl_pardiso
    IMPLICIT NONE
    INTEGER, PARAMETER           :: dp = KIND(1.0D0)
    INTEGER,       INTENT(IN)    :: ia(:)
    INTEGER,       INTENT(IN)    :: ja(:)
    REAL(KIND=DP), INTENT(IN)    :: a(:)
    REAL(KIND=DP), INTENT(INOUT) :: b(:)
    REAL(KIND=DP), INTENT(OUT)   :: x(:)
    INTEGER, INTENT(IN)          :: mnum_in
    LOGICAL, OPTIONAL            :: opt_memory_release
    INTEGER       ::  maxfct, mnum, mtype, phase, n, nrhs, error, msglvl
    INTEGER       :: i, idum(1)
    REAL(KIND=DP) :: ddum(1)
    maxfct = max_nb_of_matrices_in_pardiso
    IF (ABS(mnum_in)>maxfct) THEN
       WRITE(*,*) ' BUG in solve_pardiso, ABS(mnum_in)> maxfct', mnum_in, maxfct
       STOP
    END IF
    error  = 0 !===initialize error flag
    msglvl = 0 !=== 1 print statistical information
    nrhs = 1   !===Number of right-hand-sides
    mnum = ABS(mnum_in)
    mtype=pardiso_param(mnum)%mtype !===matrix type

    IF (PRESENT(opt_memory_release)) THEN
       IF (opt_memory_release) THEN
          phase = -1 !===memory release
          n = SIZE(ia)-1   
          CALL pardiso (my_pardiso_pt(mnum)%pt, maxfct, mnum, mtype, phase, n, a, ia, ja, &
               idum, nrhs, pardiso_param(mnum)%parm, msglvl, ddum, ddum, error)
       END IF
       RETURN
    END IF
    
    IF (mnum_in<0) THEN !===Factorization
       IF (.NOT.ALLOCATED(my_pardiso_pt(mnum)%pt)) THEN
          ALLOCATE (my_pardiso_pt(mnum)%pt(64))
       END IF
       DO i = 1, 64
          my_pardiso_pt(mnum)%pt(i)%DUMMY =  0 
       END DO
       phase = 11 !===only reordering and symbolic factorization
       n = SIZE(ia)-1   
       CALL pardiso(my_pardiso_pt(mnum)%pt, maxfct, mnum, mtype, phase, n, a, ia, ja, &
            idum, nrhs, pardiso_param(mnum)%parm, msglvl, ddum, ddum, error)
       IF (error /= 0) THEN
          WRITE(*,*) 'The following ERROR was detected: ', error
          STOP
       END IF
       phase = 22 !=== only factorization
       CALL pardiso (my_pardiso_pt(mnum)%pt, maxfct, mnum, mtype, phase, n, a, ia, ja, &
            idum, nrhs, pardiso_param(mnum)%parm, msglvl, ddum, ddum, error)
       IF (error /= 0) THEN
          WRITE(*,*) 'The following ERROR was detected: ', error
          STOP
       ENDIF
    END IF

    n = SIZE(ia)-1  
    phase = pardiso_param(mnum)%phase !===Solution
    CALL pardiso (my_pardiso_pt(mnum)%pt, maxfct, mnum, mtype, phase, n, a, ia, ja, &
         idum, nrhs, pardiso_param(mnum)%parm, msglvl, b, x, error)
    IF (error /= 0) THEN
       WRITE(*,*) 'The following ERROR was detected: ', error
       STOP
    ENDIF

  END SUBROUTINE solve_pardiso
END MODULE pardiso_solve
