INCLUDE 'mkl_pardiso.f90'
MODULE pardiso_solve_old
INTEGER :: max_nb_of_solvers_in_pardiso
CONTAINS
  SUBROUTINE solve_pardiso(a,ia,ja,b,x,mnum_in,opt_maxfct_in,opt_memory_release)
    USE mkl_pardiso
    IMPLICIT NONE
    TYPE pardiso_pt
       TYPE(MKL_PARDISO_HANDLE), ALLOCATABLE :: pt(:)
    END type pardiso_pt
    type(pardiso_pt), ALLOCATABLE, SAVE :: my_pardiso_pt(:)
    INTEGER, PARAMETER  :: dp = KIND(1.0D0)
    INTEGER,       INTENT(IN) :: ia( : )
    INTEGER,       INTENT(IN) :: ja( : )
    REAL(KIND=DP), INTENT(IN) :: a( : )
    REAL(KIND=DP), INTENT(INOUT) :: b( : )
    REAL(KIND=DP), INTENT(OUT):: x( : )
    INTEGER, INTENT(IN) :: mnum_in
    INTEGER, OPTIONAL   :: opt_maxfct_in
    LOGICAL, OPTIONAL   :: opt_memory_release

    !.. Internal solver memory pointer 
    !TYPE(MKL_PARDISO_HANDLE), ALLOCATABLE,  SAVE :: pt(:)
    INTEGER, ALLOCATABLE,                   SAVE :: iparm(:)
    INTEGER,                                SAVE :: maxfct
    LOGICAL,                                SAVE :: once=.TRUE.
    INTEGER :: mnum, mtype, phase, n, nrhs, error, msglvl
    INTEGER :: i, idum(1)
    REAL(KIND=DP) :: ddum(1)

    IF (once) THEN
       once = .FALSE.
       IF (.NOT.PRESENT(opt_maxfct_in)) THEN
          maxfct = max_nb_of_solvers_in_pardiso
       ELSE
          IF (opt_maxfct_in.LE.0) THEN
             maxfct = 1
          ELSE
             maxfct = opt_maxfct_in
          END IF
       END IF
       ALLOCATE(my_pardiso_pt(maxfct))
       ALLOCATE( iparm ( 64 ) )

       DO i = 1, 64
          iparm(i) = 0
       END DO
       iparm(1) = 1 ! no solver default
       iparm(2) = 2 ! fill-in reordering from METIS
       iparm(4) = 0 ! no iterative-direct algorithm
       iparm(5) = 0 ! no user fill-in reducing permutation
       iparm(6) = 0 ! =0 solution on the first n compoments of x
       iparm(8) = 9 ! numbers of iterative refinement steps
       iparm(10) = 7 ! perturb the pivot elements with 1E-13
       iparm(11) = 1 ! use nonsymmetric permutation and scaling MPS
       iparm(13) = 0 ! maximum weighted matching algorithm is switched-off (default for symmetric).
       ! Try iparm(13) = 1 in case of inappropriate accuracy
       iparm(14) = 0 ! Output: number of perturbed pivots
       iparm(18) = -1 ! Output: number of nonzeros in the factor LU
       iparm(19) = -1 ! Output: Mflops for LU factorization
       iparm(20) = 0 ! Output: Numbers of CG Iterations
    END IF

    IF (ABS(mnum_in)> maxfct) THEN
       WRITE(*,*) ' BUG in solve_pardiso, ABS(mnum_in)> maxfct'
       STOP
    END IF

    error  = 0 ! initialize error flag
    msglvl = 0 ! 1 print statistical information
    nrhs = 1 
    mtype= 11  ! Nonsymmetric matrix

    IF (mnum_in<0) THEN !Factorization
       mnum = ABS(mnum_in)
       IF (.NOT.ALLOCATED(my_pardiso_pt(mnum)%pt)) THEN
          ALLOCATE (my_pardiso_pt(mnum)%pt(64))
       END IF
       DO i = 1, 64
          my_pardiso_pt(mnum)%pt(i)%DUMMY =  0 
       END DO

!!$       IF (present(opt_memory_release)) THEN
!!$          IF (opt_memory_release) THEN
!!$             phase = -1 ! memory release
!!$             n = SIZE(ia)-1   
!!$             CALL pardiso (my_pardiso_pt(mnum)%pt, maxfct, mnum, mtype, phase, n, a, ia, ja, &
!!$                  idum, nrhs, iparm, msglvl, ddum, ddum, error)
!!$             write(*,*) 'releasing memory done'
!!$          END IF
!!$       END IF

       phase = 11 !===Analysis. Only reordering and symbolic factorization
       n = SIZE(ia)-1   
       CALL pardiso(my_pardiso_pt(mnum)%pt, maxfct, mnum, mtype, phase, n, a, ia, ja, &
            idum, nrhs, iparm, msglvl, ddum, ddum, error)
       !WRITE(*,*) 'Reordering completed ... '
       IF (error /= 0) THEN
          WRITE(*,*) 'The following ERROR was detected: ', error
          STOP
       END IF
       !WRITE(*,*) 'Number of nonzeros in factors = ',iparm(18)
       !WRITE(*,*) 'Number of factorization MFLOPS = ',iparm(19)
       !.. Factorization.
       phase = 22 !===Numerical factorization only
       CALL pardiso (my_pardiso_pt(mnum)%pt, maxfct, mnum, mtype, phase, n, a, ia, ja, &
            idum, nrhs, iparm, msglvl, ddum, ddum, error)
       !WRITE(*,*) 'Factorization completed ... '
       IF (error /= 0) THEN
          WRITE(*,*) 'The following ERROR was detected: ', error
          STOP
       ENDIF
    END IF

    mnum = ABS(mnum_in)
    n = SIZE(ia)-1  
    !.. Back substitution and iterative refinement
    iparm(8) = 2 !===max numbers of iterative refinement steps
    phase = 33 !===only back substitution
    CALL pardiso (my_pardiso_pt(mnum)%pt, maxfct, mnum, mtype, phase, n, a, ia, ja, &
         idum, nrhs, iparm, msglvl, b, x, error)
    !WRITE(*,*) 'Solve completed ... '
    IF (error /= 0) THEN
       WRITE(*,*) 'The following ERROR was detected: ', error
       STOP
    ENDIF


    IF (present(opt_memory_release)) THEN
       IF (opt_memory_release) THEN
          phase = -1 ! memory release
          n = SIZE(ia)-1   
          CALL pardiso (my_pardiso_pt(mnum)%pt, maxfct, mnum, mtype, phase, n, a, ia, ja, &
               idum, nrhs, iparm, msglvl, ddum, ddum, error)
          !write(*,*) 'releasing memory done'
       END IF
    END IF


  END SUBROUTINE solve_pardiso
END MODULE pardiso_solve_old
