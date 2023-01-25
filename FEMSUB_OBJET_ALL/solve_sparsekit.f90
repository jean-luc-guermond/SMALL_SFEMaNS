MODULE solve 

  IMPLICIT NONE

  PRIVATE
  PUBLIC  :: solve_spkit, cancel_mem_precond
  PUBLIC  :: matrice_bloc

  TYPE matrice_bloc
      REAL(KIND=8), POINTER, DIMENSION(:) :: aa
      INTEGER,      POINTER, DIMENSION(:) :: ia, ja
  END TYPE matrice_bloc


  !INTEGER, PARAMETER :: max_save = 20

  INTEGER, PARAMETER :: max_krylov = 100
  INTEGER, PARAMETER :: std_krylov = 15
  INTEGER, PARAMETER :: std_maxit  = 500

  REAL(KIND=8), PARAMETER:: std_eps_rel = 1.d-6
  REAL(KIND=8), PARAMETER:: std_eps_abs = 1.d-9

  
  TYPE memo_type
     
     !  ju_size  == n_syst
     !  jau_size == n_work
     !  au_size  == n_work (meth 1-> 4) nwork+n_syst (meth 5->6)

     INTEGER :: ju_size
     INTEGER :: jau_size
     INTEGER :: au_size
     INTEGER     , DIMENSION(:), POINTER :: ju, jau
     REAL(KIND=8), DIMENSION(:), POINTER :: au
  END TYPE memo_type

  !  Defined as static vector, in order to simplify acces and to improve speed

  !TYPE(memo_type), DIMENSION(max_save) , SAVE :: memo 
  TYPE(memo_type), DIMENSION(:), ALLOCATABLE, SAVE :: memo
  LOGICAL, SAVE :: once=.TRUE.

!  EXTERNAL cg

CONTAINS
  !
  !********************************************************************
  !
  SUBROUTINE solve_spkit(ipar, fpar, ia, ja, a, rhs, sol, isave, max_save_in)

    !-----------------------------------------------------------------------
    !     Program for ilu preconditioned gmres.
    !     this program solve a  linear system generated by
    !     a sparse matrix.
    !     Different methods available:
    !     ILU0, MILU0, ILUT, ILUK, ILUD  with different values of tol and lfil
    !     ( from cheaper to more expensive preconditioners)
    !     The more accurate the preconditioner the fewer iterations 
    !     are required in pgmres, in general. 
    !-----------------------------------------------------------------------
    !
    !     ipar(1) = 1     ----> ILU0
    !             = 2     ----> MILU0
    !             = 3     ----> ILUT
    !             = 4     ----> ILUTP
    !             = 5     ----> ILUK
    !             = 6     ----> ILUD
    !             = 7     ----> GMRES or CG
    !
    !     ipar(2) = lfil  ----> If lfil < 0, default value is chosen. 
    !                     ----> Needed for ILUT, ILUTP: lfil is the number 
    !                           of non zero entries/row in the 
    !                           LU decomposition (Default=15).
    !                     ----> Needed for ILUK: Level of fill used in
    !                           the LU decomposition (Default=6).
    !     ipar(3)=n_work  ----> Work space for the ILU decomposition
    !                           If n_work < 0, default value is chosen.
    !     ipar(4) = iout  ----> Output unit number for printing intermediate
    !                           results. If iout<= 0 nothing is printed.
    !     ipar(5) = im    ----> Size of Krilov space
    !     ipar(6)=maxits  ----> Maximal number of iterations 
    !     ipar(7) = 1     ----> GMRES is used 
    !             = 2     ----> CG is used
    !
    !     fpar(1)=eps_rel ----> Stopping criterion: 
    !                           ||current residual||/||initial residual|| <= eps
    !                           or ||current residual|| <= eps_absolute
    !                           Euclidian norm is used (If eps<0, default set).
    !     fpar(2)=eps_abs ----> Stopping criterion: 
    !     fpar(3) = tol   ----> Threshold for dropping term in ILU factorization
    !                           If tol < 0, default value chosen.
    !     fpar(4) = alpha ----> Used only by ILUD: Diagonal conpensation parameter 
    !                            0 =< alpha <= 1 (IF alpha < 0, default is chosen).
    !
    !     ia              ----> Pointers to the beginning of each row of a
    !     ja, a           ----> Matrix a stored in CSR.
    !     rhs             ----> Right Hand Side
    !     sol             ----> Input: initial guess; Output: solution.
    !
    !
    !     isave           ----> OPTIONAL: integer 
    !                           If isave < 0, the preconditioner is saved
    !                                         and referred to with index |isave|.
    !                           If isave > 0, the saved-preconditioner is  
    !                                         reused.
    !-----------------------------------------------------------------------

    INTEGER,      DIMENSION(:), INTENT(IN)    :: ipar, ia
    INTEGER,      DIMENSION(:), INTENT(INOUT) :: ja
    REAL(KIND=8), DIMENSION(:), INTENT(IN)    :: fpar, a, rhs
    REAL(KIND=8), DIMENSION(:), INTENT(INOUT) :: sol 
    INTEGER,      OPTIONAL    , INTENT(IN)    :: isave, max_save_in

    INTEGER      :: k, index, index_save, iout
    INTEGER      :: meth, ierr, lfil, n_syst, n_work, im, maxits, n_bloc
    REAL(KIND=8) :: tol, permtol, eps_rel, eps_abs, alpha 
    INTEGER,      DIMENSION(16)     :: ipar_loc
    REAL(KIND=8), DIMENSION(16)     :: fpar_loc
    INTEGER, SAVE                   :: max_save

    INTEGER     , DIMENSION(:), POINTER     :: jau, ju 
    REAL(KIND=8), DIMENSION(:), POINTER     :: au
    REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: vv

    EXTERNAL user_time 
    REAL(KIND=8) :: dummy, t1, t2, user_time
    LOGICAL :: test_op

    IF (once) THEN
       IF (PRESENT(max_save_in)) THEN
          max_save = max_save_in
          ALLOCATE(memo(max_save))
          memo(:) % ju_size = 0
          once = .false.
       ELSE
          max_save = 0
       END IF

    END IF

    n_syst = SIZE(ia)-1  ! Size of linear system
    index = 0

    !-----CHECK PARAMETERS IPAR - FPAR

    IF((ipar(1) <= 0) .OR. (ipar(1) > 6)) THEN
       meth = 3
       alpha = -1.
       tol =  -1.
       n_work = -1
       write(*,*) 'Warning: method not existent -> default parameters used'
       write(*,*) 'Method:', meth

    ELSE

       meth = ipar(1)
       lfil = ipar(2)      

       n_work = ipar(3)

       IF(ipar(4) /= 6) THEN
          iout = 10
          INQUIRE(unit=iout,OPENED=test_op)
          IF(.NOT.test_op) THEN 
             OPEN(iout, file = 'solve_spkit.out', status = 'unknown')
          END IF
       ELSE
          iout = ipar(4)
       END IF

       IF((ipar(5) <= 0) .OR. (ipar(5) >= max_krylov)) THEN
          WRITE(iout,*)'Krylov space set to', std_krylov
          im = std_krylov
       ELSE
          im = ipar(5)
       ENDIF

       IF(ipar(6).LE.0) THEN
          WRITE(iout,*)'MAXIMUM Iteration number set to', std_maxit
          maxits = std_maxit
       ELSE
          maxits = ipar(6)
       END IF

       IF(fpar(1).LE.0) THEN
          WRITE(iout,*)'Relative epsilon set to', std_eps_rel
          eps_rel = std_eps_rel 
       ELSE
          eps_rel = fpar(1) 
       END IF

       IF(fpar(2).LE.0) THEN
          WRITE(iout,*)'Absolute epsilon set to', std_eps_abs
          eps_abs = std_eps_abs 
       ELSE
          eps_abs = fpar(2) 
       END IF

       tol = fpar(3)
       alpha = fpar(4)

    ENDIF

    !  Test the RHS for NULL SOLUTION
    IF (DOT_PRODUCT(rhs,rhs).le.eps_abs**2) THEN
       sol =0.d0
       IF (.NOT.PRESENT(isave)) RETURN
       IF (isave >= 0) RETURN
    END IF

    ! Tolerance ratio used to determine wether or not to permute two columns

    permtol = 0.05 

    ! Permuting is done within the diagonal blocs of size n_bloc. Useful only if
    ! several degrees of freedom of the PDE are solved at once.

    n_bloc = n_syst     

    !  Time to precondition

    t1 = user_time(dummy)

    !  Check the preconditioner status

    IF(PRESENT(isave)) THEN

       IF(ABS(isave) > max_save) THEN
          WRITE(*,*)'Error: isave,', isave,' out of range ! (maximum is ', max_save, ')'
          STOP
       ENDIF

       IF(isave >= 1) THEN
          ! Load saved preconditioner

          IF(ASSOCIATED(memo(isave) % ju)) THEN
             IF(n_syst /= memo(isave) % ju_size) THEN
                WRITE(iout,*) ' SOLVE_SPK: incoherence in the size of the linear system'
                WRITE(iout,*) ' n_syst = ',n_syst, ' ju_size', memo(isave) % ju_size
                STOP
             ENDIF
             WRITE(iout,*) ' Loading preconditioner'
             n_work = memo(isave) % jau_size
             ju  => memo(isave)%ju
             jau => memo(isave)%jau
             au  => memo(isave)%au
             index = 1
          ELSE
             WRITE(iout,*) 'SOLVE_SPK: The problem referred to by isave=',isave
             WRITE(iout,*) 'has not yet been preconditioned.'
             WRITE(iout,*) 'Use a negative integer first.'
             STOP
          ENDIF

       ELSE

          ! Preconditioning

          index_save = -isave

          IF(memo(index_save) % ju_size /=0) THEN
             WRITE(iout,*)'Place already used! ',&
                          '- Deleting old PRECONDITIONER'
             CALL cancel_mem_precond(index_save)
          ENDIF
          CALL precond_it(meth)

       ENDIF

    ELSE

       ! Preconditioning

       CALL precond_it(meth)

    ENDIF

    t1 = user_time(dummy) - t1

    !     call PGMRES

    t2 = user_time(dummy)

    IF (ipar(7) .EQ. 1) THEN
!                   == PGMRES ==
       ALLOCATE(vv(n_syst*(im+1)))
       CALL pgmres(n_syst, im, rhs, sol, vv, eps_rel, eps_abs, & 
            maxits, iout, a, ja, ia, au, jau, ju, ierr)
    ELSE IF (ipar(7) .EQ. 2) THEN
!                   == Conjugate Gradient ==
       WRITE(*,*) ' NOT PROGRAMMED '
       STOP
       ipar_loc(2) = 2
       ipar_loc(3) = 1
       ipar_loc(4) = 5*n_syst ! size of work space 
       ipar_loc(5) = 2 
       ipar_loc(6) = maxits 
       fpar_loc(1) = eps_rel 
       fpar_loc(2) = eps_abs 
       ALLOCATE(vv(5*n_syst))
!       CALL runrc(n_syst, rhs, sol, ipar_loc, fpar_loc, vv, sol, &   
!                  a, ja, ia, au, jau, ju, cg)               !  guess
       IF (ipar_loc(1) .EQ. 0)  ierr = 0
       IF (ipar_loc(1) .EQ. -1) ierr = 1 
    ELSE
       WRITE(*,*) ' SOLVE_SPKIT: solver pointed to by ipar(7) not programmed yet' 
       STOP
    END IF
 

    t2 = user_time(dummy) - t2

    select case(ierr)

    case(1)
       WRITE(iout,*) ' SOLVE_SPKIT: Convergence not achieved in ', maxits,' iterations'
       stop

    case(-1)
       WRITE(iout,*) ' SOLVE_SPKIT: Initial guess was OK'

    case(0)
       WRITE(iout,*) ' SOLVE_SPKIT: Convergence is OK'

    case DEFAULT
       WRITE(iout,*) ' SOLVE_SPKIT: Undetermined problem'

    end select

    WRITE(iout,*) ' Elimination time  = ', t2, ' Total time = ', t1+t2 
    WRITE(iout,*) '     ' 

    !  Save Preconditioner

    IF (PRESENT(isave)) THEN
       IF (isave .LT. 0) THEN

       index_save = -isave
       memo(index_save) % ju_size  = n_syst
       memo(index_save) % jau_size = n_work
       memo(index_save) % au_size  = SIZE(au)
       ALLOCATE(memo(index_save) % ju(n_syst))
       ALLOCATE(memo(index_save) % jau(n_work))
       ALLOCATE(memo(index_save) % au(memo(index_save) % au_size))

       memo(index_save) % ju  = ju
       memo(index_save) % jau = jau
       memo(index_save) % au  = au

       ENDIF
    ENDIF

    IF(ASSOCIATED(ju)) THEN
       IF (index.EQ.1) THEN 
          NULLIFY(ju,jau,au)
       ELSE
          DEALLOCATE(ju,jau,au)
       END IF
    END IF

    DEALLOCATE(vv)

  CONTAINS
    !
    !----------------------------------------------------------
    !
    SUBROUTINE precond_it(meth)

      INTEGER, INTENT(IN) :: meth

      INTEGER     , DIMENSION(:), ALLOCATABLE :: iw, iperm, levs 
      REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: wk

      ALLOCATE(ju(n_syst))
      
      SELECT CASE(meth)

      CASE(1)
         WRITE (iout,*) ' +++++ ILU(0) Preconditioner ++++ '
         n_work = SIZE(a) + n_syst
         ALLOCATE(au(n_work), jau(n_work), iw(n_syst))
         CALL ilu0 (n_syst, a, ja, ia, au, jau, ju, iw, ierr)
         IF (ierr.NE.0) THEN
            WRITE(iout,*) ' SOLVE_SPKIT: zero pivot at step ',ierr
            STOP
         END IF

      CASE(2)
         WRITE (iout,*) ' +++++ MILU(0) Preconditioner ++++ '
         n_work = SIZE(a) + n_syst
         ALLOCATE(au(n_work), jau(n_work), iw(n_syst))
         CALL milu0 (n_syst, a, ja, ia, au, jau, ju, iw, ierr)
         IF (ierr.NE.0) THEN
            WRITE(iout,*) ' SOLVE_SPKIT: zero pivot at step ',ierr
            STOP
         END IF

      CASE(3)
         WRITE (iout,*) ' +++++ ILUT Preconditioner ++++ '
         IF (lfil < 0) lfil = 30 
         IF (tol  < 0) tol =  0.0001d0
         WRITE (iout,*) ' +++++ tol = ',tol,' lfil = ',lfil,'++++ '
         n_work = (2*lfil+1)*n_syst
         ALLOCATE(wk(n_syst), iw(2*n_syst))

         DO
            ALLOCATE(au(n_work), jau(n_work))
            CALL ilut (n_syst,a,ja,ia,lfil,tol,au,jau,ju,n_work,wk,iw,ierr)
            
            IF(ierr == 0) EXIT

            IF (ierr > 0) THEN
               WRITE(iout,*) ' SOLVE_SPKIT: zero pivot at step ',ierr
               STOP
            ELSE IF(ierr == -1) THEN
               WRITE(iout,*) ' SOLVE_SPKIT: input matrix is wrong '
               STOP
            ELSE IF((ierr == -2) .OR. (ierr.EQ.-3)) THEN
               WRITE(iout,*) ' SOLVE_SPKIT: Insufficient storage for LU factorisation'
               WRITE(iout,*) ' SOLVE_SPKIT: n_work is set to 2*n_work'
               n_work = 2*n_work
               DEALLOCATE(au,jau)
               CYCLE
            ELSE IF(ierr.EQ.-4) THEN
               WRITE(iout,*) ' SOLVE_SPKIT: Illegal value of lfil'
               STOP
            ELSE
               WRITE(iout,*) ' SOLVE_SPKIT: zero row encountered'
               STOP
            ENDIF
         ENDDO

      CASE(4)
         WRITE (iout,*) ' +++++ ILUTP Preconditioner ++++ '
         IF (lfil < 0) lfil = 30 
         IF (tol  < 0) tol =  0.0001d0
         WRITE (iout,*) ' +++++ tol = ',tol,' lfil = ',lfil,'++++ '
         n_work = (2*lfil+1)*n_syst
         ALLOCATE(wk(n_syst), iw(2*n_syst), iperm(2*n_syst))

         DO
            ALLOCATE(au(n_work), jau(n_work))
            CALL ilutp(n_syst,a,ja,ia,lfil,tol,permtol,n_bloc,au,jau,ju,n_work,&
                 wk,iw,iperm,ierr)

            ! TO AVOID PERMUTING THE SOLUTION VECTORS ARRAYS FOR EACH LU-SOLVE,
            ! THE MATRIX A IS PERMUTED ON RETURN. [all column indices are
            ! changed]. SIMILARLY FOR THE U MATRIX.
            ! To permute the matrix back to its original state use the loop:

            DO k=ia(1), ia(n_syst+1)-1
               ja(k) = iperm(ja(k))
            ENDDO

            IF (ierr == 0) EXIT

            IF (ierr > 0) THEN
               WRITE(iout,*) ' SOLVE_SPKIT: zero pivot at step ',ierr
               STOP
            ELSE IF(ierr == -1) THEN
               WRITE(iout,*) ' SOLVE_SPKIT: input matrix is wrong '
               STOP
            ELSE IF((ierr == -2) .OR. (ierr == -3)) THEN
               WRITE(iout,*) ' SOLVE_SPKIT: Insufficient storage for LU factorisation'
               WRITE(iout,*) ' SOLVE_SPKIT: n_work is set to 2*n_work'
               n_work = 2*n_work
               DEALLOCATE(au,jau)
               CYCLE
            ELSE IF(ierr == -4) THEN
               WRITE(iout,*) ' SOLVE_SPKIT: Illegal value of lfil'
               STOP
            ELSE 
               WRITE(iout,*) ' SOLVE_SPKIT: zero row encountered'
               STOP
            END IF
         ENDDO

      CASE(5)
         WRITE (iout,*) ' +++++ ILUK Preconditioner ++++ '
         IF ((lfil < 1) .OR. (lfil > 20)) lfil = 6
         WRITE (iout,*) ' +++++  lfil = ',lfil,'++++ '
         IF (n_work.LE.0) n_work = 60*lfil*n_syst
         ALLOCATE(ju(n_syst), wk(n_syst),iw(3*n_syst))

         DO
            ALLOCATE(au(n_work+n_syst), jau(n_work), levs(n_work))
            CALL iluk(n_syst,a,ja,ia,lfil,au,jau,ju,levs,n_work,wk,iw,ierr)

            IF (ierr == 0) EXIT

            IF (ierr > 0) THEN
               WRITE(iout,*) ' SOLVE_SPKIT: zero pivot at step ',ierr
               STOP
            ELSE IF(ierr == -1) THEN
               WRITE(iout,*) ' SOLVE_SPKIT: input matrix is wrong '
               STOP
            ELSE IF((ierr == -2) .OR. (ierr == -3)) THEN
               WRITE(iout,*) ' SOLVE_SPKIT: Insufficient storage for LU factorisation'
               WRITE(iout,*) ' SOLVE_SPKIT: n_work is set to 2*n_work'
               n_work = 2*n_work
               DEALLOCATE(au,jau,levs)
               CYCLE
            ELSE IF(ierr == -4) THEN
               WRITE(iout,*) ' SOLVE_SPKIT: Illegal value of lfil'
               STOP
            ELSE 
               WRITE(iout,*) ' SOLVE_SPKIT: zero row encountered'
               STOP
            END IF
         ENDDO

      CASE(6)
         WRITE (iout,*) ' +++++ ILUD Preconditioner ++++ '
         IF (tol < 0.d0) tol = 0.01 
         IF ((alpha < 0.d0) .OR. (alpha > 1.d0)) alpha = 0.d0 
         WRITE (iout,*) ' +++++ tol = ',tol,' alpha = ',alpha,'++++ '
         IF (n_work.LE.0) n_work = 100*n_syst
         ALLOCATE(ju(n_syst),wk(n_syst),iw(2*n_syst))

         DO
            ALLOCATE(au(n_work+n_syst),jau(n_work))
            CALL ilud(n_syst,a,ja,ia,alpha,tol,au,jau,ju,n_work,wk,iw,ierr)

            IF (ierr == 0) EXIT

            IF (ierr > 0) THEN
               WRITE(iout,*) ' SOLVE_SPKIT: zero pivot at step ',ierr
               STOP
            ELSE IF(ierr == -1) THEN
               WRITE(iout,*) ' SOLVE_SPKIT: input matrix is wrong '
               STOP
            ELSE IF(ierr == -2) THEN
               WRITE(iout,*) ' SOLVE_SPKIT: Insufficient storage for LU factorisation'
               WRITE(iout,*) ' SOLVE_SPKIT: n_work is set to 2*n_work'
               n_work = 2*n_work
               DEALLOCATE(au,jau)
               CYCLE
            ELSE 
               WRITE(iout,*) ' SOLVE_SPKIT: zero row encountered'
               STOP
            END IF
         ENDDO

      CASE DEFAULT
         WRITE(iout,*)'Preconditioning method', meth, ' not existent!'
         STOP

      END SELECT

      IF(ALLOCATED(iw))    DEALLOCATE(iw)
      IF(ALLOCATED(wk))    DEALLOCATE(wk)
      IF(ALLOCATED(iperm)) DEALLOCATE(iperm)
      IF(ALLOCATED(levs))  DEALLOCATE(levs)

    END SUBROUTINE precond_it
    !
    !----------------------------------------------------------
    !
  END SUBROUTINE solve_spkit
  !
  !********************************************************************
  !
  SUBROUTINE cancel_mem_precond(index)

    INTEGER, INTENT(IN) :: index

    INTEGER :: well_done

    IF(memo(index) % ju_size /= 0) THEN
       DEALLOCATE(memo(index) % ju, STAT = well_done)
       IF(well_done /= 0) THEN
          WRITE(*,*)'Deallocation of stored preconditioner failed! (ju)'
          STOP
       ENDIF
       DEALLOCATE(memo(index) % jau, STAT = well_done)
       IF(well_done /= 0) THEN
          WRITE(*,*)'Deallocation of stored preconditioner failed! (jau)'
          STOP
       ENDIF
       DEALLOCATE(memo(index) % au, STAT = well_done)
       IF(well_done /= 0) THEN
          WRITE(*,*)'Deallocation of stored preconditioner failed! (au)'
          STOP
       ENDIF
    ENDIF

  END SUBROUTINE cancel_mem_precond
  !
  !********************************************************************
  !
END MODULE solve
