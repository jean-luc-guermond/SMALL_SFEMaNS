!
!Authors:  Katarzyna Boronska, Jean-Luc Guermond, Copyrights 2007
!
MODULE sft_parallele

  IMPLICIT NONE
  PUBLIC :: FFT_PAR_CROSS_PROD, FFT_PAR_CROSS_PROD_MOD, ref
  PRIVATE

CONTAINS



  SUBROUTINE FFT_PAR_CROSS_PROD(V1_in, V2_in, V_out, temps, padding)
!This a de-aliased version of the code, FEB 4, 2011, JLG
    IMPLICIT NONE
    INCLUDE 'fftw3.f'
    INCLUDE 'mpif.h'
    ! Format: V_1in(1:np,1:6,1:m_max_c) 
    ! INPUT ARE COSINE AND SINE COEFFICIENTS
    ! THEY ARE PUT IN COMPLEX FORMAT: c_0 = a_0 + i*0 and c_n = (a_n-i*b_n)/2
    REAL(KIND=8), DIMENSION(:,:,:),  INTENT(IN)  :: V1_in, V2_in
    REAL(KIND=8), DIMENSION(:,:,:),  INTENT(OUT) :: V_out
    REAL(KIND=8), DIMENSION(:), OPTIONAL, INTENT(INOUT) :: temps
    LOGICAL,                    OPTIONAL, INTENT(IN)    :: padding
    ! Saved variables
    LOGICAL, SAVE   :: once=.TRUE.
    INTEGER, SAVE   :: np_ref, np, np_tot, bloc_size, nb_field, &
         m_max, m_max_c, N_r, rang, nb_procs, MPID, m_max_pad, N_r_pad
    INTEGER(KIND=8), SAVE :: fftw_plan_multi_c2r, fftw_plan_multi_r2c
    COMPLEX(KIND=8), ALLOCATABLE, DIMENSION(:,:,:), SAVE :: cu, prod_cu 
    REAL(KIND=8),    ALLOCATABLE, DIMENSION(:,:,:), SAVE :: ru, prod_ru
    ! End saved variables

    INTEGER :: i_field 
    INTEGER :: np_glob, np_loc, reste, np_alloc, nn, nb_bloc, n_sup, n_inf, n_up, n_cache
    INTEGER :: nb_valeurs, nb, nf, shiftc, shiftl, index, jindex, longueur_tranche, i, n, p, code
    REAL(KIND=8) :: t
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:,:) :: dist_field, combined_field 
    COMPLEX(KIND=8), ALLOCATABLE, DIMENSION(:,:,:) :: combined_prod_cu, dist_prod_cu, out_prod_cu

    !Vectors to speed up the format changes
    COMPLEX(KIND=8), ALLOCATABLE, DIMENSION(:,:)  :: intermediate
    !COMPLEX, ALLOCATABLE, DIMENSION(:,:)  :: intermediate !There was a bug !15/09/2010
    !Vectors to speed up the format changes

    ! FFTW parameters
    INTEGER   :: fft_dim, howmany, istride, ostride, idist, odist
    INTEGER, DIMENSION(1) :: dim, inembed, onembed
    ! Recall complexes must be rescaled
    ! End FFTW parameters

    !Temps(1) = Temps de communication
    !Temps(2) = Temps de calcul
    !Temps(3) = Temps de changement de format

    EXTERNAL hostnm
    !EXTERNAL gethostname
    character(LEN=64) :: machine

    CALL MPI_COMM_SIZE(MPI_COMM_WORLD,nb_procs,code)
    CALL MPI_COMM_RANK(MPI_COMM_WORLD,rang,code)

    IF (PRESENT(temps)) temps = 0.d0

    IF (.NOT.once) THEN
       IF (SIZE(V1_in,1).NE.np_ref) THEN
          once = .TRUE. !Something wrong happened, reinitialize
          np_ref = SIZE(V1_in,1)
       END IF
    END IF

   once = .TRUE.

! n_cache=decoupage pr MPI_ALL_TO_ALL des points meridiens pr que chaque paquet
! ds la memoire cache-bus qui fait les transferts
   CALL hostnm(machine)
   !CALL gethostname(machine,64)
    IF (machine(1:7)=='compute' .OR. machine(1:6)=='turing') THEN
       n_cache = 5000
   ELSE IF(machine(1:6) =='vargas') THEN
       !n_cache= 10000/nb_procs + 10000
       n_cache=  10000
   ELSE IF (machine(1:6)=='brazos' .OR. machine(1:2)=='c0') THEN
       IF (nb_procs>16) THEN 
          n_cache= 600000
       ELSE IF (nb_procs>8) THEN
          n_cache= 150000 ! Tuned for Brazos, redo for Vargas
       ELSE
          n_cache = 5000 !20000
       END IF
    ELSE
       WRITE(*,*) ' BUG: Machine cache not programmed yet in fft_paralell'
       WRITE(*,*) machine
       STOP
    END IF
    np_glob      = SIZE(V1_in,1)
    m_max_c = SIZE(V1_in,3)
    np_loc = n_cache/(12*m_max_c)
    nb_bloc = MAX(np_glob/np_loc,1)

100 np_loc = np_glob/nb_bloc
    reste = np_glob - np_loc*nb_bloc
    np_alloc = np_loc + reste
    reste = np_alloc*nb_bloc - np_glob
    IF (reste>np_alloc) THEN
       nb_bloc = nb_bloc - 1
       GO TO 100
    END IF
    ! nb_bloc = nbre de blocs qui decoupe le plan meridien
    ! np_alloc = nbre de points ds 1 bloc

    n_sup = 0
    DO nn= 1, nb_bloc
       IF (once) THEN
          nb_field= SIZE(V1_in,2) ! Number of fields
          m_max_c = SIZE(V1_in,3) ! Number of complex (cosines + sines) coefficients per point
          m_max = m_max_c*nb_procs! Number of comlex coefficients per point per processor
          N_r=2*m_max-1           ! Number of Real coefficients per point
          IF (MOD(nb_field,2)/=0 .OR. m_max_c==0) THEN
             WRITE(*,*) ' BUG '
             STOP
          END IF

          ! Bloc_size is the number of points that are handled by one processor 
          ! once the Fourier modes are all collected
          ! Computation of bloc_size and np_tot
          np = np_alloc
          IF (MODULO(np,nb_procs)==0) THEN
             bloc_size = np/nb_procs
          ELSE
             bloc_size = np/nb_procs + 1
          END IF
          np_tot = nb_procs*bloc_size
          ! fin de la repartition des points


          !JLG, FEB 4, 2011
          !Only ru, cu, prod_ru, prod_cu are modified
          !IF (padding) THEN
          IF (.true.) THEN
             m_max_pad = 3*m_max/2
          ELSE
             m_max_pad = m_max
          END IF
          N_r_pad=2*m_max_pad-1

          IF (ALLOCATED(ru)) DEALLOCATE(ru,cu,prod_ru,prod_cu)
          ALLOCATE(cu(m_max_pad,nb_field,bloc_size))
          ALLOCATE(ru(N_r_pad,  nb_field,bloc_size))
          ALLOCATE(prod_cu(m_max_pad,nb_field/2,bloc_size))
          ALLOCATE(prod_ru(N_r_pad,  nb_field/2,bloc_size))
          !JLG, FEB 4, 2001
          ALLOCATE(intermediate(nb_field/2,bloc_size))
          ALLOCATE(    dist_field(m_max_c,2*nb_field,np_tot))
          ALLOCATE(combined_field(m_max_c,2*nb_field,np_tot))
          ALLOCATE(dist_prod_cu(nb_field/2,bloc_size,m_max))
          ALLOCATE(combined_prod_cu(nb_field/2,bloc_size,m_max))
          ALLOCATE(out_prod_cu(m_max_c,np_tot,nb_field/2))
       END IF


       ! Packing all 3 complex components of both v1 and v2 input fields
       ! into dist_field, where the dimension indexing the nodal points varies the least rapidly,
       ! so that after distributing the data to the processes, each one will obtain a part 
       ! on nodal points
       ! TRANSPOSE pr que la variable i associee aux modes soit la 1ere sur laquelle on va faire la FFT
       t = MPI_WTIME()
       n_inf = n_sup + 1
       IF (nn == nb_bloc) THEN
          n_up  = np_glob - n_inf + 1
       ELSE
          n_up  = np_alloc
       END IF
       n_sup = n_inf + n_up - 1

       DO i = 1, m_max_c
          dist_field(i,1:nb_field,1:n_up) = TRANSPOSE(V1_in(n_inf:n_sup,:,i))
          dist_field(i,nb_field+1:2*nb_field,1:n_up) = TRANSPOSE(V2_in(n_inf:n_sup,:,i))
       END DO
       IF (PRESENT(temps)) temps(3) = temps(3) + MPI_WTIME() -t

       IF (np/=np_tot) dist_field(:,:,np+1:np_tot) = 1.d100

       longueur_tranche=bloc_size*m_max_c*nb_field*2

       t = MPI_WTIME()
       MPID=MPI_DOUBLE_PRECISION
       CALL MPI_ALLTOALL (dist_field,longueur_tranche, MPID, combined_field, longueur_tranche, &
            MPID, MPI_COMM_WORLD, code)
       IF (PRESENT(temps)) temps(1) = temps(1) + MPI_WTIME() -t

       t = MPI_WTIME()
       !JLG, FEB 4, 2011
       cu = 0.d0
       !JLG, FEB 4, 2011
       DO n = 1, bloc_size 
          DO nb = 1, nb_procs
             shiftc = (nb-1)*bloc_size  
             shiftl = (nb-1)*m_max_c
             jindex = n + shiftc
             DO nf = 1, nb_field
                ! Put real and imaginary parts in a complex 
                ! nf=1,2,3 => V1_in
                ! nf=4,5,6 => V2_in
                ! INPUT ARE COSINE AND SINE COEFFICIENTS
                ! THEY ARE PUT IN COMPLEX FORMAT: c_0 = a_0 + i*0 and c_n = (a_n-i*b_n)/2
                cu(shiftl+1:shiftl+m_max_c,nf,n) = CMPLX(combined_field(:,2*nf-1,jindex),&
                                                        -combined_field(:,2*nf,jindex),KIND=8)/2
             END DO
          END DO
       END DO
       cu(1,:,:) = 2*CMPLX(REAL(cu(1,:,:),KIND=8),0.d0,KIND=8)
       !JLG, FEB 4, 2011
      !Padding is done by initialization of cu: cu = 0
      !This is eequivalent to cu(m_max+1:m_max_pad,:,:) = 0.d0
       !JLG, FEB 4, 2011

       IF (PRESENT(temps)) temps(3) = temps(3) + MPI_WTIME() -t

       ! Set the parameters for dfftw
       fft_dim=1; istride=1; ostride=1; 
       !JLG, FEB 4, 2011
!!$       idist=N_r;   inembed(1)=N_r; DIM(1)=N_r
!!$       odist=m_max; onembed(1)=m_max
       idist=N_r_pad;   inembed(1)=N_r_pad; DIM(1)=N_r_pad
       odist=m_max_pad; onembed(1)=m_max_pad
       !JLG, FEB 4, 2011

       howmany=bloc_size*nb_field
 

       t = MPI_WTIME()
       IF (once) CALL dfftw_plan_many_dft_c2r(fftw_plan_multi_c2r, fft_dim, dim, howmany, cu, &
            onembed, ostride, odist, ru, inembed, istride, idist, FFTW_ESTIMATE);
       CALL dfftw_execute(fftw_plan_multi_c2r)

       ! CROSS PRODDUCT
       IF (nb_field==6) THEN
          prod_ru(:,1,:) = ru(:,2,:)*ru(:,6,:) - ru(:,3,:)*ru(:,5,:)
          prod_ru(:,2,:) = ru(:,3,:)*ru(:,4,:) - ru(:,1,:)*ru(:,6,:)
          prod_ru(:,3,:) = ru(:,1,:)*ru(:,5,:) - ru(:,2,:)*ru(:,4,:)
       END IF
       ! CROSS PRODUCT

       howmany = howmany/2
       IF (once) CALL dfftw_plan_many_dft_r2c(fftw_plan_multi_r2c, fft_dim, dim, howmany, prod_ru, &
            inembed, istride, idist, prod_cu, onembed, ostride, odist, FFTW_ESTIMATE);
       CALL dfftw_execute(fftw_plan_multi_r2c)
       !JLG, FEB 4, 2011
!!$       prod_cu = prod_cu/N_r !Scaling 
       prod_cu = prod_cu/N_r_pad !Scaling 
       !JLG, FEB 4, 2011
       IF (PRESENT(temps)) temps(2) = temps(2) + MPI_WTIME() -t

       !Now we need to redistribute the Fourier coefficients on each processor

       t = MPI_WTIME()
       combined_prod_cu(:,:,1)=prod_cu(1,:,:)
       DO n=2, m_max
          !combined_prod_cu(:,:,n)=prod_cu(n,:,:)
          combined_prod_cu(:,:,n)=2*CONJG(prod_cu(n,:,:))
       END DO

       IF (PRESENT(temps)) temps(3) = temps(3) + MPI_WTIME() -t

       t = MPI_WTIME()
       longueur_tranche=bloc_size*m_max_c*nb_field
       MPID=MPI_DOUBLE_PRECISION
       CALL MPI_ALLTOALL (combined_prod_cu,longueur_tranche,MPID, dist_prod_cu,longueur_tranche, &
            MPID,MPI_COMM_WORLD,code)
       IF (PRESENT(temps)) temps(1) = temps(1) + MPI_WTIME() -t
       ! dimensions: 
       t = MPI_WTIME()
       DO i = 1, m_max_c
          DO nb = 1, nb_procs
             shiftc = (nb-1)*bloc_size  
             shiftl = (nb-1)*m_max_c
             intermediate = dist_prod_cu(:,:,shiftl+i)
             DO n = 1, bloc_size
                IF (n_inf-1+n+shiftc > np_glob ) CYCLE 
                DO i_field = 1, nb_field/2
                   v_out(n_inf-1+n+shiftc, i_field*2-1, i) = REAL (intermediate(i_field,n),KIND=8)
                   v_out(n_inf-1+n+shiftc, i_field*2 , i)  = AIMAG(intermediate(i_field,n))
                END DO
             END DO
          END DO
       END DO
       IF (PRESENT(temps)) temps(3) = temps(3) + MPI_WTIME() -t

       once = .FALSE.

    END DO

    DEALLOCATE(dist_field, combined_field, dist_prod_cu, combined_prod_cu, intermediate, out_prod_cu)

  END SUBROUTINE FFT_PAR_CROSS_PROD

  SUBROUTINE ref(V1_in, V2_in, V_out, temps)
    IMPLICIT NONE

    INCLUDE 'fftw3.f'
    INCLUDE 'mpif.h'
    ! Format: V_1in(1:np,1:6,1:m_max_c) 
    ! INPUT ARE COSINE AND SINE COEFFICIENTS
    ! THEY ARE PUT IN COMPLEX FORMAT: c_0 = a_0 + i*0 and c_n = (a_n-i*b_n)/2
    REAL(KIND=8), DIMENSION(:,:,:),  INTENT(IN)  :: V1_in, V2_in
    REAL(KIND=8), DIMENSION(:,:,:),  INTENT(OUT) :: V_out
    REAL(KIND=8), DIMENSION(:), OPTIONAL, INTENT(INOUT) :: temps

    ! Saved variables
    LOGICAL, SAVE   :: once=.TRUE.
    INTEGER, SAVE   :: np, np_tot, bloc_size, nb_field, m_max, m_max_c, N_r, rang, nb_procs, MPID
    INTEGER(KIND=8), SAVE :: fftw_plan_multi_c2r, fftw_plan_multi_r2c
    COMPLEX(KIND=8), ALLOCATABLE, DIMENSION(:,:,:), SAVE :: cu, prod_cu 
    REAL(KIND=8),    ALLOCATABLE, DIMENSION(:,:,:), SAVE :: ru, prod_ru
    ! End saved variables

    INTEGER :: i_field
    INTEGER :: nb_valeurs, nb, nf, shiftc, shiftl, index, jindex, longueur_tranche, i, n, p, code
    REAL(KIND=8) :: t
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:,:), SAVE :: dist_field, combined_field 
    COMPLEX(KIND=8), ALLOCATABLE, DIMENSION(:,:,:), SAVE :: combined_prod_cu, dist_prod_cu, out_prod_cu

    !Vectors to speed up the format changes
    COMPLEX(KIND=8), ALLOCATABLE, DIMENSION(:,:)  :: intermediate
    !REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:) :: V_int
    !Vectors to speed up the format changes

    ! FFTW parameters
    INTEGER   :: fft_dim, howmany, istride, ostride, idist, odist
    INTEGER, DIMENSION(1) :: dim, inembed, onembed
    ! Recall complexes must be rescaled
    ! End FFTW parameters

    !Temps(1) = Temps de communication
    !Temps(2) = Temps de calcul
    !Temps(3) = Temps de changement de format
    IF (PRESENT(temps)) temps = 0.d0

    IF (.NOT.once) THEN
       IF ((SIZE(V1_in,1).NE.np) .OR. (SIZE(V1_in,2).NE.nb_field) .OR. (SIZE(V1_in,3).NE.m_max_c)) THEN
          once = .TRUE. !Something wrong happened, reinitialize
       END IF
    END IF

    IF (once) THEN
       CALL MPI_COMM_SIZE(MPI_COMM_WORLD,nb_procs,code)
       CALL MPI_COMM_RANK(MPI_COMM_WORLD,rang,code)

       np      = SIZE(V1_in,1) ! Number of points in the meridian plane
       nb_field= SIZE(V1_in,2) ! Number of fields
       m_max_c = SIZE(V1_in,3) ! Number of complex (cosines + sines) coefficients per point
       m_max = m_max_c*nb_procs! Number of comlex coefficients per point per processor
       N_r=2*m_max-1           ! Number of Real coefficients per point
       IF (MOD(nb_field,2)/=0 .OR. m_max_c==0) THEN
          WRITE(*,*) ' BUG '
          STOP
       END IF

       ! Bloc_size is the number of points that are handled by one processor 
       ! once the Fourier modes are all collected
       ! Computation of bloc_size and np_tot
       IF (MODULO(np,nb_procs)==0) THEN
          bloc_size = np/nb_procs
       ELSE
          bloc_size = np/nb_procs + 1
       END IF
       np_tot = nb_procs*bloc_size

       IF (ALLOCATED(ru)) DEALLOCATE(ru,cu,prod_ru,prod_cu)
       ALLOCATE(cu(m_max,nb_field,bloc_size))
       ALLOCATE(ru(N_r,  nb_field,bloc_size))
       ALLOCATE(prod_cu(m_max,nb_field/2,bloc_size))
       ALLOCATE(prod_ru(N_r,  nb_field/2,bloc_size))
       !ALLOCATE(V_int(nb_field,np))
       !ALLOCATE(V_int(np,nb_field))
    END IF

    ALLOCATE(intermediate(nb_field/2,bloc_size))
    ALLOCATE(    dist_field(m_max_c,2*nb_field,np_tot))
    ALLOCATE(combined_field(m_max_c,2*nb_field,np_tot))
    ALLOCATE(dist_prod_cu(nb_field/2,bloc_size,m_max), combined_prod_cu(nb_field/2,bloc_size,m_max))
    ALLOCATE(out_prod_cu(m_max_c,np_tot,nb_field/2))

    ! Packing all 3 complex components of both v1 and v2 input fields
    ! into dist_field, where the dimension indexing the nodal points varies the least rapidly,
    ! so that after distributing the data to the processes, each one will obtain a part 
    ! on nodal points
    t = MPI_WTIME()
    DO i = 1, m_max_c
       !V_int = V1_in(:,:,i)
       !dist_field(i,1:nb_field,:) = transpose(V_int)
       !V_int = V2_in(:,:,i)
       !dist_field(i,nb_field+1:2*nb_field,:) = transpose(V_int)
       dist_field(i,1:nb_field,1:np) = TRANSPOSE(V1_in(:,:,i))
       dist_field(i,nb_field+1:2*nb_field,1:np) = TRANSPOSE(V2_in(:,:,i))
       !DO nf = 1, nb_field/2 
       !   dist_field(i,2*nf-1,1:np) = V1_in(:,2*nf-1,i) 
       !   dist_field(i,2*nf  ,1:np) = V1_in(:,2*nf  ,i) 
       !   dist_field(i,nb_field+2*nf-1,1:np) = V2_in(:,2*nf-1,i)
       !   dist_field(i,nb_field+2*nf  ,1:np) = V2_in(:,2*nf,i)
       !END DO
    END DO
    IF (PRESENT(temps)) temps(3) = temps(3) + MPI_WTIME() -t

    IF (np/=np_tot) dist_field(:,:,np+1:np_tot) = 1.d100

    longueur_tranche=bloc_size*m_max_c*nb_field*2

    t = MPI_WTIME()
    MPID=MPI_DOUBLE_PRECISION
    CALL MPI_ALLTOALL (dist_field,longueur_tranche, MPID, combined_field, longueur_tranche, &
         MPID, MPI_COMM_WORLD, code)
    IF (PRESENT(temps)) temps(1) = temps(1) + MPI_WTIME() -t

    t = MPI_WTIME()

    DO n = 1, bloc_size 
       DO nb = 1, nb_procs
          shiftc = (nb-1)*bloc_size  
          shiftl = (nb-1)*m_max_c
          jindex = n + shiftc
          DO nf = 1, nb_field
             ! Put real and imaginary parts in a complex 
             ! for each field, nf=1,2,3,4,5,6
             ! nf=1,2,3 => V1_in
             ! nf=4,5,6 => V2_in
             ! cu(shiftl+1:shiftl+m_max_c,nf,n) = cmplx(combined_field(:,2*nf-1,jindex),combined_field(:,2*nf,jindex))
             ! INPUT ARE COSINE AND SINE COEFFICIENTS
             ! THEY ARE PUT IN COMPLEX FORMAT: c_0 = a_0 + i*0 and c_n = (a_n-i*b_n)/2
             cu(shiftl+1:shiftl+m_max_c,nf,n) = 0.5d0*CMPLX(combined_field(:,2*nf-1,jindex),-combined_field(:,2*nf,jindex),KIND=8)
          END DO
       END DO
    END DO
    cu(1,:,:) = 2.d0*CMPLX(REAL(cu(1,:,:),KIND=8),0.d0,KIND=8)
    IF (PRESENT(temps)) temps(3) = temps(3) + MPI_WTIME() -t

    ! Set the parameters for dfftw
    fft_dim=1; istride=1; ostride=1; 
    idist=N_r;   inembed(1)=N_r; DIM(1)=N_r
    odist=m_max; onembed(1)=m_max
    IF (rang==(nb_procs-1)) THEN
       howmany= (np - bloc_size*(nb_procs-1))*nb_field 
    ELSE 
       howmany=bloc_size*nb_field
    END IF

    t = MPI_WTIME()
    IF (once) CALL dfftw_plan_many_dft_c2r(fftw_plan_multi_c2r, fft_dim, dim, howmany, cu, &
         onembed, ostride, odist, ru, inembed, istride, idist, FFTW_ESTIMATE);
    CALL dfftw_execute(fftw_plan_multi_c2r)

    ! SIMPLE PRODUCT
    !DO nf = 1, nb_field/2
    !   ! ru(1:3) contains V1_in and ru(4:6) contains V2_in
    !   prod_ru(:,nf,:)  = ru(:,nf,:)*ru(:,nb_field/2+nf,:)
    !END DO
    ! END SIMPLE PRODUCT

    ! CROSS PRODDUCT
    IF (nb_field==6) THEN
       prod_ru(:,1,:) = ru(:,2,:)*ru(:,6,:) - ru(:,3,:)*ru(:,5,:)
       prod_ru(:,2,:) = ru(:,3,:)*ru(:,4,:) - ru(:,1,:)*ru(:,6,:)
       prod_ru(:,3,:) = ru(:,1,:)*ru(:,5,:) - ru(:,2,:)*ru(:,4,:)
    END IF
    ! CROSS PRODUCT

    howmany = howmany/2
    IF (once) CALL dfftw_plan_many_dft_r2c(fftw_plan_multi_r2c, fft_dim, dim, howmany, prod_ru, &
         inembed, istride, idist, prod_cu, onembed, ostride, odist, FFTW_ESTIMATE);
    CALL dfftw_execute(fftw_plan_multi_r2c)

    prod_cu = prod_cu/N_r !Scaling 
    !   prod_cu = prod_cu/REAL(N_r,KIND=8) !Scaling 
    IF (PRESENT(temps)) temps(2) = temps(2) + MPI_WTIME() -t

    !Now we need to redistribute the Fourier coefficients on each processor

    t = MPI_WTIME()
    combined_prod_cu(:,:,1)=prod_cu(1,:,:)
    DO n=2, m_max
       !combined_prod_cu(:,:,n)=prod_cu(n,:,:)
       combined_prod_cu(:,:,n)=2.d0*CONJG(prod_cu(n,:,:))
    END DO
    IF (PRESENT(temps)) temps(3) = temps(3) + MPI_WTIME() -t

    t = MPI_WTIME()
    longueur_tranche=bloc_size*m_max_c*nb_field
    MPID=MPI_DOUBLE_PRECISION
    CALL MPI_ALLTOALL (combined_prod_cu,longueur_tranche,MPID, dist_prod_cu,longueur_tranche, &
         MPID,MPI_COMM_WORLD,code)
    IF (PRESENT(temps)) temps(1) = temps(1) + MPI_WTIME() -t

    ! dimensions: 
    ! v_out(np, nb_field, m_max_c)
    t = MPI_WTIME()

    IF (.FALSE.) THEN
       DO i_field = 1, nb_field/2
          DO n = 1, bloc_size
             DO nb = 1, nb_procs
                shiftc = (nb-1)*bloc_size
                shiftl = (nb-1)*m_max_c
                out_prod_cu(:,n+shiftc,i_field) = dist_prod_cu(i_field,n,shiftl+1:shiftl+m_max_c)
             END DO
          END DO
       END DO

       DO i_field = 1, nb_field/2
          DO i = 1, m_max_c
             v_out(:, i_field*2-1, i) = REAL(out_prod_cu(i, 1:np, i_field),KIND=8)
             v_out(:, i_field*2,   i) = AIMAG(out_prod_cu(i, 1:np, i_field))
          END DO
       END DO
    ELSE
       DO i = 1, m_max_c
          DO nb = 1, nb_procs
             shiftc = (nb-1)*bloc_size  
             shiftl = (nb-1)*m_max_c
             intermediate = dist_prod_cu(:,:,shiftl+i)
             DO n = 1, bloc_size
                IF (n+shiftc > np ) CYCLE 
                DO i_field = 1, nb_field/2
                   v_out(n+shiftc, i_field*2-1, i) = REAL (intermediate(i_field,n),KIND=8)
                   v_out(n+shiftc, i_field*2 , i)  = AIMAG(intermediate(i_field,n))
                   !v_out(n+shiftc, i_field*2-1, i) = REAL(dist_prod_cu(i_field,n,shiftl+i),KIND=8)
                   !v_out(n+shiftc, i_field*2 , i)  = AIMAG(dist_prod_cu(i_field,n,shiftl+i))
                END DO
             END DO
          END DO
       END DO
    END IF

    IF (PRESENT(temps)) temps(3) = temps(3) + MPI_WTIME() -t

    DEALLOCATE(dist_field, combined_field, dist_prod_cu, combined_prod_cu, intermediate, out_prod_cu) 
    IF (once) once = .FALSE.
  END SUBROUTINE ref


  SUBROUTINE FFT_PAR_CROSS_PROD_copy(V1_in, V2_in, V_out, temps)
    IMPLICIT NONE

    INCLUDE 'fftw3.f'
    INCLUDE 'mpif.h'

    ! Format: V_1in(1:np,1:6,1:m_max_c) 
    ! INPUT ARE COSINE AND SINE COEFFICIENTS
    ! THEY ARE PUT IN COMPLEX FORMAT: c_0 = a_0 + i*0 and c_n = (a_n-i*b_n)/2
    REAL(KIND=8), DIMENSION(:,:,:),  INTENT(IN)  :: V1_in, V2_in
    REAL(KIND=8), DIMENSION(:,:,:),  INTENT(OUT) :: V_out
    REAL(KIND=8), DIMENSION(:), OPTIONAL, INTENT(INOUT) :: temps

    ! Saved variables
    LOGICAL, SAVE   :: once=.TRUE.
    INTEGER, SAVE   :: np_ref, np, np_tot, bloc_size, nb_field, m_max, m_max_c, N_r, rang, nb_procs, MPID
    INTEGER(KIND=8), SAVE :: fftw_plan_multi_c2r, fftw_plan_multi_r2c
    COMPLEX(KIND=8), ALLOCATABLE, DIMENSION(:,:,:), SAVE :: cu, prod_cu 
    REAL(KIND=8),    ALLOCATABLE, DIMENSION(:,:,:), SAVE :: ru, prod_ru
    ! End saved variables

    INTEGER :: i_field 
    INTEGER :: np_glob, np_loc, reste, np_alloc, nn, nb_bloc, n_sup, n_inf, n_up, n_cache
    INTEGER :: nb_valeurs, nb, nf, shiftc, shiftl, index, jindex, longueur_tranche, i, n, p, code
    REAL(KIND=8) :: t
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:,:) :: dist_field, combined_field 
    COMPLEX(KIND=8), ALLOCATABLE, DIMENSION(:,:,:) :: combined_prod_cu, dist_prod_cu, out_prod_cu

    !Vectors to speed up the format changes
    COMPLEX(KIND=8), ALLOCATABLE, DIMENSION(:,:)  :: intermediate
    !COMPLEX, ALLOCATABLE, DIMENSION(:,:)  :: intermediate !There was a bug !15/09/2010
    !Vectors to speed up the format changes

    ! FFTW parameters
    INTEGER   :: fft_dim, howmany, istride, ostride, idist, odist
    INTEGER, DIMENSION(1) :: dim, inembed, onembed
    ! Recall complexes must be rescaled
    ! End FFTW parameters

    !Temps(1) = Temps de communication
    !Temps(2) = Temps de calcul
    !Temps(3) = Temps de changement de format

    EXTERNAL hostnm
    !EXTERNAL gethostname
    character(LEN=64) :: machine

    CALL MPI_COMM_SIZE(MPI_COMM_WORLD,nb_procs,code)
    CALL MPI_COMM_RANK(MPI_COMM_WORLD,rang,code)

    IF (PRESENT(temps)) temps = 0.d0

    IF (.NOT.once) THEN
       IF (SIZE(V1_in,1).NE.np_ref) THEN
          once = .TRUE. !Something wrong happened, reinitialize
          np_ref = SIZE(V1_in,1)
       END IF
    END IF

   once = .TRUE.

! n_cache=decoupage pr MPI_ALL_TO_ALL des points meridiens pr que chaque paquet
! ds la memoire cache-bus qui fait les transferts
   CALL hostnm(machine)
   !CALL gethostname(machine,64)
    IF (machine(1:7)=='compute' .OR. machine(1:6)=='turing') THEN
       n_cache = 5000
   ELSE IF(machine(1:6) =='vargas') THEN
       !n_cache= 10000/nb_procs + 10000
       n_cache=  10000
   ELSE IF (machine(1:6)=='brazos' .OR. machine(1:2)=='c0') THEN
       IF (nb_procs>16) THEN 
          n_cache= 600000
       ELSE IF (nb_procs>8) THEN
          n_cache= 150000 ! Tuned for Brazos, redo for Vargas
       ELSE
          n_cache = 5000 !20000
       END IF
    ELSE
       WRITE(*,*) ' BUG: Machine cache not programmed yet in fft_paralell'
       WRITE(*,*) machine
       STOP
    END IF
    np_glob      = SIZE(V1_in,1)
    m_max_c = SIZE(V1_in,3)
    np_loc = n_cache/(12*m_max_c)
    nb_bloc = MAX(np_glob/np_loc,1)

100 np_loc = np_glob/nb_bloc
    reste = np_glob - np_loc*nb_bloc
    np_alloc = np_loc + reste
    reste = np_alloc*nb_bloc - np_glob
    IF (reste>np_alloc) THEN
       nb_bloc = nb_bloc - 1
       GO TO 100
    END IF
! nb_bloc = nbre de blocs qui decoupe le plan meridien
! np_alloc = nbre de points ds 1 bloc
    !TEST

    !nb_bloc = 1
    !np_alloc =np_glob

    n_sup = 0
    DO nn= 1, nb_bloc
       IF (once) THEN
          nb_field= SIZE(V1_in,2) ! Number of fields
          m_max_c = SIZE(V1_in,3) ! Number of complex (cosines + sines) coefficients per point
          m_max = m_max_c*nb_procs! Number of comlex coefficients per point per processor
          N_r=2*m_max-1           ! Number of Real coefficients per point
          IF (MOD(nb_field,2)/=0 .OR. m_max_c==0) THEN
             WRITE(*,*) ' BUG '
             STOP
          END IF

          ! Bloc_size is the number of points that are handled by one processor 
          ! once the Fourier modes are all collected
          ! Computation of bloc_size and np_tot
          np = np_alloc
          IF (MODULO(np,nb_procs)==0) THEN
             bloc_size = np/nb_procs
          ELSE
             bloc_size = np/nb_procs + 1
          END IF
          np_tot = nb_procs*bloc_size
! bloc_size = nbre de points par processeur
! np_tot = nbre de points ds 1 bloc -> remplace le np_alloc
! fin de la repartition des points


          IF (ALLOCATED(ru)) DEALLOCATE(ru,cu,prod_ru,prod_cu)
          ALLOCATE(cu(m_max,nb_field,bloc_size))
          ALLOCATE(ru(N_r,  nb_field,bloc_size))
          ALLOCATE(prod_cu(m_max,nb_field/2,bloc_size))
          ALLOCATE(prod_ru(N_r,  nb_field/2,bloc_size))
          ALLOCATE(intermediate(nb_field/2,bloc_size))
          ALLOCATE(    dist_field(m_max_c,2*nb_field,np_tot))
          ALLOCATE(combined_field(m_max_c,2*nb_field,np_tot))
          ALLOCATE(dist_prod_cu(nb_field/2,bloc_size,m_max))
          ALLOCATE(combined_prod_cu(nb_field/2,bloc_size,m_max))
          ALLOCATE(out_prod_cu(m_max_c,np_tot,nb_field/2))
       END IF


       ! Packing all 3 complex components of both v1 and v2 input fields
       ! into dist_field, where the dimension indexing the nodal points varies the least rapidly,
       ! so that after distributing the data to the processes, each one will obtain a part 
       ! on nodal points
! TRANSPOSE pr que la variable i associee aux modes soit la 1ere sur laquelle on va faire la FFT
       t = MPI_WTIME()
       n_inf = n_sup + 1
       IF (nn == nb_bloc) THEN
          n_up  = np_glob - n_inf + 1
       ELSE
          n_up  = np_alloc
       END IF
       n_sup = n_inf + n_up - 1

       DO i = 1, m_max_c
          dist_field(i,1:nb_field,1:n_up) = TRANSPOSE(V1_in(n_inf:n_sup,:,i))
          dist_field(i,nb_field+1:2*nb_field,1:n_up) = TRANSPOSE(V2_in(n_inf:n_sup,:,i))
       END DO
       IF (PRESENT(temps)) temps(3) = temps(3) + MPI_WTIME() -t

       IF (np/=np_tot) dist_field(:,:,np+1:np_tot) = 1.d100

       longueur_tranche=bloc_size*m_max_c*nb_field*2

       t = MPI_WTIME()
       MPID=MPI_DOUBLE_PRECISION
       CALL MPI_ALLTOALL (dist_field,longueur_tranche, MPID, combined_field, longueur_tranche, &
            MPID, MPI_COMM_WORLD, code)
       IF (PRESENT(temps)) temps(1) = temps(1) + MPI_WTIME() -t

       t = MPI_WTIME()
       DO n = 1, bloc_size 
          DO nb = 1, nb_procs
             shiftc = (nb-1)*bloc_size  
             shiftl = (nb-1)*m_max_c
             jindex = n + shiftc
             DO nf = 1, nb_field
                ! Put real and imaginary parts in a complex 
                ! for each field, nf=1,2,3,4,5,6
                ! nf=1,2,3 => V1_in
                ! nf=4,5,6 => V2_in
                ! INPUT ARE COSINE AND SINE COEFFICIENTS
                ! THEY ARE PUT IN COMPLEX FORMAT: c_0 = a_0 + i*0 and c_n = (a_n-i*b_n)/2
                cu(shiftl+1:shiftl+m_max_c,nf,n) = CMPLX(combined_field(:,2*nf-1,jindex),&
                                                        -combined_field(:,2*nf,jindex),KIND=8)/2
             END DO
          END DO
       END DO
       cu(1,:,:) = 2*CMPLX(REAL(cu(1,:,:),KIND=8),0.d0,KIND=8)
! agir ici pr faire du padding cu_pad = cu + des zeros
! attention a la dimension de cu_pad!

       IF (PRESENT(temps)) temps(3) = temps(3) + MPI_WTIME() -t

       ! Set the parameters for dfftw
       fft_dim=1; istride=1; ostride=1; 
       idist=N_r;   inembed(1)=N_r; DIM(1)=N_r
       odist=m_max; onembed(1)=m_max
!       IF (rang==(nb_procs-1)) THEN
!          !howmany= (np - bloc_size*(nb_procs-1))*nb_field  !JLFG, There was a bug here Sept, 16, 2010
!          howmany= (np_tot - bloc_size*(nb_procs-1))*nb_field 
!       ELSE 
!          howmany=bloc_size*nb_field
!       END IF
       !JLG: Sept, 16, 2010
       howmany=bloc_size*nb_field
       !JLG: Sept, 16, 2010

       t = MPI_WTIME()
       IF (once) CALL dfftw_plan_many_dft_c2r(fftw_plan_multi_c2r, fft_dim, dim, howmany, cu, &
            onembed, ostride, odist, ru, inembed, istride, idist, FFTW_ESTIMATE);
       CALL dfftw_execute(fftw_plan_multi_c2r)

       ! CROSS PRODDUCT
       IF (nb_field==6) THEN
          prod_ru(:,1,:) = ru(:,2,:)*ru(:,6,:) - ru(:,3,:)*ru(:,5,:)
          prod_ru(:,2,:) = ru(:,3,:)*ru(:,4,:) - ru(:,1,:)*ru(:,6,:)
          prod_ru(:,3,:) = ru(:,1,:)*ru(:,5,:) - ru(:,2,:)*ru(:,4,:)
       END IF
       ! CROSS PRODUCT

       howmany = howmany/2
       IF (once) CALL dfftw_plan_many_dft_r2c(fftw_plan_multi_r2c, fft_dim, dim, howmany, prod_ru, &
            inembed, istride, idist, prod_cu, onembed, ostride, odist, FFTW_ESTIMATE);
       CALL dfftw_execute(fftw_plan_multi_r2c)

       prod_cu = prod_cu/N_r !Scaling 
       !   prod_cu = prod_cu/REAL(N_r,KIND=8) !Scaling 
       IF (PRESENT(temps)) temps(2) = temps(2) + MPI_WTIME() -t

       !Now we need to redistribute the Fourier coefficients on each processor

       t = MPI_WTIME()
       combined_prod_cu(:,:,1)=prod_cu(1,:,:)
       DO n=2, m_max
          !combined_prod_cu(:,:,n)=prod_cu(n,:,:)
          combined_prod_cu(:,:,n)=2*CONJG(prod_cu(n,:,:))
       END DO
! agir ici pr faire du de-padding combined_prod_cu_pad = combined_prod_cu + des zeros
! recuperer ensuite combined_prod_cu
!
       IF (PRESENT(temps)) temps(3) = temps(3) + MPI_WTIME() -t

       t = MPI_WTIME()
       longueur_tranche=bloc_size*m_max_c*nb_field
       MPID=MPI_DOUBLE_PRECISION
       CALL MPI_ALLTOALL (combined_prod_cu,longueur_tranche,MPID, dist_prod_cu,longueur_tranche, &
            MPID,MPI_COMM_WORLD,code)
       IF (PRESENT(temps)) temps(1) = temps(1) + MPI_WTIME() -t
       ! dimensions: 
       ! v_out(np, nb_field, m_max_c)
       t = MPI_WTIME()
       DO i = 1, m_max_c
          DO nb = 1, nb_procs
             shiftc = (nb-1)*bloc_size  
             shiftl = (nb-1)*m_max_c
             intermediate = dist_prod_cu(:,:,shiftl+i)
             DO n = 1, bloc_size
                IF (n_inf-1+n+shiftc > np_glob ) CYCLE 
                DO i_field = 1, nb_field/2
                   v_out(n_inf-1+n+shiftc, i_field*2-1, i) = REAL (intermediate(i_field,n),KIND=8)
                   v_out(n_inf-1+n+shiftc, i_field*2 , i)  = AIMAG(intermediate(i_field,n))
                   !v_out(n+shiftc, i_field*2-1, i) = REAL(dist_prod_cu(i_field,n,shiftl+i),KIND=8)
                   !v_out(n+shiftc, i_field*2 , i)  = AIMAG(dist_prod_cu(i_field,n,shiftl+i))
                END DO
             END DO
          END DO
       END DO
       IF (PRESENT(temps)) temps(3) = temps(3) + MPI_WTIME() -t

       once = .FALSE.

    END DO

    DEALLOCATE(dist_field, combined_field, dist_prod_cu, combined_prod_cu, intermediate, out_prod_cu)

  END SUBROUTINE FFT_PAR_CROSS_PROD_copy



  SUBROUTINE fft_par_cross_prod_mod(V1_in, V2_in, V_out, pb_type, temps)
    IMPLICIT NONE
    TYPE dyn_FFTW
       INTEGER :: np, np_tot, bloc_size, nb_field, m_max, m_max_c, N_r
       INTEGER(KIND=8) :: fftw_plan_multi_c2r, fftw_plan_multi_r2c
       COMPLEX(KIND=8), POINTER, DIMENSION(:,:,:) :: cu, prod_cu 
       REAL(KIND=8),    POINTER, DIMENSION(:,:,:) :: ru, prod_ru
    END TYPE dyn_FFTW                     ! instead of ALLOCATABLE
      
    INCLUDE 'fftw3.f'
    INCLUDE 'mpif.h'
    ! Format: V_1in(1:np,1:6,1:m_max_c) 
    ! INPUT ARE COSINE AND SINE COEFFICIENTS
    ! THEY ARE PUT IN COMPLEX FORMAT: c_0 = a_0 + i*0 and c_n = (a_n-i*b_n)/2
    REAL(KIND=8), DIMENSION(:,:,:),  INTENT(IN)  :: V1_in, V2_in
    REAL(KIND=8), DIMENSION(:,:,:),  INTENT(OUT) :: V_out
    INTEGER, INTENT(IN)                                 :: pb_type 
    REAL(KIND=8), DIMENSION(:), OPTIONAL, INTENT(INOUT) :: temps


    ! Saved variables
    TYPE(dyn_FFTW), DIMENSION(2), TARGET, SAVE :: FFTW_SAVE 
    ! End saved variables


    INTEGER    :: np, np_tot, bloc_size, nb_field, m_max, m_max_c, N_r
    INTEGER(KIND=8) :: fftw_plan_multi_c2r, fftw_plan_multi_r2c
    LOGICAL :: once
    INTEGER :: i_field, rang, nb_procs, code, MPID
    INTEGER :: nb_valeurs, nb, nf, shiftc, shiftl, index, jindex, longueur_tranche, i, n, p
    REAL(KIND=8) :: t
    !COMPLEX(KIND=8), ALLOCATABLE, DIMENSION(:,:,:) :: cu, prod_cu 
    !REAL(KIND=8),    ALLOCATABLE, DIMENSION(:,:,:) :: ru, prod_ru
    REAL(KIND=8),    ALLOCATABLE, DIMENSION(:,:,:) :: dist_field, combined_field 
    COMPLEX(KIND=8), ALLOCATABLE, DIMENSION(:,:,:) :: combined_prod_cu, dist_prod_cu, out_prod_cu

    !Vectors to speed up the format changes
    COMPLEX(KIND=8), ALLOCATABLE, DIMENSION(:,:)  :: intermediate
    !REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:) :: V_int
    !Vectors to speed up the format changes

    ! FFTW parameters
    INTEGER   :: fft_dim, howmany, istride, ostride, idist, odist
    INTEGER, DIMENSION(1) :: dim, inembed, onembed
    ! Recall complexes must be rescaled
    ! End FFTW parameters

    CALL MPI_COMM_SIZE(MPI_COMM_WORLD,nb_procs,code)
    CALL MPI_COMM_RANK(MPI_COMM_WORLD,rang,code)

    !Temps(1) = Temps de communication
    !Temps(2) = Temps de calcul
    !Temps(3) = Temps de changement de format
    IF (PRESENT(temps)) temps = 0.d0

    IF (SIZE(V1_in,1) .NE. FFTW_SAVE(pb_type)%np) THEN
       once = .TRUE.
    ELSE
       once = .FALSE.
    END IF

    IF (once) THEN 
       np      = SIZE(V1_in,1) ! Number of points in the meridian plane
       nb_field= SIZE(V1_in,2) ! Number of fields
       m_max_c = SIZE(V1_in,3) ! Number of complex (cosines + sines) coefficients per point
       m_max = m_max_c*nb_procs! Number of comlex coefficients per point per processor
       N_r=2*m_max-1           ! Number of Real coefficients per point
       IF (MOD(nb_field,2)/=0 .OR. m_max_c==0) THEN
          WRITE(*,*) ' BUG '
          STOP
       END IF

       ! Bloc_size is the number of points that are handled by one processor 
       ! once the Fourier modes are all collected
       ! Computation of bloc_size and np_tot
       IF (MODULO(np,nb_procs)==0) THEN
          bloc_size = np/nb_procs
       ELSE
          bloc_size = np/nb_procs + 1
       END IF
       np_tot = nb_procs*bloc_size
       ALLOCATE(FFTW_SAVE(pb_type)%cu(m_max,nb_field,bloc_size))
       ALLOCATE(FFTW_SAVE(pb_type)%ru(N_r,  nb_field,bloc_size))
       ALLOCATE(FFTW_SAVE(pb_type)%prod_cu(m_max,nb_field/2,bloc_size))
       ALLOCATE(FFTW_SAVE(pb_type)%prod_ru(N_r,  nb_field/2,bloc_size))
    ELSE
       IF ((SIZE(V1_in,1).NE.FFTW_SAVE(pb_type)%np) .OR. (SIZE(V1_in,2).NE.FFTW_SAVE(pb_type)%nb_field) &
            .OR. (SIZE(V1_in,3).NE.FFTW_SAVE(pb_type)%m_max_c)) THEN
          WRITE(*,*) ' BUG IN FFT, Something wrong happened'
          STOP
       END IF
       np       = FFTW_SAVE(pb_type)%np
       nb_field = FFTW_SAVE(pb_type)%nb_field
       m_max_c  = FFTW_SAVE(pb_type)%m_max_c
       m_max    = FFTW_SAVE(pb_type)%m_max
       N_r      = FFTW_SAVE(pb_type)%N_r
       bloc_size= FFTW_SAVE(pb_type)%bloc_size
       np_tot   = FFTW_SAVE(pb_type)%np_tot
       fftw_plan_multi_c2r = FFTW_SAVE(pb_type)%fftw_plan_multi_c2r
       fftw_plan_multi_r2c = FFTW_SAVE(pb_type)%fftw_plan_multi_r2c
    END IF

    !t = MPI_WTIME()
    !ALLOCATE(V_int(nb_field,np))
    !ALLOCATE(V_int(np,nb_field))

    ALLOCATE(intermediate(nb_field/2,bloc_size))
    ALLOCATE(    dist_field(m_max_c,2*nb_field,np_tot))
    ALLOCATE(combined_field(m_max_c,2*nb_field,np_tot))
    ALLOCATE(dist_prod_cu(nb_field/2,bloc_size,m_max), combined_prod_cu(nb_field/2,bloc_size,m_max))
    ALLOCATE(out_prod_cu(m_max_c,np_tot,nb_field/2))
    !ALLOCATE(cu(m_max,nb_field,bloc_size))
    !ALLOCATE(ru(N_r,  nb_field,bloc_size))
    !ALLOCATE(prod_cu(m_max,nb_field/2,bloc_size))
    !ALLOCATE(prod_ru(N_r,  nb_field/2,bloc_size))
    !WRITE(*,*) ' time to allocate', MPI_WTIME() -t

    ! Packing all 3 complex components of both v1 and v2 input fields
    ! into dist_field, where the dimension indexing the nodal points varies the least rapidly,
    ! so that after distributing the data to the processes, each one will obtain a part 
    ! on nodal points
    t = MPI_WTIME()
    DO i = 1, m_max_c
       !V_int = V1_in(:,:,i)
       !dist_field(i,1:nb_field,:) = transpose(V_int)
       !V_int = V2_in(:,:,i)
       !dist_field(i,nb_field+1:2*nb_field,:) = transpose(V_int)
       dist_field(i,1:nb_field,1:np) = TRANSPOSE(V1_in(:,:,i))
       dist_field(i,nb_field+1:2*nb_field,1:np) = TRANSPOSE(V2_in(:,:,i))
       !DO nf = 1, nb_field/2 
       !   dist_field(i,2*nf-1,1:np) = V1_in(:,2*nf-1,i) 
       !   dist_field(i,2*nf  ,1:np) = V1_in(:,2*nf  ,i) 
       !   dist_field(i,nb_field+2*nf-1,1:np) = V2_in(:,2*nf-1,i)
       !   dist_field(i,nb_field+2*nf  ,1:np) = V2_in(:,2*nf,i)
       !END DO
    END DO
    IF (PRESENT(temps)) temps(3) = temps(3) + MPI_WTIME() -t

    IF (np/=np_tot) dist_field(:,:,np+1:np_tot) = 1.d100

    longueur_tranche=bloc_size*m_max_c*nb_field*2

    t = MPI_WTIME()
    MPID=MPI_DOUBLE_PRECISION
    CALL MPI_ALLTOALL (dist_field,longueur_tranche, MPID, combined_field, longueur_tranche, &
         MPID, MPI_COMM_WORLD, code)
    IF (PRESENT(temps)) temps(1) = temps(1) + MPI_WTIME() -t

    t = MPI_WTIME()

    DO n = 1, bloc_size 
       DO nb = 1, nb_procs
          shiftc = (nb-1)*bloc_size  
          shiftl = (nb-1)*m_max_c
          jindex = n + shiftc
          DO nf = 1, nb_field
             ! Put real and imaginary parts in a complex 
             ! for each field, nf=1,2,3,4,5,6
             ! nf=1,2,3 => V1_in
             ! nf=4,5,6 => V2_in
             ! cu(shiftl+1:shiftl+m_max_c,nf,n) = cmplx(combined_field(:,2*nf-1,jindex),combined_field(:,2*nf,jindex))
             ! INPUT ARE COSINE AND SINE COEFFICIENTS
             ! THEY ARE PUT IN COMPLEX FORMAT: c_0 = a_0 + i*0 and c_n = (a_n-i*b_n)/2
             FFTW_SAVE(pb_type)%cu(shiftl+1:shiftl+m_max_c,nf,n) = 0.5d0*CMPLX(combined_field(:,2*nf-1,jindex),-combined_field(:,2*nf,jindex),KIND=8)
          END DO
       END DO
    END DO
    FFTW_SAVE(pb_type)%cu(1,:,:) = 2.d0*CMPLX(REAL(FFTW_SAVE(pb_type)%cu(1,:,:),KIND=8),0.d0,KIND=8)
    IF (PRESENT(temps)) temps(3) = temps(3) + MPI_WTIME() -t

    ! Set the parameters for dfftw
    fft_dim=1; istride=1; ostride=1; 
    idist=N_r;   inembed(1)=N_r; DIM(1)=N_r
    odist=m_max; onembed(1)=m_max
    IF (rang==(nb_procs-1)) THEN
       howmany= (np - bloc_size*(nb_procs-1))*nb_field 
    ELSE 
       howmany=bloc_size*nb_field
    END IF

    t = MPI_WTIME()
    IF (once) CALL dfftw_plan_many_dft_c2r(fftw_plan_multi_c2r, fft_dim, dim, howmany, FFTW_SAVE(pb_type)%cu, &
         onembed, ostride, odist, FFTW_SAVE(pb_type)%ru, inembed, istride, idist, FFTW_ESTIMATE);
    CALL dfftw_execute(fftw_plan_multi_c2r)

    ! SIMPLE PRODUCT
    !DO nf = 1, nb_field/2
    !   ! ru(1:3) contains V1_in and ru(4:6) contains V2_in
    !   prod_ru(:,nf,:)  = ru(:,nf,:)*ru(:,nb_field/2+nf,:)
    !END DO
    ! END SIMPLE PRODUCT

    ! CROSS PRODDUCT
    IF (nb_field==6) THEN
       FFTW_SAVE(pb_type)%prod_ru(:,1,:) = FFTW_SAVE(pb_type)%ru(:,2,:)*FFTW_SAVE(pb_type)%ru(:,6,:) - FFTW_SAVE(pb_type)%ru(:,3,:)*FFTW_SAVE(pb_type)%ru(:,5,:)
       FFTW_SAVE(pb_type)%prod_ru(:,2,:) = FFTW_SAVE(pb_type)%ru(:,3,:)*FFTW_SAVE(pb_type)%ru(:,4,:) - FFTW_SAVE(pb_type)%ru(:,1,:)*FFTW_SAVE(pb_type)%ru(:,6,:)
       FFTW_SAVE(pb_type)%prod_ru(:,3,:) = FFTW_SAVE(pb_type)%ru(:,1,:)*FFTW_SAVE(pb_type)%ru(:,5,:) - FFTW_SAVE(pb_type)%ru(:,2,:)*FFTW_SAVE(pb_type)%ru(:,4,:)
    END IF
    ! CROSS PRODUCT

    howmany = howmany/2
    IF (once) CALL dfftw_plan_many_dft_r2c(fftw_plan_multi_r2c, fft_dim, dim, howmany, FFTW_SAVE(pb_type)%prod_ru, &
         inembed, istride, idist, FFTW_SAVE(pb_type)%prod_cu, onembed, ostride, odist, FFTW_ESTIMATE);
    CALL dfftw_execute(fftw_plan_multi_r2c)

    FFTW_SAVE(pb_type)%prod_cu = FFTW_SAVE(pb_type)%prod_cu/N_r !Scaling 
    !   prod_cu = prod_cu/REAL(N_r,KIND=8) !Scaling 
    IF (PRESENT(temps)) temps(2) = temps(2) + MPI_WTIME() -t


    !Now we need to redistribute the Fourier coefficients on each processor
    t = MPI_WTIME()
    combined_prod_cu(:,:,1)=FFTW_SAVE(pb_type)%prod_cu(1,:,:)
    DO n=2, m_max
       !combined_prod_cu(:,:,n)=prod_cu(n,:,:)
       combined_prod_cu(:,:,n)=2.d0*CONJG(FFTW_SAVE(pb_type)%prod_cu(n,:,:))
    END DO
    IF (PRESENT(temps)) temps(3) = temps(3) + MPI_WTIME() -t

    t = MPI_WTIME()
    longueur_tranche=bloc_size*m_max_c*nb_field
    MPID=MPI_DOUBLE_PRECISION
    CALL MPI_ALLTOALL (combined_prod_cu,longueur_tranche,MPID, dist_prod_cu,longueur_tranche, &
         MPID,MPI_COMM_WORLD,code)
    IF (PRESENT(temps)) temps(1) = temps(1) + MPI_WTIME() -t

    ! dimensions: 
    ! v_out(np, nb_field, m_max_c)
    t = MPI_WTIME()

    IF (.FALSE.) THEN
       DO i_field = 1, nb_field/2
          DO n = 1, bloc_size
             DO nb = 1, nb_procs
                shiftc = (nb-1)*bloc_size
                shiftl = (nb-1)*m_max_c
                out_prod_cu(:,n+shiftc,i_field) = dist_prod_cu(i_field,n,shiftl+1:shiftl+m_max_c)
             END DO
          END DO
       END DO

       DO i_field = 1, nb_field/2
          DO i = 1, m_max_c
             v_out(:, i_field*2-1, i) = REAL(out_prod_cu(i, 1:np, i_field),KIND=8)
             v_out(:, i_field*2,   i) = AIMAG(out_prod_cu(i, 1:np, i_field))
          END DO
       END DO
    ELSE
       DO i = 1, m_max_c
          DO nb = 1, nb_procs
             shiftc = (nb-1)*bloc_size  
             shiftl = (nb-1)*m_max_c
             intermediate = dist_prod_cu(:,:,shiftl+i)
             DO n = 1, bloc_size
                IF (n+shiftc > np ) CYCLE 
                DO i_field = 1, nb_field/2
                   v_out(n+shiftc, i_field*2-1, i) = REAL (intermediate(i_field,n),KIND=8)
                   v_out(n+shiftc, i_field*2 , i)  = AIMAG(intermediate(i_field,n))
                   !v_out(n+shiftc, i_field*2-1, i) = REAL(dist_prod_cu(i_field,n,shiftl+i),KIND=8)
                   !v_out(n+shiftc, i_field*2 , i)  = AIMAG(dist_prod_cu(i_field,n,shiftl+i))
                END DO
             END DO
          END DO
       END DO
    END IF

    IF (PRESENT(temps)) temps(3) = temps(3) + MPI_WTIME() -t

    DEALLOCATE(dist_field, combined_field, dist_prod_cu, combined_prod_cu, intermediate) 

    IF (once) THEN
       FFTW_SAVE(pb_type)%np = np
       FFTW_SAVE(pb_type)%nb_field = nb_field
       FFTW_SAVE(pb_type)%m_max_c = m_max_c
       FFTW_SAVE(pb_type)%m_max = m_max
       FFTW_SAVE(pb_type)%N_r = N_r
       FFTW_SAVE(pb_type)%bloc_size = bloc_size
       FFTW_SAVE(pb_type)%np_tot = np_tot
       FFTW_SAVE(pb_type)%fftw_plan_multi_c2r=fftw_plan_multi_c2r
       FFTW_SAVE(pb_type)%fftw_plan_multi_r2c=fftw_plan_multi_r2c
    END IF

  END SUBROUTINE FFT_PAR_CROSS_PROD_mod


END MODULE sft_parallele
