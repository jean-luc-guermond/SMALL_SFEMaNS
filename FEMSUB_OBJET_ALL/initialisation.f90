!
!
!Authors Jean-Luc Guermond, Caroline Nore, Copyrights 2005
!Revised June 2008, Jean-Luc Guermond
!
MODULE initialisation

  USE def_type_mesh
  IMPLICIT NONE

  PUBLIC:: initial, navier_stokes, maxwell, mhd, post_proc_test, sauvegarde
  PUBLIC:: prodmat_maxwell, prodmat_maxwell_int_by_parts
  PRIVATE 

  !Champs pour Navier-Stokes-------------------------------------------------
  TYPE(mesh_type), TARGET                         :: pp_mesh, vv_mesh    
  REAL(KIND=8), TARGET, ALLOCATABLE, DIMENSION(:,:,:)          :: un, un_m1
  ! (noeuds,type,mode) composante du champ de vitesse a deux instants sur vv_mesh
  REAL(KIND=8), TARGET, ALLOCATABLE, DIMENSION(:,:,:)          :: pn, pn_m1
  ! (noeuds,type,mode) composante du champ de pression a deux instants sur pp_mesh
  REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:,:)     :: incpn, incpn_m1
  !---------------------------------------------------------------------------

  !Champs pour Maxwell--------------------------------------------------------
  REAL(KIND=8), TARGET, ALLOCATABLE, DIMENSION(:,:,:):: Hn, Hn1, Hext, phin, phin1, vel
  REAL(KIND=8), TARGET, ALLOCATABLE, DIMENSION(:) :: sigma_field, mu_H_field
  TYPE(mesh_type), TARGET                         :: H_mesh, phi_mesh
  TYPE(interface_type), TARGET                    :: interface_H_mu, interface_H_phi

  TYPE(periodic_type)                             :: H_phi_per
  TYPE(periodic_type)                             :: vvrt_per
  TYPE(periodic_type)                             :: vvz_per
  TYPE(periodic_type)                             :: pp_per

  !---------------------------------------------------------------------------

  !Variables de couplage------------------------------------------------------
  REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:,:)     :: Hn_p2
  CHARACTER(len=3)                                :: type_pb
  !October 7, 2008, JLG
  INTEGER, ALLOCATABLE, DIMENSION(:)              :: jj_v_to_H
  !October 7, 2008, JLG
  !---------------------------------------------------------------------------

  !liste des modes------------------------------------------------------------
  INTEGER,      TARGET, ALLOCATABLE, DIMENSION(:) :: list_mode
  INTEGER                                         :: code, rang, nb_procs, m_max, m_max_c
  !---------------------------------------------------------------------------

  !Noms reserves--------------------------------------------------------------
  INTEGER                                         :: nb_iteration 
  REAL(KIND=8)                                    :: dt, time, Re, Rem
  REAL(KIND=8)                                    :: mu_phi, R_fourier
  REAL(KIND=8), DIMENSION(3)                      :: stab 
  INTEGER                                         :: index_fourier
  CHARACTER(len=64)                               :: data_fichier 
  LOGICAL                                         :: second_order_ext_pressure, dom_H_larger_dom_ns 
  !---------------------------------------------------------------------------

  !Gestion locale-------------------------------------------------------------
  INTEGER,                   DIMENSION(2)         :: nb_syst_lin_ns, nb_syst_lin_mx
  !---------------------------------------------------------------------------

  !Nom du fichier pour restart-----------------------------------------------
  CHARACTER(len=64)                               :: directory, file_name
  !---------------------------------------------------------------------------

  !-------------END OF DECLARATIONS------------------------------------------

CONTAINS 
  !---------------------------------------------------------------------------


  !---------------------------------------------------------------------------
  SUBROUTINE initial(vv_mesh_out, pp_mesh_out, H_mesh_out, phi_mesh_out,&
       interface_H_phi_out, interface_H_mu_out, list_mode_out, & 
       un_out, pn_out, Hn_out, phin_out, vel_out, mu_H_field_out, sigma_field_out, mu_phi_out, &
       Re_out, Rem_out, time_out, dt_out, nb_iteration_out, m_max_c_out, type_pb_out, &
       test_de_convergence_out) 
    USE boundary 
    IMPLICIT NONE
    TYPE(mesh_type), POINTER                :: pp_mesh_out, vv_mesh_out
    TYPE(mesh_type), POINTER                :: H_mesh_out, phi_mesh_out
    TYPE(interface_type), POINTER           :: interface_H_mu_out, interface_H_phi_out
    INTEGER,      POINTER,  DIMENSION(:)    :: list_mode_out
    REAL(KIND=8), POINTER,  DIMENSION(:,:,:):: un_out, pn_out, Hn_out, phin_out, vel_out
    REAL(KIND=8), POINTER,  DIMENSION(:)    :: sigma_field_out, mu_H_field_out
    REAL(KIND=8)                            :: mu_phi_out, Re_out, Rem_out, time_out, dt_out
    INTEGER                                 :: nb_iteration_out, m_max_c_out
    CHARACTER(len=3)                        :: type_pb_out
    LOGICAL                                 :: test_de_convergence_out

    CALL INIT
    vv_mesh_out => vv_mesh
    pp_mesh_out => pp_mesh
    H_mesh_out => H_mesh
    phi_mesh_out => phi_mesh
    interface_H_mu_out => interface_H_mu
    interface_H_phi_out => interface_H_phi
    list_mode_out => list_mode
    un_out   => un
    pn_out   => pn
    Hn_out   => Hn
    phin_out => phin
    vel_out => vel
    mu_H_field_out => mu_H_field
    sigma_field_out => sigma_field 
    mu_phi_out = mu_phi
    Re_out   = Re
    Rem_out  = Rem
    time_out = time
    dt_out   = dt
    nb_iteration_out = nb_iteration
    m_max_c_out = m_max_c
    type_pb_out = type_pb
    test_de_convergence_out = test_de_convergence
  END SUBROUTINE initial
  !---------------------------------------------------------------------------


  !---------------------------------------------------------------------------
  SUBROUTINE navier_stokes(time_in, un_out, pn_out)
    USE subroutine_ns
    IMPLICIT NONE
    REAL(KIND=8), INTENT(in)                               :: time_in
    REAL(KIND=8), POINTER,  DIMENSION(:,:,:)               :: un_out, pn_out

    time = time_in
    CALL three_level_ns(time, dt, Re, list_mode, pp_mesh, &
         vv_mesh, incpn_m1, incpn, pn_m1, pn, un_m1, un, nb_syst_lin_ns, &
         vvrt_per, vvz_per, pp_per, data_fichier, Hn_p2, second_order_ext_pressure)
    un_out => un 
    pn_out => pn

  END SUBROUTINE navier_stokes
  !---------------------------------------------------------------------------


  !---------------------------------------------------------------------------
  SUBROUTINE maxwell(time_in, Hn_out, phin_out)
    USE update_maxwell 
    USE boundary
    IMPLICIT NONE
    REAL(KIND=8), INTENT(in)                               :: time_in
    REAL(KIND=8), POINTER,  DIMENSION(:,:,:)               :: Hn_out, phin_out
    INTEGER :: i

    time = time_in

    CALL maxwell_decouple(H_mesh, phi_mesh, interface_H_phi, interface_H_mu, Hn, phin, Hn1, phin1, vel, &
         stab, sigma_field, R_fourier, index_fourier, mu_H_field, mu_phi, &
         time, dt, Rem, list_mode, nb_syst_lin_mx, H_phi_per, data_fichier)

    Hn_out => Hn 
    phin_out => phin

  END SUBROUTINE maxwell
  !---------------------------------------------------------------------------


  !---------------------------------------------------------------------------
  SUBROUTINE mhd(time_in, un_out, pn_out, Hn_out, phin_out, select_mode_ns, select_mode_mxw)
    USE subroutine_ns
    USE update_maxwell 
    USE boundary
    IMPLICIT NONE
    REAL(KIND=8), INTENT(in)                               :: time_in
    REAL(KIND=8), POINTER,  DIMENSION(:,:,:):: un_out, pn_out, Hn_out, phin_out
    INTEGER, DIMENSION(:), OPTIONAL, INTENT(IN) :: select_mode_ns, select_mode_mxw
    INTEGER :: i, k, m, vel_np, ms, ms1, ms2, j1, j2, ns, m1, m2

    IF (PRESENT(select_mode_mxw)) THEN
       Hn(:,:,select_mode_mxw) = 0.d0
       phin(:,:,select_mode_mxw) = 0.d0
       Hn1(:,:,select_mode_mxw) = 0.d0
       phin1(:,:,select_mode_mxw) = 0.d0
    END IF

    time = time_in
    DO i = 1, m_max_c
       DO k = 1, 6
          Hext(:,k,i) = (2*Hn(:,k,i)-Hn1(:,k,i))*SQRT(mu_H_field)
       END DO
    END DO

    !IF (H_mesh%gauss%n_w==vv_mesh%gauss%n_w) THEN
    !   Hn_p2 = Hext(1:vv_mesh%np,:,:)
    !ELSE !Injection champ magnetique P1 sur maillage P2
    !   CALL inject_champ(H_mesh%jj, vv_mesh%jj, Hext, Hn_p2)
    !END IF
    !October 7, 28 JLG
    Hn_p2(jj_v_to_H,:,:) = Hext(1:SIZE(jj_v_to_H),:,:)
    DO ms = 1, interface_H_mu%mes
       ms2 = interface_H_mu%mesh2(ms)
       ms1 = interface_H_mu%mesh1(ms)
       m2 = H_mesh%neighs(ms2)
       m1 = H_mesh%neighs(ms1)
       IF (m1> vv_mesh%me .OR. m2> vv_mesh%me) CYCLE
       DO ns = 1, H_mesh%gauss%n_ws
          j1 =  interface_H_mu%jjs1(ns,ms)
          j2 =  interface_H_mu%jjs2(ns,ms)
          IF (jj_v_to_H(j1)/=jj_v_to_H(j2)) THEN
             write(*,*) jj_v_to_H(j1), jj_v_to_H(j2)
             WRITE(*,*) ' BUG'
             STOP
          END IF
          Hn_p2(jj_v_to_H(j1),:,:) = (Hext(j1,:,:)+Hext(j2,:,:))/2
       END DO
    END DO

    IF (H_mesh%gauss%n_w/=vv_mesh%gauss%n_w) THEN
       DO m = 1, vv_mesh%me
          Hn_p2(vv_mesh%jj(4,m),  :,:) = (Hn_p2(vv_mesh%jj(2,m),:,:) + Hn_p2(vv_mesh%jj(3,m),:,:))/2
          Hn_p2(vv_mesh%jj(5,m),  :,:) = (Hn_p2(vv_mesh%jj(3,m),:,:) + Hn_p2(vv_mesh%jj(1,m),:,:))/2
          Hn_p2(vv_mesh%jj(6,m),  :,:) = (Hn_p2(vv_mesh%jj(1,m),:,:) + Hn_p2(vv_mesh%jj(2,m),:,:))/2 
       END DO
    END IF
    !October 7, JLG

    !May 20, 2009, JLG
    !Moved this line before calling three_level_ns 
    !plus put incpn and incpn_m1 to zero for modes in slect_mode_ns
    IF (PRESENT(select_mode_ns)) THEN
       un(:,:,select_mode_ns) = 0.d0
       pn(:,:,select_mode_ns) = 0.d0
       incpn(:,:,select_mode_ns) = 0.d0
       un_m1(:,:,select_mode_ns) = 0.d0
       pn_m1(:,:,select_mode_ns) = 0.d0
       incpn_m1(:,:,select_mode_ns) = 0.d0
    END IF
    !May 20, 2009, JLG
    CALL three_level_ns(time, dt, Re, list_mode, pp_mesh, &
         vv_mesh, incpn_m1, incpn, pn_m1, pn, un_m1, un, nb_syst_lin_ns, &
         vvrt_per, vvz_per, pp_per, data_fichier, Hn_p2, second_order_ext_pressure)


    !October 7, 2008, JLG
    !IF (H_mesh%gauss%n_w==vv_mesh%gauss%n_w) THEN
    !vel(1:vv_mesh%np,:,:) = un !BUG corrected June 4, 2008, JLG: 2*un-un_m1
    !vel_np = vv_mesh%np
    !ELSE !Projection du champ de vitesse P2 sur le maillage P1
    !BUG corrected June 4, 2008, JLG
    !vel_np = pp_mesh%np
    !CALL project_champ(un, vel(1:vel_np,:,:), vv_mesh, pp_mesh)
    !END IF
    vel_np = SIZE(jj_v_to_H)
    vel(1:vel_np,:,:) = un(jj_v_to_H,:,:)
    !October 7, 2008, JLG

    IF (dom_H_larger_dom_ns) THEN !We extend vel
       DO i = 1, m_max_c
          DO k= 1, 6 !The user has to code extension_vel 
             vel(vel_np+1:,k,i) = extension_vel(k, H_mesh, list_mode(i), time, vel_np+1)
          END DO
       END DO
    END IF

    CALL maxwell_decouple(H_mesh, phi_mesh, interface_H_phi, interface_H_mu, Hn, phin, Hn1, phin1, vel, &
         stab, sigma_field, R_fourier, index_fourier, mu_H_field, mu_phi, &
         time, dt, Rem, list_mode, nb_syst_lin_mx, H_phi_per, data_fichier)

    un_out => un 
    pn_out => pn
    Hn_out => Hn 
    phin_out => phin

  END SUBROUTINE mhd
  !---------------------------------------------------------------------------

  SUBROUTINE post_proc_test(time_in)
    USE fem_tn_NS_MHD
    USE tn_parallele
    USE sub_plot
    USE boundary
    USE boundary_anal

    IMPLICIT NONE
    REAL(KIND=8), INTENT(IN) :: time_in
    INTEGER :: i, k, code
    REAL(KIND=8) :: err, norm

    time = time_in
    IF(test_conv_ns) THEN
       DO i = 1, m_max_c
          DO k= 1, 6
             un_m1(:,k,i) = un(:,k,i) - vv_exact(k,vv_mesh%rr,list_mode(i),time)
          END DO
          DO k= 1, 2
             pn_m1(:,k,i) = pn(:,k,i) - pp_exact(k,pp_mesh%rr,list_mode(i),time)
          END DO
       END DO
       IF (rang==0) THEN
          WRITE(rang+10,*) 'CHAMP de Vitesse #####################'
          WRITE(rang+10,*) 'Erreur L1   vitesse   = ', norme_max_champ_par(un_m1(:,:,:))
          WRITE(rang+10,*) 'erreur L2   vitesse   = ', norme_l2_champ_par(vv_mesh, list_mode, un_m1)
          WRITE(rang+10,*) 'erreur H1   vitesse   = ', norme_H1_champ_par(vv_mesh, list_mode, un_m1)
          WRITE(rang+10,*) 'norme  L2 divergence  = ', norme_div_par(vv_mesh, list_mode, un)
          WRITE(rang+10,*) 'CHAMP de pression #####################'
          WRITE(rang+10,*) 'Erreur L1   pression  = ', norme_max_champ_par(pn_m1(:,:,:))
          WRITE(rang+10,*) 'erreur L2   pression  = ', norme_l2_champ_par(pp_mesh, list_mode, pn_m1)
          !CALL plot_scalar_field(pp_mesh%jj, pp_mesh%rr, pn_m1(:,1,1),'errp11.plt')
          !CALL plot_scalar_field(pp_mesh%jj, pp_mesh%rr, pn_m1(:,2,1),'errp21.plt')
       END IF
    END IF

    IF (test_conv_maxwell) THEN
       CALL norme_interface(H_mesh,phi_mesh,interface_H_phi,mu_H_field,mu_phi,Hn(:,:,1),phin(:,:,1),list_mode(1),err)
       WRITE(*,*) ' erreur on the interface H_phi final', err
       CALL norme_interface_H_mu(H_mesh,interface_H_mu,mu_H_field,Hn(:,:,1),list_mode(1),err)
       WRITE(*,*) ' erreur on the interface H_mu final', err
       DO i = 1, m_max_c
          DO k =1, 6
             Hn1(:,k,i) = Hn(:,k,i) - Hexact(k, H_mesh%rr, list_mode(i), mu_H_field, time)
          END DO
          DO k =1, 2
             phin1(:,k,i) = phin(:,k,i) - Phiexact(k, phi_mesh%rr, list_mode(i), mu_phi, time)
          END DO
       END DO
       !IF (rang==0)  CALL plot_scalar_field(H_mesh%jj, H_mesh%rr, Hn1(:,1,1),'err_Hn1.plt')
       !IF (rang==0)  CALL plot_scalar_field(H_mesh%jj, H_mesh%rr, Hn1(:,2,1),'err_Hn2.plt')
       !IF (rang==0)  CALL plot_scalar_field(H_mesh%jj, H_mesh%rr, Hn1(:,3,1),'err_Hn3.plt')
       !IF (rang==0)  CALL plot_scalar_field(H_mesh%jj, H_mesh%rr, Hn1(:,4,1),'err_Hn4.plt')
       !IF (rang==0)  CALL plot_scalar_field(H_mesh%jj, H_mesh%rr, Hn1(:,5,1),'err_Hn5.plt')
       !IF (rang==0)  CALL plot_scalar_field(H_mesh%jj, H_mesh%rr, Hn1(:,6,1),'err_Hn6.plt')
       !IF (rang==0)  CALL plot_scalar_field(phi_mesh%jj, phi_mesh%rr, phin1(:,1,1),'err_phin1.plt')
       !IF (rang==0)  CALL plot_scalar_field(phi_mesh%jj, phi_mesh%rr, phin1(:,2,1),'err_phin2.plt')
       WRITE(rang+10,*) 'CHAMP Hn #####################'
       norm = norme_L2_champ_par(H_mesh,list_mode,Hn)
       err  = norme_L2_champ_par(H_mesh,list_mode,Hn1)
       WRITE(rang+10,*) 'erreur finale Hn   L2       = ', err/norm
       norm = norm + norme_H1_champ_par(H_mesh,list_mode,Hn)
       !June 17 2008 norm = norme_H1_champ_par(H_mesh,list_mode,Hn)
       err  = norme_curl_par(H_mesh, list_mode, Hn1)
       WRITE(rang+10,*) 'erreur finale rot(Hn)   L2  = ', err/norm
       DO i = 1, m_max_c
          DO k = 1, 6
             Hn1(:,k,i) = mu_H_field(:)*Hn(:,k,i)
          ENDDO
       ENDDO
       norm = norme_L2_champ_par(H_mesh,list_mode,Hn1) &
            + norme_H1_champ_par(H_mesh,list_mode,Hn1)
       err = norme_div_par(H_mesh,list_mode, Hn1)
       WRITE(rang+10,*) 'erreur finale div(mu Hn)   L2  = ', err/norm
       !
       WRITE(rang+10,*) 'CHAMP phin #####################'
       norm = norme_H1_champ_par(phi_mesh,list_mode,phin)
       err  = norme_H1_champ_par(phi_mesh,list_mode,phin1)
       WRITE(rang+10,*) 'erreur final phin  H1  = ', err/norm
    END IF


    !STOP

  END SUBROUTINE post_proc_test
  !---------------------------------------------------------------------------


  !----------------SAUVEGARDE-------------------------------------------------
  SUBROUTINE sauvegarde(it, freq_restart) 
    USE restart   
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: it, freq_restart

    IF (type_pb=='nst' .OR. type_pb=='mhd') THEN
       CALL write_restart_ns(vv_mesh, pp_mesh, time, list_mode, un, un_m1, pn, pn_m1, &
            incpn, incpn_m1, file_name, it, freq_restart)
    END IF
    IF (type_pb=='mxw' .OR. type_pb=='mhd' .OR. type_pb=='mxx') THEN
       CALL write_restart_maxwell(H_mesh, phi_mesh, time, list_mode, Hn, Hn1, &
            phin, phin1, file_name, it, freq_restart)
    END IF
  END SUBROUTINE sauvegarde
  !---------------------------------------------------------------------------

  SUBROUTINE INIT
    !==================

    USE chaine_caractere
    USE periodic
    USE prep_maill
    USE prep_mesh_interface
    USE subroutine_ns
    USE restart   
    USE boundary
    USE boundary_anal
    USE sub_plot

    IMPLICIT NONE
    INCLUDE 'mpif.h'
    INTEGER, DIMENSION(:), ALLOCATABLE      :: list_dom_ns, list_dom_H, list_dom_H_in, list_dom_phi, list_inter_H_phi
    INTEGER, DIMENSION(:), ALLOCATABLE      :: list_inter_mu, H_in_to_new
    INTEGER                                 :: nb_dom_ns, nb_dom_H, nb_dom_phi, nb_inter
    INTEGER                                 :: nb_inter_mu
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:) :: sigma, mu_H
    INTEGER                                 :: type_fe_H, type_fe_phi
    INTEGER, DIMENSION(:), ALLOCATABLE      :: list_mode_lect
    LOGICAL                                 :: select_mode
    INTEGER                                 :: nb_select_mode
    CHARACTER(len=64)                       :: repertoire, fichier, data_directory, data_file
    INTEGER                                 :: nb_bord_periodic
    LOGICAL                                 :: irestart_u=.FALSE., irestart_h=.FALSE.
    LOGICAL                                 :: iformatted, periodic_r
    INTEGER                                 :: k, kp, m, n, i, d_end, f_end, numero_du_test, vel_np 
    INTEGER                                 :: code, rang, nb_procs
    REAL(KIND=8)                            :: time_u, time_h, error

    CALL MPI_COMM_SIZE(MPI_COMM_WORLD,nb_procs,code)
    CALL MPI_COMM_RANK(MPI_COMM_WORLD,rang,code)

    !-------------CONTROL OF DATA INPUTS-------------------------------------------
    OPEN(UNIT = 21, FILE = 'data', FORM = 'formatted', STATUS = 'unknown')
    !------------------------------------------------------------------------------

    !-------------DECIDE WHETHER DEBUGGING OR NOT----------------------------------
    OPEN(UNIT = 21, FILE = 'data', FORM = 'formatted', STATUS = 'unknown')
    CALL read_until(21, 'data_debug')
    READ (21,*) test_de_convergence
    IF (test_de_convergence) THEN 
       READ(21,*) data_directory 
       READ(21,*) numero_du_test
       IF (numero_du_test==1) THEN
          test_conv_ns=.TRUE.
          test_conv_maxwell=.FALSE.
          no_test = 1
          data_file='data_1'
       ELSE IF (numero_du_test==2) THEN
          test_conv_ns=.TRUE.
          test_conv_maxwell=.FALSE.
          no_test = 2
          data_file='data_2'
       ELSE IF (numero_du_test==3) THEN
          test_conv_ns=.FALSE.
          test_conv_maxwell=.TRUE.
          no_test = 1
          data_file='data_3'
       ELSE IF (numero_du_test==4) THEN
          test_conv_ns=.FALSE.
          test_conv_maxwell=.TRUE.
          no_test = 2
          data_file='data_4'
       ELSE IF (numero_du_test==5) THEN
          test_conv_ns=.FALSE.
          test_conv_maxwell=.TRUE.
          no_test = 1
          data_file='data_5'
       ELSE IF (numero_du_test==6) THEN !Durand, H/H/phi
          test_conv_ns=.FALSE.
          test_conv_maxwell=.TRUE.
          no_test = 3
          data_file='data_6'
       ELSE IF (numero_du_test==7) THEN !Durand, phi/H/phi
          test_conv_ns=.FALSE.
          test_conv_maxwell=.TRUE.
          no_test = 4
          data_file='data_7'
       ELSE
          WRITE(*,*) 'BUG: test non prevu ', numero_du_test
          STOP
       END IF
    ELSE
       data_directory = '.'
       data_file='data'
    END IF
    CLOSE(21)

    d_end = last_c_leng (64, data_directory)
    f_end = last_c_leng (64, data_file)
    data_fichier = data_directory(1:d_end)//'/'//data_file(1:f_end) 
    WRITE(*,*) ' data_fichier ', data_fichier
    !------------------------------------------------------------------------------

    !--------------PARAMETRES DU CALCUL--------------------------------------------
    OPEN(UNIT = 21, FILE = data_fichier, FORM = 'formatted', STATUS = 'unknown')
    CALL read_until(21, 'data_scheme_time')
    READ (21, *) type_pb
    READ (21, *) dt, nb_iteration
    IF (type_pb=='nst') THEN
       READ (21 ,*) irestart_u
    ELSE IF (type_pb=='mhd' .OR. type_pb=='mxw') THEN 
       READ (21 ,*) irestart_u, irestart_h
    ELSE
       WRITE(*, *) ' BUG type of probleme not yet defined ', type_pb
       STOP
    END IF
    IF (type_pb == 'mxw' .AND. irestart_u) THEN
       !It means we read the velocity in a file, i.e. we do not need Vexact
       type_pb = 'mxx' 
    END IF
    !------------------------------------------------------------------------------


    !-------------LE MAILLAGE------------------------------------------------------
    CALL read_until(21, 'data_mesh')
    READ (21, *) iformatted
    IF (test_de_convergence) THEN
       READ (21, *) directory, file_name
       directory = data_directory
    ELSE
       READ (21, *) directory, file_name
    END IF
    CALL read_until(21, 'data_periodic')
    READ (21, *) nb_bord_periodic 
    IF (nb_bord_periodic < 1) THEN
       ns_periodic=.FALSE.; mxw_periodic=.FALSE.
       vvrt_per%n_bord = 0; vvz_per%n_bord = 0; pp_per%n_bord = 0
       H_phi_per%n_bord = 0
    ELSE
       ns_periodic=.TRUE.; mxw_periodic=.TRUE.
    END IF

    !------------------------------------------------------------------------------

    IF (type_pb=='nst' .OR. type_pb=='mhd' .OR. type_pb=='mxx') THEN
       !-------------ORGANISATION DU MAILLAGE NAVIER_STOKES---------------------------
       CALL read_until(21, 'mesh_navier_stokes')
       READ (21, *)  nb_dom_ns 
       ALLOCATE(list_dom_ns(nb_dom_ns))
       READ (21, *)  list_dom_ns
       CLOSE(21)

       CALL load_mesh_free_format_ordered(directory, file_name, list_dom_ns, 1, pp_mesh, iformatted, edge_stab=.TRUE.)
       CALL load_mesh_free_format_ordered(directory, file_name, list_dom_ns, 2, vv_mesh, iformatted, edge_stab=.TRUE.)

       IF (ns_periodic) THEN
          repertoire = data_directory; fichier=data_file
          CALL prep_periodic(repertoire, fichier, pp_mesh, pp_per)
          CALL prep_periodic(repertoire, fichier, vv_mesh, vvz_per)
          CALL prep_periodic_bloc(repertoire, fichier, vv_mesh, vvrt_per, 2)
       END IF

       !------------------------------------------------------------------------------


       !--------------PARAMETRES DU CALCUL--------------------------------------------
       OPEN(UNIT = 21, FILE = data_fichier, FORM = 'formatted', STATUS = 'unknown')
       CALL read_until(21, 'data_scheme_navier_stokes')
       READ (21, *) Re, second_order_ext_pressure
       CLOSE(21)
       !------------------------------------------------------------------------------
    END IF


    IF (type_pb=='mxw' .OR. type_pb=='mhd' .OR. type_pb=='mxx') THEN
       !-------------ORGANISATION DU MAILLAGE MAXWELL---------------------------------
       OPEN(UNIT = 21, FILE = data_fichier, FORM = 'formatted', STATUS = 'unknown')
       CALL read_until(21, 'mesh_maxwell')
       READ (21, *) type_fe_H, type_fe_phi
       READ (21, *) nb_dom_H  ! number of sub_domains for H
       ALLOCATE(list_dom_H_in(nb_dom_H), list_dom_H(nb_dom_H), H_in_to_new(nb_dom_H)) ! JLG/AR Nov 17 2008
       READ (21, *) list_dom_H_in
       READ (21, *) nb_dom_phi  ! number of sub_domains for phi
       ALLOCATE(list_dom_phi(nb_dom_phi))
       READ (21, *) list_dom_phi  
       READ (21, *) nb_inter  ! number of interfaces between H and phi
       ALLOCATE(list_inter_H_phi(nb_inter))
       READ (21, *) list_inter_H_phi
       READ (21, *) mu_phi
       ALLOCATE(mu_H(nb_dom_H))
       READ (21, *) mu_H
       ALLOCATE(sigma(nb_dom_H))
       READ (21, *) sigma

       CALL read_until(21, 'interface_mu')
       READ (21, *) nb_inter_mu
       ALLOCATE(list_inter_mu(nb_inter_mu))
       READ (21, *) list_inter_mu
       CLOSE(21)

       IF (type_pb=='mhd' .OR.  type_pb=='mxx') THEN
          IF (SIZE(list_dom_H) < SIZE(list_dom_ns)) THEN
             WRITE(*,*) ' BUG: NS must be a subset of Maxwell ' 
             STOP
          END IF
          DO k = 1, nb_dom_ns
             ! JLG/AR Nov 17 2008 
             IF (MINVAL(ABS(list_dom_H_in - list_dom_ns(k))) /= 0) THEN
                WRITE(*,*) ' BUG : NS must be a subset of Maxwell '
                STOP
             END IF
             DO kp = 1, nb_dom_H
                IF (list_dom_H_in(kp) == list_dom_ns(k)) EXIT  
             END DO
             H_in_to_new(k) = kp
             ! JLG/AR Nov 17 2008 
             list_dom_H(k) = list_dom_ns(k)
          END DO
          m = nb_dom_ns 
          DO k = 1, nb_dom_H
             IF (MINVAL(ABS(list_dom_H_in(k) - list_dom_ns)) == 0) CYCLE 
             m = m + 1
             ! JLG/AR Nov 17 2008 
             H_in_to_new(m) = k
             ! JLG/AR Nov 17 2008 
             list_dom_H(m) = list_dom_H_in(k)
          END DO
          IF (m/=nb_dom_H) THEN
             WRITE(*,*) ' BUG : m/=nb_dom_H ' 
             STOP
          END IF
          IF (nb_dom_H > nb_dom_ns) THEN
             dom_H_larger_dom_ns = .TRUE.
          ELSE
             dom_H_larger_dom_ns = .FALSE.
          END IF
       ELSE
          ! JLG/AR Nov 17 2008 
          DO k = 1, nb_dom_H
             H_in_to_new(k) = k
          END DO
          ! JLG/AR Nov 17 2008 
          list_dom_H = list_dom_H_in
       END IF

       CALL load_dg_mesh_free_format(directory, file_name, list_dom_H,  list_inter_mu, type_fe_H, H_mesh,   iformatted)
       ! June 13 2008
       !WRITE(*,*) ' ATTENTION:::::::: load_dg_mesh_free_format ENLEVE'
       !CALL load_mesh_free_format_ordered(directory, file_name, list_dom_H, type_fe_H, H_mesh, iformatted)  
       !WRITE(*,*) ' ATTENTION:::::::: load_dg_mesh_free_format ENLEVE'
       ! June 13 2008
       CALL load_mesh_free_format_ordered(directory, file_name, list_dom_phi, &
            type_fe_phi, phi_mesh, iformatted)  

       !June 17 2008, for test 6
       H_mesh_analytic => H_mesh
       !June 17 2008

       IF (type_pb=='mhd' .OR.  type_pb=='mxx') THEN
          ! Verify if the two meshes coincide on the NS domain
          error = 0.d0
          DO k = 1, 2
             DO n = 1, H_mesh%gauss%n_w
                error = error + MAXVAL(ABS(vv_mesh%rr(k,vv_mesh%jj(n,:))-H_mesh%rr(k,H_mesh%jj(n,1:vv_mesh%me))))
             END DO
          END DO
          IF (error/MAXVAL(ABS(H_mesh%rr(1,1) -H_mesh%rr(1,:))) .GE. 5.d-14) THEN
             WRITE(*,*) ' BUG : the two meshes do not coincide on the NS domain',&
                  error/MAXVAL(H_mesh%rr(1,1) -H_mesh%rr(1,:))
             STOP
          END IF

          error = error + MAXVAL(ABS(vv_mesh%rr(1,vv_mesh%jj(4,1:vv_mesh%me)) &
               -(H_mesh%rr(1,H_mesh%jj(2,1:vv_mesh%me))+H_mesh%rr(1,H_mesh%jj(3,1:vv_mesh%me)))/2))&
               + MAXVAL(ABS(vv_mesh%rr(1,vv_mesh%jj(5,:)) &
               -(H_mesh%rr(1,H_mesh%jj(3,1:vv_mesh%me))+H_mesh%rr(1,H_mesh%jj(1,1:vv_mesh%me)))/2))&
               + MAXVAL(ABS(vv_mesh%rr(1,vv_mesh%jj(6,:)) &
               -(H_mesh%rr(1,H_mesh%jj(1,1:vv_mesh%me))+H_mesh%rr(1,H_mesh%jj(2,1:vv_mesh%me)))/2))&
               + MAXVAL(ABS(vv_mesh%rr(2,vv_mesh%jj(4,:)) &
               -(H_mesh%rr(2,H_mesh%jj(2,1:vv_mesh%me))+H_mesh%rr(2,H_mesh%jj(3,1:vv_mesh%me)))/2))&
               + MAXVAL(ABS(vv_mesh%rr(2,vv_mesh%jj(5,:)) &
               -(H_mesh%rr(2,H_mesh%jj(3,1:vv_mesh%me))+H_mesh%rr(2,H_mesh%jj(1,1:vv_mesh%me)))/2))&
               + MAXVAL(ABS(vv_mesh%rr(2,vv_mesh%jj(6,:)) &
               -(H_mesh%rr(2,H_mesh%jj(1,1:vv_mesh%me))+H_mesh%rr(2,H_mesh%jj(2,1:vv_mesh%me)))/2))

          IF (error/MAXVAL(ABS(H_mesh%rr(1,1) -H_mesh%rr(1,:))) .GE. 5.d-14) THEN
             WRITE(*,*) ' WARNING: The two meshes do not coincide on the NS domain.'
             WRITE(*,*) ' WARNING: Either you use curved elements P2 elements or BUG, ', &
                  error/MAXVAL(ABS(H_mesh%rr(1,1) -H_mesh%rr(1,:)))
          END IF
       END IF

       CALL load_interface(H_mesh, H_mesh, list_inter_mu, interface_H_mu, .FALSE.)
       CALL load_interface(H_mesh, phi_mesh, list_inter_H_phi, interface_H_phi, .TRUE.)

       !TTTTTTTTTTTT
       !deallocate(sigma)
       !allocate(sigma(H_mesh%me))
       !sigma = 0
       !sigma(H_mesh%neighs) = H_mesh%sides 
       !sigma(H_mesh%neighs(interface_H_mu%mesh1)) = H_mesh%sides(interface_H_mu%mesh1)
       !sigma(H_mesh%neighs(interface_H_mu%mesh2)) = H_mesh%sides(interface_H_mu%mesh2)
       !sigma(H_mesh%neighs(interface_H_phi%mesh1)) = H_mesh%sides(interface_H_phi%mesh1)
       !IF (rang==0)  CALL plot_const_p1_label(H_mesh%jj, H_mesh%rr, sigma,'sides_H.plt')
       !deallocate(sigma)
       !allocate(sigma(phi_mesh%me))
       !sigma = 0
       !sigma(phi_mesh%neighs) = phi_mesh%sides
       !sigma(phi_mesh%neighs(interface_H_phi%mesh2)) = phi_mesh%sides(interface_H_phi%mesh2)
       !IF (rang==0)  CALL plot_const_p1_label(phi_mesh%jj, phi_mesh%rr, sigma,'sides.plt')
       !stop
       !TTTTTTTTTTTTTTTTT
       !IF (rang==0)  CALL plot_scalar_field(phi_mesh%jj, phi_mesh%rr, phi_mesh%rr(1,:),'r.plt')

       IF (mxw_periodic) THEN
          repertoire = data_directory; fichier=data_file
          CALL prep_periodic_mhd_bc(repertoire, fichier, H_mesh, phi_mesh, H_phi_per)
       END IF

       !------------------------------------------------------------------------------


       !-------------SIGMA POUR MAXWELL-----------------------------------------------
       ALLOCATE(sigma_field(H_mesh%me)) !sigmafield defined on elements of H_mesh
       DO m = 1, H_mesh%me
          DO k=1, nb_dom_H
             IF (H_mesh%i_d(m) == list_dom_H(k)) THEN
                !JLG/AR Nov 17 2008
                !sigma_field(m) = sigma(k)
                sigma_field(m) = sigma(H_in_to_new(k))
                !JLG/AR Nov 17 2008
             ENDIF
          ENDDO
       END DO
       !------------------------------------------------------------------------------

       !-------------MU---------------------------------------------------------------
       ALLOCATE(mu_H_field(H_mesh%np)) !mu_H_field defined at nodes of H_mesh
       DO m = 1, H_mesh%me
          DO k=1, nb_dom_H
             IF (H_mesh%i_d(m) == list_dom_H(k)) THEN
                !JLG/AR Nov 17 2008
                !mu_H_field(H_mesh%jj(:,m)) = mu_H(k)
                mu_H_field(H_mesh%jj(:,m)) = mu_H(H_in_to_new(k))
                !JLG/AR Nov 17 2008
             ENDIF 
          ENDDO
       END DO
       DEALLOCATE(mu_H)
       !------------------------------------------------------------------------------


       !--------------PARAMETRES DU CALCUL--------------------------------------------
       OPEN(UNIT = 21, FILE = data_fichier, FORM = 'formatted', STATUS = 'unknown')
       CALL read_until(21, 'data_scheme_maxwell')
       READ (21, *) Rem
       READ (21, *) R_fourier, index_fourier
       READ (21, *) stab
       !stab(1) = stab_div; stab(2) = stab_colle_H_phi; stab(3) = stab_colle_H_mu
       CLOSE(21)
       !------------------------------------------------------------------------------
    END IF

    !------------CHOIX DES MODES---------------------------------------------------
    OPEN(UNIT = 21, FILE = data_fichier, FORM = 'formatted', STATUS = 'unknown')
    CALL read_until(21, 'select_mode')
    READ (21, *) m_max
    READ (21, *) select_mode, nb_select_mode
    IF (select_mode) THEN
       ALLOCATE(list_mode_lect(nb_select_mode))
       READ(21,*) list_mode_lect
       IF (MOD(nb_select_mode,nb_procs)/= 0) THEN
          WRITE(*,*) ' BUG dans list_mode '
          STOP
       ELSE
          m_max = nb_select_mode ! On ecrase m_max
          m_max_c = m_max/nb_procs
          ALLOCATE(list_mode(m_max_c))
          DO i = 1, m_max_c
             !List_mode contient les modes physiques
             list_mode(i) = list_mode_lect(i + rang*m_max_c) 
          END DO
          DEALLOCATE(list_mode_lect)
       END IF
    ELSE
       m_max_c = m_max/nb_procs
       IF (MOD(m_max,nb_procs)/=0 .OR. m_max_c==0) THEN
          WRITE(*,*) ' BUG '
          STOP
       END IF
       ALLOCATE(list_mode(m_max_c))
       DO i = 1, m_max_c
          !List_mode contient les modes physiques
          list_mode(i) = i + rang*m_max_c - 1
       END DO
    END IF
    CLOSE(21)
    !------------------------------------------------------------------------------


    !-------------COHERENCE DES DONNEES--------------------------------------------
    IF (type_pb=='mhd' .OR. type_pb=='mxx') THEN
       DO m = 1, vv_mesh%me
          IF (MAXVAL(ABS(H_mesh%rr(:,H_mesh%jj(1:3,m))-vv_mesh%rr(:,vv_mesh%jj(1:3,m))))/=0.d0) THEN
             WRITE(*,*) ' BUG Pb de maillage '
             STOP
          END IF
       END DO
    END IF
    IF (ALLOCATED(sigma)) DEALLOCATE(sigma,list_dom_H, list_dom_H_in, list_dom_phi, list_inter_H_phi)    
    IF (ALLOCATED(list_dom_ns)) DEALLOCATE(list_dom_ns)
    !------------------------------------------------------------------------------


    !------------ARRAY ALLOCATION FOR NAVIER STOKES--------------------------------
    IF (type_pb=='nst' .OR. type_pb=='mhd' .OR. type_pb=='mxx') THEN
       ALLOCATE(un_m1   (vv_mesh%np, 6, m_max_c))
       ALLOCATE(un      (vv_mesh%np, 6, m_max_c))
       ALLOCATE(pn_m1   (pp_mesh%np, 2, m_max_c))
       ALLOCATE(pn      (pp_mesh%np, 2, m_max_c))
       ALLOCATE(incpn_m1(pp_mesh%np, 2, m_max_c))
       ALLOCATE(incpn   (pp_mesh%np, 2, m_max_c))
    END IF
    !------------------------------------------------------------------------------


    !------------ARRAY ALLOCATION--------------------------------------------------
    IF (type_pb=='mxw' .OR. type_pb=='mhd' .OR. type_pb=='mxx' ) THEN
       ALLOCATE(Hn1  (H_mesh%np,  6,  m_max_c))
       ALLOCATE(Hn   (H_mesh%np,  6,  m_max_c))
       ALLOCATE(phin1(phi_mesh%np,2,  m_max_c))   
       ALLOCATE(phin (phi_mesh%np,2,  m_max_c))
    END IF
    !------------------------------------------------------------------------------


    !------------ALLOCATION DES VARIABLES DE COUPLAGE------------------------------
    IF (type_pb=='mxw' .OR. type_pb=='mhd' .OR. type_pb=='mxx') THEN
       ALLOCATE(vel(H_mesh%np, 6, m_max_c))
       ! October 7, 2008, JLG
       IF (type_pb=='mhd' .OR. type_pb=='mxx') THEN
          !------------CALCUL DU TABLEAU jj_v_to_H---------------------------------
          !Fev 17, 2009, JLG + CN
          ALLOCATE(jj_v_to_H(MAXVAL(H_mesh%jj(:,1:vv_mesh%me))))
          DO m = 1, vv_mesh%me
             jj_v_to_H(H_mesh%jj(:,m)) = vv_mesh%jj(1:H_mesh%gauss%n_w,m)
          END DO
          !Fev 17, 2009, JLG + CN
          !------------------------------------------------------------------------
       END IF
       ! October 7, 2008, JLG
    END IF
    IF (type_pb=='mhd') THEN
       ALLOCATE(Hn_p2(vv_mesh%np, 6, m_max_c))
       ALLOCATE(Hext(H_mesh%np, 6, m_max_c))
    ELSE
       ALLOCATE(Hn_p2(1, 1, 1))
    END IF
    !------------------------------------------------------------------------------

    !-----------GESTION DU NOMBRE DES SYSTEMES LINEAIRES---------------------------
    IF (type_pb=='nst') THEN
       nb_syst_lin_ns(1) = 4*m_max_c + 1
       nb_syst_lin_ns(2) = 0
       nb_syst_lin_mx = 0
    ELSE IF (type_pb=='mxw' .OR. type_pb=='mxx') THEN
       nb_syst_lin_mx(1) = 2*m_max_c
       nb_syst_lin_mx(2) = 0
       nb_syst_lin_ns = 0
    ELSE IF (type_pb=='mhd') THEN
       nb_syst_lin_ns(1) = 4*m_max_c + 1 + 2*m_max_c
       nb_syst_lin_ns(2) = 0
       nb_syst_lin_mx(1) = nb_syst_lin_ns(1)
       nb_syst_lin_mx(2) = 4*m_max_c + 1
    ELSE
       WRITE(*,*) ' BUG type_pb'
       STOP
    END IF
    !------------------------------------------------------------------------------


    !------------INITIALISATION de NAVIER STOKES-----------------------------------
    IF (type_pb=='nst' .OR. type_pb=='mhd' .OR. type_pb=='mxx') THEN
       IF (irestart_u) THEN
          CALL read_restart_ns(vv_mesh, pp_mesh, time_u, list_mode, un, un_m1, pn, pn_m1, &
               incpn, incpn_m1, file_name)
       ELSE 
          CALL init_up(vv_mesh,pp_mesh,time_u,dt,list_mode,un_m1,un,pn_m1,pn,incpn_m1,incpn)
       END IF
    END IF
    !------------------------------------------------------------------------------
    !Cas dynamo cinematique initialisee a l'aide d'un champ analytique 
    IF (type_pb=='mxw') THEN
       DO i = 1, m_max_c       !Initialization of vel
          vel(:,:,i) = Vexact(list_mode(i), H_mesh)
       END DO
    ENDIF

    !Cas dynamo cinematique initialisee a l'aide d'un champ sortant de NS
    IF (type_pb=='mxx') THEN
       !October 7, 2008, JLG
       !IF (H_mesh%gauss%n_w==vv_mesh%gauss%n_w) THEN
       !   vel_np = vv_mesh%np
       !   vel(1:vv_mesh%np,:,:) = un
       !ELSE !Projection du champ de vitesse P2 sur le maillage P1
       !vel_np = pp_mesh%np
       !CALL project_champ(un, vel, vv_mesh, pp_mesh)
       !END IF
       vel_np = SIZE(jj_v_to_H)
       vel(1:vel_np,:,:) = un(jj_v_to_H,:,:)
       !October 7, 2008, JLG

       IF (dom_H_larger_dom_ns) THEN !We extend vel
          DO i = 1, m_max_c
             DO k= 1, 6 !The user has to code extension_vel 
                vel(vel_np+1:,k,i) = extension_vel(k, H_mesh, list_mode(i), time_u, vel_np+1)
             END DO
          END DO
       END IF

       DEALLOCATE(un_m1)
       DEALLOCATE(un)
       DEALLOCATE(pn_m1)
       DEALLOCATE(pn)
       DEALLOCATE(incpn_m1)
       DEALLOCATE(incpn)
    ENDIF


    !------------INITIALISATION de MAXWELL-----------------------------------
    IF (type_pb=='mxw' .OR. type_pb=='mhd' .OR. type_pb=='mxx') THEN
       IF (irestart_h) THEN
          CALL read_restart_maxwell(H_mesh, phi_mesh, time_h, list_mode, Hn, Hn1, &
               phin, phin1, file_name)
       ELSE
          CALL init_maxwell(H_mesh,phi_mesh,time_h,dt,mu_H_field,mu_phi,list_mode,&
               Hn1,Hn,phin1,phin)
       END IF
    END IF
    !------------------------------------------------------------------------------

    !------------------------------------------------------------------------------
    !On se met sur le temps de celui que l'on a relu si irestart_u /= irestart_h     
    IF (irestart_h .AND. irestart_u) THEN
       time = MAX(time_u,time_h)
    ELSE IF (irestart_h) THEN
       time = time_h
    ELSE IF (irestart_u) THEN
       time = time_u
    ELSE
       !April 17th, 2008, JLG: there was a bug here
       !time = MIN(time_u,time_h)
       time = 0.d0
    END IF
    WRITE(*,*) 'INITIALISATION DONE'
  END SUBROUTINE INIT



    !---------------------------------------------------------------------------
  SUBROUTINE prodmat_navier_stokes(vect_in, vect_out, ndim, i)
    USE subroutine_ns

    IMPLICIT NONE
    INTEGER :: ndim
    REAL(KIND=8), DIMENSION(ndim) :: vect_in, vect_out
    INTEGER :: i

    INTEGER :: type, i_deb, i_fin

    WRITE(*,*) ' Prodmat_navier_stokes not yet programmed '
    STOP

  END SUBROUTINE prodmat_navier_stokes

  !---------------------------------------------------------------------------
  SUBROUTINE prodmat_maxwell_int_by_parts(vect_in, vect_out, ndim, i)
    USE update_maxwell

    IMPLICIT NONE
    INTEGER :: ndim
    REAL(KIND=8), DIMENSION(ndim) :: vect_in, vect_out
    INTEGER :: i

    INTEGER :: type, i_deb, i_fin

    time = 0.d0
    DO type = 1, 6
       i_deb = (type-1)*H_mesh%np+1
       i_fin = i_deb + H_mesh%np -1
       IF (MODULO(type,2)==0 .AND. list_mode(i)==0) THEN
          Hn(:,type,i)   = 0.d0
       ELSE
          Hn(:,type,i)   = vect_in(i_deb:i_fin)
       END IF
    END DO
    DO type = 1, 2
       phin(:,type,i) =0.d0
    END DO

    DO type = 1, 6
       i_deb = 6*H_mesh%np + (type-1)*H_mesh%np+1
       i_fin = i_deb + H_mesh%np -1
       IF (MODULO(type,2)==0 .AND. list_mode(i)==0) THEN
          Hn1(:,type,i)   = 0.d0
       ELSE
          Hn1(:,type,i)   = vect_in(i_deb:i_fin)
       END IF
    END DO
    DO type = 1, 2
       phin1(:,type,i) =0.d0
    END DO

    CALL maxwell_decouple(H_mesh, phi_mesh, interface_H_phi, interface_H_mu, &
         Hn, phin, Hn1, phin1, vel, &
         stab, sigma_field, R_fourier, index_fourier, mu_H_field, mu_phi, &
         time, dt, Rem, list_mode, nb_syst_lin_mx, H_phi_per, data_fichier)

    DO type = 1, 6
       i_deb = (type-1)*H_mesh%np+1
       i_fin = i_deb + H_mesh%np -1
       vect_out(i_deb:i_fin) = Hn(:,type,i) 
    END DO

    DO type = 1, 6
       i_deb = 6*H_mesh%np + (type-1)*H_mesh%np+1
       i_fin = i_deb + H_mesh%np -1
       vect_out(i_deb:i_fin) = Hn1(:,type,i) 
    END DO

  END SUBROUTINE prodmat_maxwell_int_by_parts

  SUBROUTINE prodmat_maxwell(vect_in, vect_out, ndim, i)
    USE update_maxwell

    IMPLICIT NONE
    INTEGER :: ndim
    REAL(KIND=8), DIMENSION(ndim) :: vect_in, vect_out
    INTEGER :: i

    INTEGER :: type, i_deb, i_fin

    time = 0.d0
    DO type = 1, 6
       i_deb = (type-1)*H_mesh%np+1
       i_fin = i_deb + H_mesh%np -1
       IF (MODULO(type,2)==0 .AND. list_mode(i)==0) THEN
          Hn(:,type,i)   = 0.d0
       ELSE
          Hn(:,type,i)   = vect_in(i_deb:i_fin)
       END IF
    END DO
    DO type = 1, 2
       i_deb = 6*H_mesh%np + (type-1)*phi_mesh%np+1
       i_fin = i_deb + phi_mesh%np -1
       IF (MODULO(type,2)==0 .AND. list_mode(i)==0) THEN
          phin(:,type,i) =0.d0
       ELSE
          phin(:,type,i) = vect_in(i_deb:i_fin)
       END IF
     END DO

     DO type = 1, 6
       i_deb = 6*H_mesh%np + 2*phi_mesh%np + (type-1)*H_mesh%np+1
       i_fin = i_deb + H_mesh%np -1
       IF (MODULO(type,2)==0 .AND. list_mode(i)==0) THEN
          Hn1(:,type,i)   = 0.d0
       ELSE
          Hn1(:,type,i)   = vect_in(i_deb:i_fin)
       END IF
    END DO
     DO type = 1, 2
       i_deb = 12*H_mesh%np + 2*phi_mesh%np + (type-1)*phi_mesh%np+1
       i_fin = i_deb + phi_mesh%np -1
       IF (MODULO(type,2)==0 .AND. list_mode(i)==0) THEN
          phin1(:,type,i) =0.d0
       ELSE
          phin1(:,type,i) = vect_in(i_deb:i_fin)
       END IF
     END DO

    CALL maxwell_decouple(H_mesh, phi_mesh, interface_H_phi, interface_H_mu, &
         Hn, phin, Hn1, phin1, vel, &
         stab, sigma_field, R_fourier, index_fourier, mu_H_field, mu_phi, &
         time, dt, Rem, list_mode, nb_syst_lin_mx, H_phi_per, data_fichier)

    DO type = 1, 6
       i_deb = (type-1)*H_mesh%np+1
       i_fin = i_deb + H_mesh%np -1
       vect_out(i_deb:i_fin) = Hn(:,type,i) 
    END DO
    DO type = 1, 2
       i_deb = 6*H_mesh%np + (type-1)*phi_mesh%np+1
       i_fin = i_deb + phi_mesh%np -1
       vect_out(i_deb:i_fin) = phin(:,type,i) 
     END DO
     DO type = 1, 6
       i_deb = 6*H_mesh%np + 2*phi_mesh%np + (type-1)*H_mesh%np+1
       i_fin = i_deb + H_mesh%np -1
       vect_out(i_deb:i_fin) = Hn1(:,type,i) 
    END DO
     DO type = 1, 2
       i_deb = 12*H_mesh%np + 2*phi_mesh%np + (type-1)*phi_mesh%np+1
       i_fin = i_deb + phi_mesh%np -1
       vect_out(i_deb:i_fin) = phin1(:,type,i)
     END DO

  END SUBROUTINE prodmat_maxwell

    !---------------------------------------------------------------------------


END MODULE initialisation
