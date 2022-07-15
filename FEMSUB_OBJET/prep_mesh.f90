!===
!Authors: Jean-Luc Guermond, Lugi Quartapelle, Copyright 1994
!Author:  Jean-Luc Guermond, Copyright December 24th 2018
!===
MODULE prep_maill

  IMPLICIT NONE

  PUBLIC :: load_mesh, load_mesh_formatted, load_mesh_free_format, &
       load_dg_mesh_free_format, load_mesh_free_format_ordered, prep_interfaces, &
       load_mesh_free_format_iso
  PRIVATE

CONTAINS

 SUBROUTINE load_mesh_free_format_iso(dir, fil, list_dom, type_fe, mesh, mesh_formatted, edge_stab)
    USE def_type_mesh
    USE chaine_caractere
    USE mod_gauss_points_2d
    USE Dir_nodes
    IMPLICIT NONE
    CHARACTER(len=200),     INTENT(IN) :: dir, fil
    INTEGER, DIMENSION(:), INTENT(IN) :: list_dom
    INTEGER,               INTENT(IN) :: type_fe
    TYPE(mesh_type)                   :: mesh_p1
    TYPE(mesh_type)                   :: mesh
    LOGICAL,               INTENT(IN) :: mesh_formatted !formatted <=> mesh_formatted=.true.
    LOGICAL, OPTIONAL,     INTENT(IN) :: edge_stab
    INTEGER, DIMENSION(3)             :: a_d
    INTEGER, DIMENSION(2)             :: a_ds
    INTEGER, ALLOCATABLE, DIMENSION(:)   :: nouv_nd, nouv_el, nouv_els, &
         ancien_nd, ancien_el, ancien_els
    INTEGER, ALLOCATABLE, DIMENSION(:,:) :: jj_lect, neigh_lect, jjs_lect
    INTEGER, ALLOCATABLE, DIMENSION(:)   :: i_d_lect, sides_lect, neighs_lect
    LOGICAL, ALLOCATABLE, DIMENSION(:)   :: virgin_nd, virgin_el
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:) :: rr_lect
    LOGICAL :: test
    INTEGER :: mnouv, nnouv, i, nb_of_edges, n_a_d, nn
    INTEGER :: n, m, mop, edge, ms, msnouv, neighs1, neighs2
    INTEGER :: d_end, f_end
    INTEGER :: np, nw, me, nws, mes, kd, nwneigh, nw_new, nws_new
    CHARACTER(len=20)                 :: text
    CHARACTER(len=2)                  :: truc

    text = 'Mesh'   
    d_end = last_c_leng (20, text)
    DO n = 1, SIZE(list_dom)   
       d_end = last_c_leng (20, text)
       WRITE(truc,'(i2)') list_dom(n)
       f_end = start_of_string (truc)
       text = text(1:d_end)//'_'//truc(f_end:)
    END DO

    d_end = last_c_leng (20, text)
    IF (type_fe==1) THEN
       text = text(1:d_end)//'_FE_1'
    ELSE IF (type_fe==2) THEN
       text = text(1:d_end)//'_FE_2'
    ELSE IF (type_fe==3) THEN
       text = text(1:d_end)//'_FE_3'
    ELSE
       WRITE(*,*) ' BUG load_mesh_free_format_iso, type_fe not defined'
       STOP
    END IF

    WRITE (*,*) 'Loading mesh-file ...'
    IF (mesh_formatted) THEN
       OPEN(30,FILE=TRIM(ADJUSTL(dir))//'/'//TRIM(ADJUSTL(fil)),FORM='formatted')
    ELSE
       OPEN(30,FILE=TRIM(ADJUSTL(dir))//'/'//TRIM(ADJUSTL(fil)),FORM='unformatted')
    END IF
    OPEN(UNIT=20,FILE=text, FORM='formatted', STATUS='unknown')

    !  READ GRID DATA AND ARRAY ALLOCATION ----------------------------------------
 
    !===Read P1 mesh
    IF (mesh_formatted) THEN
       READ  (30, *)  np,  nw,  me,  nws,  mes
    ELSE
       READ(30)  np,  nw,  me,  nws,  mes
    END IF

    IF (nw==3 .AND. nws==2) THEN ! Decide about space dimension
       kd = 2; nwneigh = 3
    ELSE IF (nw==4 .AND. nws==3) THEN
       kd = 3; nwneigh = 4
    ELSE 
       WRITE(*,*) ' Finite element not yet programmed ', nw, nws
       STOP
    END IF

    ALLOCATE (jj_lect(nw,me),neigh_lect(nwneigh,me),i_d_lect(me))
    ALLOCATE (nouv_nd(np),   ancien_nd(np), virgin_nd(np), &
         nouv_el(0:me), ancien_el(me), virgin_el(me))

    nouv_el = 0

    IF (mesh_formatted) THEN
       DO m = 1, me
          READ(30,*) jj_lect(:,m), neigh_lect(:,m), i_d_lect(m)
       END DO
    ELSE
       READ(30) jj_lect, neigh_lect, i_d_lect
    END IF

    !===Change enumeration
    virgin_nd = .TRUE.
    virgin_el = .TRUE.
    mnouv = 0
    nnouv = 0
    DO m = 1, me
       IF (MINVAL(ABS(list_dom-i_d_lect(m))) /= 0)  CYCLE ! i_d(m) dans la liste
       virgin_el(m) = .FALSE.
       mnouv = mnouv + 1  ! Nouvel element
       nouv_el(m) = mnouv
       ancien_el(mnouv) = m
       DO n = 1, nw; i = jj_lect(n,m)
          IF (virgin_nd(i)) THEN ! Nouveau point
             virgin_nd(i) = .FALSE.
             nnouv = nnouv + 1
             nouv_nd(i) = nnouv
             ancien_nd(nnouv) = i
          END IF
       END DO
    END DO
    mesh%me = mnouv
    mesh_p1%np = nnouv
    ALLOCATE(mesh_p1%jj(nw,mesh%me), mesh_p1%neigh(nwneigh,mesh%me), mesh%i_d(mesh%me))
    DO m = 1, mesh%me
       mesh_p1%jj(:,m) = nouv_nd(jj_lect(:,ancien_el(m)))
       mesh_p1%neigh(:,m) = nouv_el(neigh_lect(:,ancien_el(m)))
       mesh%i_d(m)     = i_d_lect(ancien_el(m))
    END DO
    !===End change enumeration
    
    ALLOCATE (jjs_lect(nws,mes), neighs_lect(mes), sides_lect(mes))
    IF (mesh_formatted) THEN
       DO ms = 1, mes
          READ(30,*) jjs_lect(:,ms), neighs_lect(ms), sides_lect(ms)
       END DO
    ELSE
       READ(30) jjs_lect, neighs_lect, sides_lect
    END IF

    !===Change enumeration
    ALLOCATE(nouv_els(mes), ancien_els(mes))
    DEALLOCATE(virgin_el)
    ALLOCATE(virgin_el(mes))
    virgin_el = .TRUE.
    msnouv = 0
    DO ms = 1, mes
       neighs1 = neighs_lect(ms)
       DO n = 1, nw
          IF (MINVAL(ABS(jjs_lect(:,ms)-jj_lect(n,neighs1))) /= 0) EXIT
       END DO
       neighs2 = neigh_lect(n,neighs1)
       test = .FALSE.
       IF (MINVAL(ABS(list_dom-i_d_lect(neighs1))) == 0) THEN
          test=.TRUE.
       ELSE IF (neighs2 /= 0) THEN
          IF (MINVAL(ABS(list_dom-i_d_lect(neighs2))) == 0) THEN
             test=.TRUE.
             neighs_lect(ms) = neighs2 ! On change de cote
          END IF
       END IF
       IF (.NOT.test) CYCLE
       virgin_el(ms) = .FALSE.
       msnouv = msnouv + 1
       nouv_els(ms) = msnouv
       ancien_els(msnouv) = ms
    END DO
    mesh%mes = msnouv
    ALLOCATE (mesh_p1%jjs(nws,mesh%mes), mesh_p1%neighs(mesh%mes), &
         mesh%sides(mesh%mes))
    DO ms = 1, mesh%mes
       mesh_p1%jjs(:,ms) = nouv_nd(jjs_lect(:,ancien_els(ms)))
       mesh_p1%neighs(ms) = nouv_el(neighs_lect(ancien_els(ms)))
       mesh%sides(ms)  = sides_lect(ancien_els(ms))
    END DO
    !===Check number of edges
    edge = 0
    DO m = 1, mesh%me
       DO n = 1, 3
          mop = mesh_p1%neigh(n,m)
          IF (mop==0) CYCLE !===Edge on boundary
          edge = edge + 1 !===New edge
       END DO
    END DO
    edge = edge/2 !===Number of internal edges
    nb_of_edges = edge + mesh%mes
    IF (edge/=(3*mesh%me - mesh%mes)/2) THEN
       WRITE(*,*) ' BUG, edge/=(3*mesh%me + mesh%mes)/2'
       WRITE(*,*) ' edge ', edge, (3*mesh%me - mesh%mes)/2, MINVAL(mesh_p1%neigh),MAXVAL(mesh_p1%neigh)
       WRITE(*,*) ' mesh%mes ', mesh%mes, ' mesh%me ', mesh%me
       STOP
    END IF
    !===End check number of edges
    !===Change enumeration

    ALLOCATE(rr_lect(kd,np))
    IF (mesh_formatted) THEN
       DO n = 1, np
          READ(30,*) rr_lect(:,n)
       END DO
    ELSE
       READ(30) rr_lect
    END IF
    ALLOCATE(mesh_p1%rr(kd,mesh_p1%np))
    mesh_p1%rr = rr_lect(:,ancien_nd(1:mesh_p1%np))

    !===Make sure indexing in done with the lowest to highest index convention
    jj_lect(:,:mesh%me)    = mesh_p1%jj
    jjs_lect(:,:mesh%mes)  = mesh_p1%jjs
    neigh_lect(:,:mesh%me) = mesh_p1%neigh
    DO m = 1, mesh%me !===loop on the elements
       CALL tri_jlg(jj_lect(:,m),a_d,n_a_d)
       DO n = 1, nw
          i = mesh_p1%jj(n,m)
          DO nn = 1, nw
             IF (a_d(nn) == i) THEN
                mesh_p1%neigh(nn,m)=neigh_lect(n,m)
                EXIT
             END IF
          END DO
       END DO
       mesh_p1%jj(:,m)=a_d
    END DO
    DO ms = 1, mes !===loop on the elements
       CALL tri_jlg(jjs_lect(:,ms),a_ds,n_a_d)
       mesh_p1%jjs(:,ms)=a_ds
    END DO
    !===End make sure indexing in done with the lowest to highest index convention
    
    DEALLOCATE(jj_lect, neigh_lect, i_d_lect)
    DEALLOCATE(jjs_lect, neighs_lect, sides_lect)
    DEALLOCATE(rr_lect, virgin_el, virgin_nd)
    DEALLOCATE(nouv_nd,nouv_el,nouv_els,ancien_nd,ancien_el,ancien_els)
    !===End of P1 grid reading

    !===Prepare actual mesh (works in 2D only)
    IF (kd==3) THEN
       WRITE(*,*) 'k_d==3 not programmed yet'
       STOP
    END IF
    mesh%np = mesh_p1%np + (type_fe-1)*nb_of_edges
    IF (type_fe.GE.3)  mesh%np = mesh%np + (type_fe-2)*mesh%me 
    nw_new = (type_fe+1)*(type_fe+2)/2
    nws_new = type_fe+1
    ALLOCATE(mesh%jj(nw_new,mesh%me))
    ALLOCATE(mesh%neigh(nw,mesh%me))
    ALLOCATE(mesh%jjs(nws_new,mesh%mes))
    ALLOCATE(mesh%neighs(mesh%mes))
    ALLOCATE(mesh%rr(kd,mesh%np))

    !===Create final mesh
    mesh%neigh=mesh_p1%neigh
    mesh%neighs=mesh_p1%neighs
    IF (type_fe==1) THEN
       mesh%jj = mesh_p1%jj
       mesh%jjs = mesh_p1%jjs
       mesh%rr = mesh_p1%rr
    ELSE
       CALL create_iso_grid(mesh_p1%jj, mesh_p1%jjs, mesh_p1%rr, mesh_p1%neigh, mesh_p1%neighs,&
            mesh%jj,  mesh%jjs,  mesh%rr, type_fe)
    END IF
    DEALLOCATE(mesh_p1%jj, mesh_p1%jjs, mesh_p1%rr, mesh_p1%neigh, mesh_p1%neighs)
    !===End Prepare actual mesh
    
    ALLOCATE(mesh%iis(nws_new,mesh%mes))
    CALL dirichlet_nodes(mesh%jjs, SPREAD(1,1,mesh%mes), SPREAD(.TRUE.,1,1), mesh%j_s)
    CALL surf_nodes_i(mesh%jjs, mesh%j_s,  mesh%iis)
    mesh%nps = SIZE(mesh%j_s)

    WRITE (20,*)  'np ', mesh%np, 'me ',  mesh%me, 'mes ',mesh%mes,'nps ', mesh%nps 

    !===Create structure for edge stabilization
    IF (PRESENT(edge_stab)) THEN
       mesh%edge_stab = edge_stab
       IF (mesh%edge_stab) THEN
          CALL prep_interfaces(mesh) ! JLG April 2009
       ELSE
          mesh%mi = 0
       END IF
    ELSE
       mesh%edge_stab = .FALSE.
       mesh%mi = 0
    END IF
    !===End create structure for edge stabilization

    !===Create Gauss points
    IF (kd==2) THEN
       CALL gauss_points_2d(mesh,type_fe)
    ELSE 
       WRITE(*,*) ' Finite element not yet programmed '
       STOP
    END IF
    !===End Create Gauss points
    
    CLOSE(20)
    CLOSE(30)
  END SUBROUTINE load_mesh_free_format_iso

  SUBROUTINE create_iso_grid(jj_in, jjs_in, rr_in, m_op, neigh_el,&
       jj_f,  jjs_f,  rr_f, type_fe)
    !===jj(:, :)    nodes of the  volume_elements of the input grid
    !===jjs(:, :)    nodes of the surface_elements of the input grid
    !===rr(:, :)    cartesian coordinates of the nodes of the input grid
    !===m_op(:,:)   volume element opposite to each node 
    !===neigh_el(:) volume element ajacent to the surface element 
    !===jj_f(:, :)  nodes of the  volume_elements of the output p2 grid
    !===jjs_f(:, :)  nodes of the surface_elements of the output p2 grid
    !===rr_f(:, :)  cartesian coordinates of the nodes of the output p2 grid
    IMPLICIT NONE
    INTEGER,      DIMENSION(:,:), INTENT(INOUT) :: jj_in, jjs_in, m_op
    INTEGER,      DIMENSION(:),   INTENT(IN)    :: neigh_el
    REAL(KIND=8), DIMENSION(:,:), INTENT(IN)    :: rr_in
    INTEGER,                      INTENT(IN)    :: type_fe 
    INTEGER,      DIMENSION(:,:), INTENT(OUT)   :: jj_f, jjs_f
    REAL(KIND=8), DIMENSION(:,:), INTENT(OUT)   :: rr_f
    LOGICAL, DIMENSION(:),   ALLOCATABLE        :: virgin
    INTEGER, DIMENSION(:,:), ALLOCATABLE        :: j_mid, jjs_mid
    INTEGER      :: np, me, mes, nw, nws, kd, n, m, k, l, n_dof
    INTEGER      :: n1, n2, n3, n4, ms, n_start, n_end
    INTEGER      :: n_k1, n_k2, m_op_k, kk, i, mm, ms_bord
    REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: r_mid
    INTEGER,      DIMENSION(type_fe+1)      :: ns3
    REAL(KIND=8), DIMENSION(type_fe+1)      :: scos
    REAL(KIND=8) :: epsilon=1.d-13, dist, d1, d2, s1, s2, s3, shalf, ref, scc
    INTEGER      :: ns, ns1, index, nb_angle, f_dof
    LOGICAL      :: iso
    
    nb_angle = 0
    nw  = SIZE(jj_in,1)   !===nodes in each volume element (3 in 2D)
    me  = SIZE(jj_in,2)   !===number of cells
    kd  = SIZE(rr_in,1)   !===space dimensions
    np  = SIZE(rr_in,2)   !===number of P1 vertices
    mes = SIZE(jjs_in,2)
    nws = SIZE(jjs_in,1)
    f_dof = type_fe-1
    ALLOCATE(virgin(me), j_mid(nw*f_dof,me), jjs_mid(f_dof,mes), r_mid(kd))

    IF (kd == 3) THEN 
       WRITE(*,*) ' CREATE_GRID_Pk: 3D case not programmed yet !'
       STOP
    END IF

    !===GENERATION OF THE Pk GRID
    rr_f(:, 1:np) = rr_in
    jj_f(1:nw, :) = jj_in
    virgin = .TRUE.
    n_dof = np
    DO m = 1, me !===loop on the elements
       DO k = 1, nw !===loop on the nodes (sides) of the element
          m_op_k = m_op(k,m)
          n_k1 = MODULO(k,nw) + 1
          n_k2 = MODULO(k+1,nw) + 1
          n1 = jj_in(n_k1,m)
          n2 = jj_in(n_k2,m)
          IF (n1<n2) THEN !===Go from lowest global index to highest global index
             n_start = n1
             n_end   = n2
          ELSE
             n_start = n2
             n_end   = n1
          END IF
          IF (m_op_k == 0) THEN  !===the side is on the boundary
             iso = .TRUE.
          ELSE
             iso = .FALSE.
          END IF

          IF (iso) THEN
             DO ms = 1, SIZE(jjs_in,2)
                DO ns = 1, SIZE(jjs_in,1) 
                   dist = SQRT(SUM((rr_in(:,n1)-rr_in(:,jjs_in(ns,ms)))**2))   
                   IF (dist.LE.epsilon) THEN
                      ns1 = MODULO(ns,   SIZE(jjs_in,1)) + 1 
                      dist = SQRT(SUM((rr_in(:,n2)-rr_in(:,jjs_in(ns1,ms)))**2))
                      IF (dist.LE.epsilon) THEN
                         ms_bord = ms
                         GO TO 100
                      END IF
                   END IF
                END DO
             END DO
             WRITE(*,*) ' BUG in create_iso_grid'
             !===Algorithm not designed yet for internal interfaces
             STOP

100          index = 1
             ref = SQRT(SUM((rr_in(:,n1)-rr_in(:,n2))**2))
             DO ms = 1, SIZE(jjs_in,2)
                IF (ms==ms_bord) CYCLE  
                d1 = SQRT(SUM((rr_in(:,n1)-rr_in(:,jjs_in(1,ms)))**2))/ref
                d2 = SQRT(SUM((rr_in(:,n1)-rr_in(:,jjs_in(2,ms)))**2))/ref
                IF (d1.LE.epsilon .AND. d2.GT.2*epsilon) THEN
                   scc=SUM((rr_in(:,n1)-rr_in(:,jjs_in(2,ms)))*              (rr_in(:,n1)-rr_in(:,n2)))/ &
                        (SQRT(SUM((rr_in(:,n1)-rr_in(:,jjs_in(2,ms)))**2))*SQRT(SUM((rr_in(:,n1)-rr_in(:,n2))**2)))
                   IF (index.GE.3) THEN
                      IF (scc .GE. MINVAL(scos)) CYCLE
                      index = 2
                   END IF
                   ns3(index)  = jjs_in(2,ms)
                   scos(index) = scc
                   index = index + 1
                ELSE IF (d2.LE.epsilon .AND. d1.GT.2*epsilon) THEN
                   scc=SUM((rr_in(:,n1)-rr_in(:,jjs_in(1,ms)))*              (rr_in(:,n1)-rr_in(:,n2)))/ &
                        (SQRT(SUM((rr_in(:,n1)-rr_in(:,jjs_in(1,ms)))**2))*SQRT(SUM((rr_in(:,n1)-rr_in(:,n2))**2)))
                   IF (index.GE.3) THEN
                      IF (scc .GE. MINVAL(scos)) CYCLE
                      index = 2
                   END IF
                   ns3(index)  = jjs_in(1,ms)
                   scos(index) = scc
                   index = index + 1
                END IF
                d1 = SQRT(SUM((rr_in(:,n2)-rr_in(:,jjs_in(1,ms)))**2))/ref
                d2 = SQRT(SUM((rr_in(:,n2)-rr_in(:,jjs_in(2,ms)))**2))/ref
                IF (d1.LE.epsilon .AND. d2.GT.2*epsilon) THEN
                   scc=SUM((rr_in(:,n2)-rr_in(:,jjs_in(2,ms)))*              (rr_in(:,n2)-rr_in(:,n1)))/ &
                        (SQRT(SUM((rr_in(:,n2)-rr_in(:,jjs_in(2,ms)))**2))*SQRT(SUM((rr_in(:,n2)-rr_in(:,n1))**2)))
                   IF (index.GE.3) THEN
                      IF (scc .GE. MINVAL(scos)) CYCLE
                      index = 2
                   END IF
                   ns3(index)  = jjs_in(2,ms)
                   scos(index) = scc
                   index = index + 1
                ELSE IF (d2.LE.epsilon .AND. d1.GT.2*epsilon) THEN
                   scc=SUM((rr_in(:,n2)-rr_in(:,jjs_in(1,ms)))*              (rr_in(:,n2)-rr_in(:,n1)))/ &
                        (SQRT(SUM((rr_in(:,n2)-rr_in(:,jjs_in(1,ms)))**2))*SQRT(SUM((rr_in(:,n2)-rr_in(:,n1))**2)))
                   IF (index.GE.3) THEN
                      IF (scc .GE. MINVAL(scos)) CYCLE
                      index = 2
                   END IF
                   ns3(index)  = jjs_in(1,ms)
                   scos(index) = scc
                   index = index + 1
                END IF
             END DO

             IF (index.LT.2) THEN
                WRITE(*,*) SIZE(jjs_in,2), ms_bord
                WRITE(*,*) ' BUG: bad index', rr_in(1,jjs_in(1,ms_bord)), rr_in(1,jjs_in(2,ms_bord))
                WRITE(*,*) ' BUG: bad index', rr_in(2,jjs_in(1,ms_bord)), rr_in(2,jjs_in(2,ms_bord))
                STOP
             END IF
             IF (ABS(scos(1)) > ABS(scos(2))) THEN
                n3 = ns3(1)
             ELSE
                n3 = ns3(2)
             END IF
             IF (MINVAL(ABS(scos)) < 0.95) THEN 
                nb_angle = nb_angle + 1
             END IF
             d1 = SQRT(SUM((rr_in(:,n1)-rr_in(:,n3))**2))
             d2 = SQRT(SUM((rr_in(:,n2)-rr_in(:,n3))**2))
             IF (d1 .LT. d2) THEN 
                n4 = n2
                n2 = n1
                n1 = n4
             END IF
             d1 = SQRT(SUM((rr_in(:,n1)-rr_in(:,n3))**2)) 
             d2 = SQRT(SUM((rr_in(:,n1)-rr_in(:,n2))**2))
             s3 = d1 / d2
             s2 = 1.d0
             s1 = 0.d0
             shalf = 0.5d0
             r_mid = rr_in(:,n1)*(shalf - s2)*(shalf - s3)/((s1 - s2)*(s1 - s3)) & 
                  + rr_in(:,n2)*(shalf - s3)*(shalf - s1)/((s2 - s3)*(s2 - s1)) & 
                  + rr_in(:,n3)*(shalf - s1)*(shalf - s2)/((s3 - s1)*(s3 - s2))
             DO l = 1, f_dof 
                n_dof = n_dof + 1 !===New index created
                j_mid((k-1)*f_dof+l, m) = n_dof
                IF (n1<n2) THEN
                   shalf = l/dble(type_fe)
                ELSE
                   shalf = 1-l/dble(type_fe)
                END IF
                rr_f(:, n_dof)  = rr_in(:,n1)*(shalf - s2)*(shalf - s3)/((s1 - s2)*(s1 - s3)) & 
                     + rr_in(:,n2)*(shalf - s3)*(shalf - s1)/((s2 - s3)*(s2 - s1)) & 
                     + rr_in(:,n3)*(shalf - s1)*(shalf - s2)/((s3 - s1)*(s3 - s2))
             END DO
             !===End of iso-grid
             !===Surface elements of the grid are defined later
          ELSE !===the side is internal
             IF (virgin(m_op_k)) THEN !===This side is new
                DO l = 1, f_dof 
                   n_dof = n_dof + 1 !===New index created
                   j_mid((k-1)*f_dof+l, m) = n_dof
                   rr_f(:, n_dof) = rr_in(:, n_start) + l*(rr_in(:, n_end)-rr_in(:, n_start))/type_fe
                END DO
             ELSE !===the side has been already considered
                mm = m_op_k
                DO i = 1, nw
                   IF (m_op(i,mm) == m) THEN
                      kk = i
                      EXIT
                   END IF
                ENDDO
                DO l = 1, f_dof 
                   j_mid((k-1)*f_dof+l, m) = j_mid((kk-1)*f_dof+l, mm) !===New index created
                END DO
             ENDIF
          ENDIF
       ENDDO
       virgin(m) = .FALSE.
       
    ENDDO

    !===connectivity array for iso grid
    DO m = 1, me
       DO n = 1, f_dof*nw
          jj_f(nw+n, m) = j_mid(n,m)
       END DO
       IF (type_fe==3) THEN
          n_dof = n_dof + 1
          jj_f(10,m) = n_dof
          rr_f(:,n_dof) = (rr_in(:,jj_in(1,m))+rr_in(:,jj_in(2,m))+rr_in(:,jj_in(3,m)))/3
       END IF
    END DO

    !==connectivity array the surface elements of the iso grid
    DO ms = 1, mes
       m = neigh_el(ms)
       DO n = 1, kd+1 !===Simplices
          IF (MINVAL(ABS(jj_in(n,m)-jjs_in(:,ms)))/=0) THEN
             kk = n
             EXIT
          END IF
       ENDDO
       DO l = 1, f_dof
          jjs_mid(l,ms) = j_mid((kk-1)*f_dof+l, m) !===New index created
       END DO
    ENDDO
    jjs_f(1:nws, :) = jjs_in
    jjs_f(nws+1:, :) = jjs_mid
    
    DEALLOCATE(virgin,j_mid,jjs_mid,r_mid)
  END SUBROUTINE  create_iso_grid

  SUBROUTINE tri_jlg(a,a_d,n_a_d)
    !  sort in ascending order of the integer array  a  and generation
    !  of the integer array  a_d  whose first  n_a_d  leading entries
    !  contain different values in ascending order, while all the
    !  remaining entries are set to zero
    !  sorting by Shell's method.
    IMPLICIT NONE
    INTEGER, DIMENSION(:), INTENT(INOUT) :: a
    INTEGER, DIMENSION(:), INTENT(OUT)   :: a_d
    INTEGER,               INTENT(OUT)   :: n_a_d
    INTEGER :: n, na, inc, i, j, k, ia
    na = SIZE(a)
    !  sort phase
    IF (na == 0) THEN
       n_a_d = 0
       RETURN
    ENDIF
    inc = 1
    DO WHILE (inc <= na)
       inc = inc * 3
       inc = inc + 1
    ENDDO
    DO WHILE (inc > 1)
       inc = inc/3
       DO i = inc + 1, na
          ia = a(i)
          j = i
          DO WHILE (a(j-inc) > ia)
             a(j) = a(j-inc)
             j = j - inc
             IF (j <= inc) EXIT
          ENDDO
          a(j) = ia
       ENDDO
    ENDDO
    !  compression phase
    n = 1
    a_d(n) = a(1)
    DO k = 2, na
       IF (a(k) > a(k-1)) THEN
          n = n + 1
          a_d(n) = a(k)
       ELSE
          WRITE(*,*) 'We have a problem in the compression phase of tri_jlg', k, k-1
       ENDIF
    ENDDO
    n_a_d = n
    a_d(n_a_d + 1 : na) = 0
  END SUBROUTINE tri_jlg
  
  !------------------------------------------------------------------------------
  SUBROUTINE load_dg_mesh_free_format(dir, fil, list_dom, list_inter, type_fe, &
       mesh, mesh_formatted)
    USE def_type_mesh
    USE chaine_caractere
    USE mod_gauss_points_2d
    USE Dir_nodes
    IMPLICIT NONE
    CHARACTER(len=200),     INTENT(IN) :: dir, fil
    INTEGER, DIMENSION(:), INTENT(IN) :: list_dom, list_inter
    INTEGER,               INTENT(IN) :: type_fe
    TYPE(mesh_type)                   :: mesh
    LOGICAL,               INTENT(IN) :: mesh_formatted !formatted <=> mesh_formatted=.true.
    INTEGER, ALLOCATABLE, DIMENSION(:)   :: nouv_nd, nouv_el, nouv_els
    INTEGER, ALLOCATABLE, DIMENSION(:,:) :: jj_lect, neigh_lect, jjs_lect
    INTEGER, ALLOCATABLE, DIMENSION(:)   :: i_d_lect, sides_lect, neighs_lect, stat
    LOGICAL, ALLOCATABLE, DIMENSION(:)   :: virgin_nd, virgin_ms
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:) :: rr_lect
    LOGICAL :: t1, t2
    INTEGER :: mnouv, nnouv, i, dom
    INTEGER :: n, m, mop, ms, msnouv, neighs1, neighs2
    INTEGER :: d_end, f_end
    INTEGER :: np, nw, me, nws, mes, kd, nwneigh
    CHARACTER(len=40)                 :: text
    CHARACTER(len=2)                  :: truc

    text = 'Mesh'   
    d_end = last_c_leng (20, text)
    DO n = 1, SIZE(list_dom)   
       d_end = last_c_leng (20, text)
       WRITE(truc,'(i2)') list_dom(n)
       f_end = start_of_string (truc)
       text = text(1:d_end)//'_'//truc(f_end:)
    END DO

    d_end = last_c_leng (20, text)
    IF (type_fe==1) THEN
       text = text(1:d_end)//'_FE_1'
    ELSE
       text = text(1:d_end)//'_FE_2'
    END IF

    WRITE (*,*) 'Loading mesh-file ...'
    IF (mesh_formatted) THEN
       OPEN(30,FILE=TRIM(ADJUSTL(dir))//'/'//TRIM(ADJUSTL(fil)),FORM='formatted')
    ELSE
       OPEN(30,FILE=TRIM(ADJUSTL(dir))//'/'//TRIM(ADJUSTL(fil)),FORM='unformatted')
    END IF
    OPEN(UNIT=20,FILE=text, FORM='formatted', STATUS='unknown')

    !  READ GRID DATA AND ARRAY ALLOCATION ----------------------------------------

    IF (type_fe == 2) THEN ! Skip P1 data if needed
       IF (mesh_formatted) THEN
          READ(30,*) np, nw, me, nws, mes
          DO m = 1, me
             READ(30,*)
          END DO
          DO ms = 1, mes
             READ(30,*)
          END DO
          DO n = 1, np
             READ(30,*)
          END DO
       ELSE
          READ(30) np, nw, me, nws, mes
          READ(30)
          READ(30)
          READ(30)
       END IF
    END IF

    IF (mesh_formatted) THEN
       READ(30, *)  np,  nw,  me,  nws,  mes
    ELSE
       READ(30)  np,  nw,  me,  nws,  mes
    END IF

    IF (nw==3 .AND. nws==2) THEN ! Decide about space dimension
       kd = 2; nwneigh = 3
    ELSE IF (nw==6 .AND. nws==3) THEN
       kd = 2; nwneigh = 3
    ELSE IF (nw==4 .AND. nws==3) THEN
       kd = 3; nwneigh = 4
    ELSE IF (nw==10 .AND. nws==6) THEN
       kd = 3; nwneigh = 4
    ELSE 
       WRITE(*,*) ' Finite element not yet programmed '
       WRITE(*,*) nw, nws
       STOP
    END IF

    ALLOCATE (jj_lect(nw,me),neigh_lect(nwneigh,me),i_d_lect(me))
    ALLOCATE (jjs_lect(nws,mes), neighs_lect(mes), sides_lect(mes))
    ALLOCATE(rr_lect(kd,np))

    IF (mesh_formatted) THEN
       DO m = 1, me
          READ(30,*) jj_lect(:,m), neigh_lect(:,m), i_d_lect(m)
       END DO
       DO ms = 1, mes
          READ(30,*) jjs_lect(:,ms), neighs_lect(ms), sides_lect(ms)
       END DO
       DO n = 1, np
          READ(30,*) rr_lect(:,n)
       END DO
    ELSE
       READ(30) jj_lect, neigh_lect, i_d_lect
       READ(30) jjs_lect, neighs_lect, sides_lect
       READ(30) rr_lect
    END IF


    !---Renumerotation------------------------------------------------------------
    ! Identify the status of faces
    ! stat = 1 (interface to be forgotten), stat = 2 (boundary), stat = 3 (real interface)
    ALLOCATE (stat(mes))
    DO ms = 1, mes
       neighs1 = neighs_lect(ms)
       IF (neighs1==0) THEN
          WRITE(*,*) ' BUG in prep_mesh, neighs1=0 '
          STOP
       END IF
       IF (MINVAL(ABS(i_d_lect(neighs1) - list_dom))==0) THEN
          t1 = .TRUE.
       ELSE
          t1 = .FALSE.
       END IF
       DO n = 1, nw
          IF (MINVAL(ABS(jjs_lect(:,ms)-jj_lect(n,neighs1))) /= 0) EXIT ! n not on the interface
       END DO
       neighs2 = neigh_lect(n,neighs1)
       IF (neighs2==0) THEN
          IF (t1) THEN
             stat(ms) = 2  ! face on the boundary of the domain of interest 
          ELSE
             stat(ms) = 1  ! face does not touch the domain of interest 
          END IF
          CYCLE
       END IF
       ! neighs2 /=0 
       IF (MINVAL(ABS(i_d_lect(neighs2) - list_dom))==0) THEN 
          t2 = .TRUE.
       ELSE
          t2 = .FALSE.
       END IF

       IF (t1) THEN
          IF (t2) THEN
             IF (SIZE(list_inter)==0) THEN
                stat(ms) = 1 ! no inteface to treat
             ELSE IF (MINVAL(ABS(sides_lect(ms)-list_inter))==0) THEN
                stat(ms) = 3 ! real interface
             ELSE
                stat(ms) = 1 ! interface to be forgotten 
             END IF
          ELSE
             stat(ms) = 2 ! face at the boundary of the domain of interest
          END IF
       ELSE
          IF (t2) THEN
             stat(ms) = 2 ! face at the boundary of the domain of interest
          ELSE
             stat(ms) = 1 ! on an interface of no interest
          END IF
       END IF

    END DO

    ALLOCATE (nouv_nd(np),  nouv_els(mes), virgin_nd(np), virgin_ms(mes), nouv_el(me))
    nouv_nd = -1000
    virgin_nd = .TRUE.
    virgin_ms = .TRUE.
    mnouv = 0
    msnouv= 0
    nnouv = 0
    DO dom = 1, SIZE(list_dom)
       DO m = 1, me ! Count new nodes from domain: i_d=dom
          IF (list_dom(dom) /= i_d_lect(m))  CYCLE ! i_d(m) pas dans la liste
          mnouv = mnouv + 1  ! Nouvel element
          nouv_el(m) = mnouv
          DO n = 1, nw; i = jj_lect(n,m)
             IF (virgin_nd(i)) THEN ! Nouveau point
                virgin_nd(i) = .FALSE.
                nnouv = nnouv + 1
             END IF
          END DO
       END DO

       DO ms = 1, mes
          IF (stat(ms) /= 1 .AND. stat(ms) /=2 .AND. stat(ms) /=3) THEN
             WRITE(*,*) ' BUG in prep_mesh, stat out of bounds '
             STOP
          END IF

          !I test if ms touches the current domain of interest: i_d = dom
          IF (stat(ms) == 1) CYCLE
          neighs1 = neighs_lect(ms)
          DO n = 1, nw
             IF (MINVAL(ABS(jjs_lect(:,ms)-jj_lect(n,neighs1))) /= 0) EXIT 
             ! exit when n is not on the interface
          END DO
          neighs2 = neigh_lect(n,neighs1)
          IF (i_d_lect(neighs1) /= list_dom(dom)) THEN 
             IF (neighs2 == 0) CYCLE ! face on the boundary and does not touch dom yet
             IF (i_d_lect(neighs2) /= list_dom(dom)) CYCLE ! dom is on neither sides
          END IF
          !End test if ms touches the domain of interest

          IF (virgin_ms(ms)) THEN !New interface
             virgin_ms(ms) = .FALSE.
             msnouv = msnouv + 1
          END IF
          IF (stat(ms) ==3) THEN 
             ! Nodes and sides on the interface are virgin again
             virgin_nd(jjs_lect(:,ms)) = .TRUE. ! interface nodes are virgin again
             virgin_ms(ms) = .TRUE.
          END IF
       END DO
    END DO
    mesh%me  = mnouv
    mesh%np  = nnouv
    mesh%mes = msnouv


    ALLOCATE(mesh%jj(nw,mesh%me), mesh%neigh(nwneigh,mesh%me), mesh%i_d(mesh%me))
    ALLOCATE(mesh%jjs(nws,mesh%mes), mesh%neighs(mesh%mes),  mesh%sides(mesh%mes))
    ALLOCATE(mesh%rr(kd,mesh%np))

    virgin_nd = .TRUE.
    virgin_ms = .TRUE.
    nnouv = 0
    msnouv = 0
    DO dom = 1, SIZE(list_dom)
       !Loop on me and get nouv_el and nouv_nd right
       DO m = 1, me
          IF (list_dom(dom) /=i_d_lect(m))  CYCLE ! i_d(m) pas dans la liste
          DO n = 1, nw; i = jj_lect(n,m)
             IF (virgin_nd(i)) THEN ! Nouveau point
                virgin_nd(i) = .FALSE.
                nnouv = nnouv + 1
                nouv_nd(i) = nnouv
             END IF
          END DO
       END DO

       !Loop again on me and update
       DO m = 1, me
          IF (list_dom(dom) /=i_d_lect(m))  CYCLE ! i_d(m) pas dans la liste
          DO n = 1, nw; i = jj_lect(n,m)
             IF (n .LE. nwneigh) THEN
                mop = neigh_lect(n,m) 
                IF (mop .LE. 0) THEN
                   mesh%neigh(n,nouv_el(m)) = 0
                ELSE IF (MINVAL(ABS(list_dom - i_d_lect(mop))) == 0) THEN
                   mesh%neigh(n,nouv_el(m)) = nouv_el(mop)
                ELSE
                   mesh%neigh(n,nouv_el(m)) = 0
                END IF
             END IF
             mesh%rr(:,nouv_nd(i)) = rr_lect(:,i)
          END DO
          mesh%i_d(nouv_el(m)) = i_d_lect(m)
          mesh%jj(:,nouv_el(m)) = nouv_nd(jj_lect(:,m))
       END DO

       !Loop on mes and get neighs_lect and nouv_els right
       DO ms = 1, mes
          !I test if ms touches the current domain of interest: i_d = dom
          IF (stat(ms) == 1) CYCLE
          neighs1 = neighs_lect(ms)
          DO n = 1, nw
             IF (MINVAL(ABS(jjs_lect(:,ms)-jj_lect(n,neighs1))) /= 0) EXIT 
             ! exit when n is not on the interface
          END DO
          neighs2 = neigh_lect(n,neighs1)
          IF (i_d_lect(neighs1) /= list_dom(dom)) THEN 
             IF (neighs2 == 0) CYCLE ! face on the boundary and does not touch dom yet
             IF (i_d_lect(neighs2) /= list_dom(dom)) CYCLE ! dom is on neither sides
          END IF
          !End test if ms touches the domain of interest

          IF (virgin_ms(ms)) THEN !New interface
             virgin_ms(ms) = .FALSE.
             msnouv = msnouv + 1
             nouv_els(ms) = msnouv
          END IF
          IF (stat(ms) ==3) THEN 
             ! Nodes and sides on the interface are virgin again
             virgin_nd(jjs_lect(:,ms)) = .TRUE. ! interface nodes are virgin again
             virgin_ms(ms) = .TRUE.
          END IF

          ! Swapping problem
          neighs1 = neighs_lect(ms)
          IF ((ABS(i_d_lect(neighs1) - list_dom(dom)))==0) THEN
             t1 = .TRUE.
          ELSE
             t1 = .FALSE.
          END IF
          DO n = 1, nw
             IF (MINVAL(ABS(jjs_lect(:,ms)-jj_lect(n,neighs1))) /= 0) EXIT ! n not on the interface
          END DO
          neighs2 = neigh_lect(n,neighs1)
          IF (neighs2==0) THEN
             CYCLE
          END IF
          ! neighs2 /=0 
          IF ((ABS(i_d_lect(neighs2) - list_dom(dom)))==0) THEN 
             t2 = .TRUE.
          ELSE
             t2 = .FALSE.
          END IF
          IF (.NOT.t1 .AND. t2) THEN
             neighs_lect(ms) = neighs2 !get things right (swap neighs) 
          END IF
       END DO

       !Loop again on mes and update
       DO ms = 1, mes

          !I test if ms touches the current domain of interest: i_d = dom
          IF (stat(ms) == 1) CYCLE
          neighs1 = neighs_lect(ms)
          DO n = 1, nw
             IF (MINVAL(ABS(jjs_lect(:,ms)-jj_lect(n,neighs1))) /= 0) EXIT 
             ! exit when n is not on the interface
          END DO
          neighs2 = neigh_lect(n,neighs1)
          IF (i_d_lect(neighs1) /= list_dom(dom)) THEN 
             IF (neighs2 == 0) CYCLE ! face on the boundary and does not touch dom yet
             IF (i_d_lect(neighs2) /= list_dom(dom)) CYCLE ! dom is on neither sides
          END IF
          !End test if ms touches the domain of interest

          mesh%jjs(:,nouv_els(ms))  = nouv_nd(jjs_lect(:,ms))
          mesh%neighs(nouv_els(ms)) = nouv_el(neighs_lect(ms))
          mesh%sides(nouv_els(ms))  = sides_lect(ms)
          IF (stat(ms)==3) THEN ! side is an interface to be kept
             neighs1 = neighs_lect(ms)
             DO n = 1, nw
                IF (MINVAL(ABS(jjs_lect(:,ms)-jj_lect(n,neighs1))) /= 0) EXIT ! n not on the interface
             END DO
             mesh%neigh(n,nouv_el(neighs1)) = 0
             neighs2 = neigh_lect(n,neighs1)
             DO n = 1, nw
                IF (MINVAL(ABS(jjs_lect(:,ms)-jj_lect(n,neighs2))) /= 0) EXIT ! n not on the interface
             END DO
             mesh%neigh(n,nouv_el(neighs2)) = 0
          END IF
       END DO

    END DO
    !---Fin renumerotation--------------------------------------------------------
    !TEST June 11, 2008
    !write(*,*) 'ok'
    !allocate(const(mesh%me))
    !const = 3.d0
    !do ms = 1, mesh%mes
    !   const(mesh%neighs(ms)) = mesh%i_d(mesh%neighs(ms))
    !end do
    !write(*,*) 'ok'
    !write(*,*) size(mesh%jj,2), size(const)
    !call plot_const_p1_label(mesh%jj, mesh%rr, const, 't.plt')

    !DO ms = 1, mesh%mes
    !   m = mesh%neighs(ms)
    !   !Il faut savoir quelle face est la bonne
    !   DO n = 1, 3
    !      IF (MINVAL(ABS(mesh%jjs(:,ms)-mesh%jj(n,m)))==0) CYCLE
    !      face = n 
    !   END DO
    !   n1 = MODULO(face,3)+1; n2 = MODULO(face+1,3)+1
    !   IF (mesh%jjs(1,ms)==mesh%jj(n1,m) .AND. mesh%jjs(2,ms)==mesh%jj(n2,m)) THEN
    !      nface(1)=n1; nface(2)=n2
    !      ngauss(1)=1; ngauss(2)=2
    !   ELSE IF (mesh%jjs(1,ms)==mesh%jj(n2,m) .AND. mesh%jjs(2,ms)==mesh%jj(n1,m)) THEN
    !      nface(1)=n2; nface(2)=n1
    !      ngauss(1)=2; ngauss(2)=1
    !   ELSE
    !      WRITE(*,*) mesh%np
    !      WRITE(*,*)  mesh%rr(:,mesh%jj(n1,m))
    !      WRITE(*,*)  mesh%rr(:,mesh%jj(n2,m))  
    !      WRITE(*,*)  mesh%rr(:,mesh%jj(face,m))
    !      WRITE(*,*)  mesh%rr(:,mesh%jjs(1,ms))
    !      WRITE(*,*)  mesh%rr(:,mesh%jjs(2,ms))
    !      STOP
    !   END IF
    !END DO
    !stop
    !write(100,*) mesh%neigh
    !write(100,*) mesh%jj
    !write(100,*) mesh%i_d
    !write(100,*) mesh%jjs
    !write(100,*) mesh%neighs
    !write(100,*) mesh%mes
    !write(100,*) mesh%sides
    !stop
    !TEST June 11, 2008


    DEALLOCATE(jj_lect, neigh_lect, i_d_lect)
    DEALLOCATE(jjs_lect, neighs_lect, sides_lect)
    DEALLOCATE(rr_lect, virgin_nd)
    DEALLOCATE(nouv_nd,nouv_el,nouv_els)
    !  END OF GRID READING --------------------------------------------------------

    ALLOCATE(mesh%iis(nws,mesh%mes))
    CALL dirichlet_nodes(mesh%jjs, SPREAD(1,1,mesh%mes), SPREAD(.TRUE.,1,1), mesh%j_s)
    CALL surf_nodes_i(mesh%jjs, mesh%j_s,  mesh%iis)
    mesh%nps = SIZE(mesh%j_s)

    WRITE (20,*)  'np_lect',np,'nw_lect',nw,'nws_lect',nws,'me_lect',me,&
         'mes_lect',mes
    WRITE (20,*)  'np ', mesh%np, 'me ',  mesh%me, &
         'mes ',mesh%mes,'nps ', mesh%nps 

    !   CALL Gauss_gen(mesh%np, mesh%me, mesh%nps, mesh%mes, &
    !                    mesh%jj, mesh%jjs, mesh%rr)
    !===Create Gauss points
    IF (kd==2) THEN
       CALL gauss_points_2d(mesh,type_fe)
    ELSE 
       WRITE(*,*) ' Finite element not yet programmed '
       STOP
    END IF
    !===End Create Gauss points

    CLOSE(20)
    CLOSE(30)

  END SUBROUTINE load_dg_mesh_free_format

  !------------------------------------------------------------------------------

  SUBROUTINE load_mesh_free_format_ordered(dir, fil, list_dom, type_fe, &
       mesh, mesh_formatted, edge_stab)

    USE def_type_mesh
    USE chaine_caractere
    USE mod_gauss_points_2d
    USE Dir_nodes

    IMPLICIT NONE
    CHARACTER(len=200),     INTENT(IN) :: dir, fil
    INTEGER, DIMENSION(:), INTENT(IN) :: list_dom
    INTEGER,               INTENT(IN) :: type_fe
    TYPE(mesh_type)                   :: mesh
    LOGICAL,               INTENT(IN) :: mesh_formatted !formatted <=> mesh_formatted=.true.
    LOGICAL, OPTIONAL,     INTENT(IN) :: edge_stab

    INTEGER, ALLOCATABLE, DIMENSION(:)   :: nouv_nd, nouv_el, nouv_els, &
         ancien_nd, ancien_el, ancien_els
    INTEGER, ALLOCATABLE, DIMENSION(:,:) :: jj_lect, neigh_lect, jjs_lect
    INTEGER, ALLOCATABLE, DIMENSION(:)   :: i_d_lect, sides_lect, neighs_lect
    LOGICAL, ALLOCATABLE, DIMENSION(:)   :: virgin_nd, virgin_el
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:) :: rr_lect

    LOGICAL :: test
    INTEGER :: mnouv, nnouv, i, dom
    INTEGER :: n, m, ms, msnouv, neighs1, neighs2
    INTEGER :: d_end, f_end
    INTEGER :: np, nw, me, nws, mes, kd, nwneigh
    CHARACTER(len=20)                 :: text
    CHARACTER(len=2)                  :: truc

    text = 'Mesh'   
    d_end = last_c_leng (20, text)
    DO n = 1, SIZE(list_dom)   
       d_end = last_c_leng (20, text)
       WRITE(truc,'(i2)') list_dom(n)
       f_end = start_of_string (truc)
       text = text(1:d_end)//'_'//truc(f_end:)
    END DO

    d_end = last_c_leng (20, text)
    IF (type_fe==1) THEN
       text = text(1:d_end)//'_FE_1'
    ELSE
       text = text(1:d_end)//'_FE_2'
    END IF

    WRITE (*,*) 'Loading mesh-file ...'
    IF (mesh_formatted) THEN
       OPEN(30,FILE=TRIM(ADJUSTL(dir))//'/'//TRIM(ADJUSTL(fil)),FORM='formatted')
    ELSE
       OPEN(30,FILE=TRIM(ADJUSTL(dir))//'/'//TRIM(ADJUSTL(fil)),FORM='unformatted')
    END IF
    OPEN(UNIT=20,FILE=text, FORM='formatted', STATUS='unknown')

    !     READ GRID DATA AND ARRAY ALLOCATION ----------------------------------------

    IF (type_fe == 2) THEN    ! Skip P1 data if needed
       IF (mesh_formatted) THEN
          READ(30,*) np, nw, me, nws, mes
       ELSE
          READ(30) np, nw, me, nws, mes
       END IF

       IF (mesh_formatted) THEN
          DO m = 1, me
             READ(30,*)
          END DO
       ELSE
          READ(30)
       END IF

       IF (mesh_formatted) THEN
          DO ms = 1, mes
             READ(30,*)
          END DO
       ELSE
          READ(30)
       END IF

       IF (mesh_formatted) THEN
          DO n = 1, np
             READ(30,*)
          END DO
       ELSE
          READ(30)
       END IF
    END IF

    IF (mesh_formatted) THEN
       READ  (30, *)  np,  nw,  me,  nws,  mes
    ELSE
       READ(30)  np,  nw,  me,  nws,  mes
    END IF

    IF (nw==3 .AND. nws==2) THEN ! Decide about space dimension
       kd = 2; nwneigh = 3
    ELSE IF (nw==6 .AND. nws==3) THEN
       kd = 2; nwneigh = 3
    ELSE IF (nw==4 .AND. nws==3) THEN
       kd = 3; nwneigh = 4
    ELSE IF (nw==10 .AND. nws==6) THEN
       kd = 3; nwneigh = 4
    ELSE 
       WRITE(*,*) ' Finite element not yet programmed '
       WRITE(*,*) kd, nw, nws
       STOP
    END IF

    ALLOCATE (jj_lect(nw,me),neigh_lect(nwneigh,me),i_d_lect(me))
    ALLOCATE (nouv_nd(np),   ancien_nd(np), virgin_nd(np), &
         nouv_el(0:me), ancien_el(me), virgin_el(me))

    nouv_el = 0

    IF (mesh_formatted) THEN
       DO m = 1, me
          READ(30,*) jj_lect(:,m), neigh_lect(:,m), i_d_lect(m)
       END DO
    ELSE
       READ(30) jj_lect, neigh_lect, i_d_lect
    END IF



    !---  Renumerotation------------------------------------------------------------
    virgin_nd = .TRUE.
    virgin_el = .TRUE.
    mnouv = 0
    nnouv = 0
    DO dom = 1, SIZE(list_dom)
       DO m = 1, me
          IF (ABS(list_dom(dom)-i_d_lect(m)) /= 0)  CYCLE ! i_d(m) dans la liste
          virgin_el(m) = .FALSE.
          mnouv = mnouv + 1   ! Nouvel element
          nouv_el(m) = mnouv
          ancien_el(mnouv) = m
          DO n = 1, nw; i = jj_lect(n,m)
             IF (virgin_nd(i)) THEN ! Nouveau point
                virgin_nd(i) = .FALSE.
                nnouv = nnouv + 1
                nouv_nd(i) = nnouv
                ancien_nd(nnouv) = i
             END IF
          END DO
       END DO
    END DO
    mesh%me = mnouv
    mesh%np = nnouv

    ALLOCATE(mesh%jj(nw,mesh%me), mesh%neigh(nwneigh,mesh%me), mesh%i_d(mesh%me))

    DO m = 1, mesh%me
       mesh%jj(:,m) = nouv_nd(jj_lect(:,ancien_el(m)))
       mesh%neigh(:,m) = nouv_el(neigh_lect(:,ancien_el(m)))
       mesh%i_d(m) = i_d_lect(ancien_el(m))
    END DO

    !---  Fin renumerotation--------------------------------------------------------
    ALLOCATE (jjs_lect(nws,mes), neighs_lect(mes), sides_lect(mes))

    IF (mesh_formatted) THEN
       DO ms = 1, mes
          READ(30,*) jjs_lect(:,ms), neighs_lect(ms), sides_lect(ms)
       END DO
    ELSE
       READ(30) jjs_lect, neighs_lect, sides_lect
    END IF

    !---  Renumerotation------------------------------------------------------------
    ALLOCATE(nouv_els(mes), ancien_els(mes))
    DEALLOCATE(virgin_el)
    ALLOCATE(virgin_el(mes))
    virgin_el = .TRUE.
    msnouv = 0
    DO dom = 1, SIZE(list_dom)
       DO ms = 1, mes
          neighs1 = neighs_lect(ms)
          DO n = 1, nw
             IF (MINVAL(ABS(jjs_lect(:,ms)-jj_lect(n,neighs1))) /= 0) EXIT
          END DO
          neighs2 = neigh_lect(n,neighs1)
          test = .FALSE.
          IF (ABS(list_dom(dom)-i_d_lect(neighs1)) == 0) THEN
             test=.TRUE.
          ELSE IF (neighs2 /= 0) THEN
             IF (ABS(list_dom(dom)-i_d_lect(neighs2)) == 0) THEN
                test=.TRUE.
                neighs_lect(ms) = neighs2 ! On change de cote
             END IF
          END IF

          IF (.NOT.test) CYCLE
          !11 June 2007 
          IF (.NOT.virgin_el(ms)) CYCLE 
          !11 June 2007 
          virgin_el(ms) = .FALSE.
          msnouv = msnouv + 1
          nouv_els(ms) = msnouv
          ancien_els(msnouv) = ms
       END DO
    END DO
    mesh%mes = msnouv

    ALLOCATE (mesh%jjs(nws,mesh%mes), mesh%neighs(mesh%mes), &
         mesh%sides(mesh%mes))

    DO ms = 1, mesh%mes
       mesh%jjs(:,ms) = nouv_nd(jjs_lect(:,ancien_els(ms)))
       mesh%neighs(ms) = nouv_el(neighs_lect(ancien_els(ms)))
       mesh%sides(ms)  = sides_lect(ancien_els(ms))
    END DO

    !---  Fin renumerotation--------------------------------------------------------

    ALLOCATE(rr_lect(kd,np))

    IF (mesh_formatted) THEN
       DO n = 1, np
          READ(30,*) rr_lect(:,n)
       END DO
    ELSE
       READ(30) rr_lect
    END IF

    ALLOCATE(mesh%rr(kd,mesh%np))

    mesh%rr = rr_lect(:,ancien_nd(1:mesh%np))

    DEALLOCATE(jj_lect, neigh_lect, i_d_lect)
    DEALLOCATE(jjs_lect, neighs_lect, sides_lect)
    DEALLOCATE(rr_lect, virgin_el, virgin_nd)
    DEALLOCATE(nouv_nd,nouv_el,nouv_els,ancien_nd,ancien_el,ancien_els)
    !     END OF GRID READING --------------------------------------------------------

    ALLOCATE(mesh%iis(nws,mesh%mes))
    CALL dirichlet_nodes(mesh%jjs, SPREAD(1,1,mesh%mes), SPREAD(.TRUE.,1,1), mesh%j_s)
    CALL surf_nodes_i(mesh%jjs, mesh%j_s,  mesh%iis)
    mesh%nps = SIZE(mesh%j_s)

    WRITE (20,*)  'np_lect',np,'nw_lect',nw,'nws_lect',nws,'me_lect',me,&
         'mes_lect',mes
    WRITE (20,*)  'np ', mesh%np, 'me ',  mesh%me, &
         'mes ',mesh%mes,'nps ', mesh%nps 

    !     CALL Gauss_gen(mesh%np, mesh%me, mesh%nps, mesh%mes, &
    !     mesh%jj, mesh%jjs, mesh%rr)
    IF (PRESENT(edge_stab)) THEN
       mesh%edge_stab = edge_stab
       IF (mesh%edge_stab) CALL prep_interfaces(mesh) ! JLG Fevrier 2010
    ELSE
       mesh%edge_stab = .FALSE.
       mesh%mi = 0
    END IF
    !===Create Gauss points
    IF (kd==2) THEN
       CALL gauss_points_2d(mesh,type_fe)
    ELSE 
       WRITE(*,*) ' Finite element not yet programmed '
       STOP
    END IF
    !===End Create Gauss points

    CLOSE(20)
    CLOSE(30)

  END SUBROUTINE load_mesh_free_format_ordered

  !------------------------------------------------------------------------------

  SUBROUTINE load_mesh_free_format(dir, fil, list_dom, type_fe, mesh, mesh_formatted, edge_stab)
    USE def_type_mesh
    USE chaine_caractere
    USE mod_gauss_points_2d
    USE Dir_nodes
    IMPLICIT NONE
    CHARACTER(len=200),     INTENT(IN) :: dir, fil
    INTEGER, DIMENSION(:), INTENT(IN) :: list_dom
    INTEGER,               INTENT(IN) :: type_fe
    TYPE(mesh_type)                   :: mesh
    LOGICAL,               INTENT(IN) :: mesh_formatted !formatted <=> mesh_formatted=.true.
    LOGICAL, OPTIONAL,     INTENT(IN) :: edge_stab
    INTEGER, ALLOCATABLE, DIMENSION(:)   :: nouv_nd, nouv_el, nouv_els, &
         ancien_nd, ancien_el, ancien_els
    INTEGER, ALLOCATABLE, DIMENSION(:,:) :: jj_lect, neigh_lect, jjs_lect
    INTEGER, ALLOCATABLE, DIMENSION(:)   :: i_d_lect, sides_lect, neighs_lect
    LOGICAL, ALLOCATABLE, DIMENSION(:)   :: virgin_nd, virgin_el
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:) :: rr_lect
    LOGICAL :: test
    INTEGER :: mnouv, nnouv, i, edge
    INTEGER :: n, m, mop, ms, msnouv, neighs1, neighs2
    INTEGER :: d_end, f_end
    INTEGER :: np, nw, me, nws, mes, kd, nwneigh
    CHARACTER(len=20)                 :: text
    CHARACTER(len=2)                  :: truc

    text = 'Mesh'   
    d_end = last_c_leng (20, text)
    DO n = 1, SIZE(list_dom)   
       d_end = last_c_leng (20, text)
       WRITE(truc,'(i2)') list_dom(n)
       f_end = start_of_string (truc)
       text = text(1:d_end)//'_'//truc(f_end:)
    END DO

    d_end = last_c_leng (20, text)
    IF (type_fe==1) THEN
       text = text(1:d_end)//'_FE_1'
    ELSE
       text = text(1:d_end)//'_FE_2'
    END IF

    WRITE (*,*) 'Loading mesh-file ...'
    IF (mesh_formatted) THEN
       OPEN(30,FILE=TRIM(ADJUSTL(dir))//'/'//TRIM(ADJUSTL(fil)),FORM='formatted')
    ELSE
       OPEN(30,FILE=TRIM(ADJUSTL(dir))//'/'//TRIM(ADJUSTL(fil)),FORM='unformatted')
    END IF
    OPEN(UNIT=20,FILE=text, FORM='formatted', STATUS='unknown')

    !  READ GRID DATA AND ARRAY ALLOCATION ----------------------------------------

    IF (type_fe == 2) THEN ! Skip P1 data if needed
       IF (mesh_formatted) THEN
          READ(30,*) np, nw, me, nws, mes
       ELSE
          READ(30) np, nw, me, nws, mes
       END IF

       IF (mesh_formatted) THEN
          DO m = 1, me
             READ(30,*)
          END DO
       ELSE
          READ(30)
       END IF

       IF (mesh_formatted) THEN
          DO ms = 1, mes
             READ(30,*)
          END DO
       ELSE
          READ(30)
       END IF

       IF (mesh_formatted) THEN
          DO n = 1, np
             READ(30,*)
          END DO
       ELSE
          READ(30)
       END IF
    END IF

    IF (mesh_formatted) THEN
       READ  (30, *)  np,  nw,  me,  nws,  mes
    ELSE
       READ(30)  np,  nw,  me,  nws,  mes
    END IF

    IF (nw==3 .AND. nws==2) THEN ! Decide about space dimension
       kd = 2; nwneigh = 3
    ELSE IF (nw==6 .AND. nws==3) THEN
       kd = 2; nwneigh = 3
    ELSE IF (nw==4 .AND. nws==3) THEN
       kd = 3; nwneigh = 4
    ELSE IF (nw==10 .AND. nws==6) THEN
       kd = 3; nwneigh = 4
    ELSE 
       WRITE(*,*) ' Finite element not yet programmed ', nw, nws
       STOP
    END IF

    ALLOCATE (jj_lect(nw,me),neigh_lect(nwneigh,me),i_d_lect(me))
    ALLOCATE (nouv_nd(np),   ancien_nd(np), virgin_nd(np), &
         nouv_el(0:me), ancien_el(me), virgin_el(me))

    nouv_el = 0

    IF (mesh_formatted) THEN
       DO m = 1, me
          READ(30,*) jj_lect(:,m), neigh_lect(:,m), i_d_lect(m)
       END DO
    ELSE
       READ(30) jj_lect, neigh_lect, i_d_lect
    END IF

    !---Renumerotation------------------------------------------------------------
    virgin_nd = .TRUE.
    virgin_el = .TRUE.
    mnouv = 0
    nnouv = 0
    DO m = 1, me
       IF (MINVAL(ABS(list_dom-i_d_lect(m))) /= 0)  CYCLE ! i_d(m) dans la liste
       virgin_el(m) = .FALSE.
       mnouv = mnouv + 1  ! Nouvel element
       nouv_el(m) = mnouv
       ancien_el(mnouv) = m
       DO n = 1, nw; i = jj_lect(n,m)
          IF (virgin_nd(i)) THEN ! Nouveau point
             virgin_nd(i) = .FALSE.
             nnouv = nnouv + 1
             nouv_nd(i) = nnouv
             ancien_nd(nnouv) = i
          END IF
       END DO
    END DO

    mesh%me = mnouv
    mesh%np = nnouv

    ALLOCATE(mesh%jj(nw,mesh%me), mesh%neigh(nwneigh,mesh%me), mesh%i_d(mesh%me))

    DO m = 1, mesh%me
       mesh%jj(:,m) = nouv_nd(jj_lect(:,ancien_el(m)))
       mesh%neigh(:,m) = nouv_el(neigh_lect(:,ancien_el(m)))
       mesh%i_d(m) = i_d_lect(ancien_el(m))
    END DO

    !---Fin renumerotation--------------------------------------------------------
    ALLOCATE (jjs_lect(nws,mes), neighs_lect(mes), sides_lect(mes))

    IF (mesh_formatted) THEN
       DO ms = 1, mes
          READ(30,*) jjs_lect(:,ms), neighs_lect(ms), sides_lect(ms)
       END DO
    ELSE
       READ(30) jjs_lect, neighs_lect, sides_lect
    END IF

    !---Renumerotation------------------------------------------------------------
    ALLOCATE(nouv_els(mes), ancien_els(mes))
    DEALLOCATE(virgin_el)
    ALLOCATE(virgin_el(mes))
    virgin_el = .TRUE.
    msnouv = 0
    DO ms = 1, mes
       neighs1 = neighs_lect(ms)
       DO n = 1, nw
          IF (MINVAL(ABS(jjs_lect(:,ms)-jj_lect(n,neighs1))) /= 0) EXIT
       END DO
       neighs2 = neigh_lect(n,neighs1)
       test = .FALSE.
       IF (MINVAL(ABS(list_dom-i_d_lect(neighs1))) == 0) THEN
          test=.TRUE.
       ELSE IF (neighs2 /= 0) THEN
          IF (MINVAL(ABS(list_dom-i_d_lect(neighs2))) == 0) THEN
             test=.TRUE.
             neighs_lect(ms) = neighs2 !===Change side
          END IF
       END IF

       IF (.NOT.test) CYCLE
       virgin_el(ms) = .FALSE.
       msnouv = msnouv + 1
       nouv_els(ms) = msnouv
       ancien_els(msnouv) = ms
    END DO

    mesh%mes = msnouv

    ALLOCATE (mesh%jjs(nws,mesh%mes), mesh%neighs(mesh%mes), &
         mesh%sides(mesh%mes))

    DO ms = 1, mesh%mes
       mesh%jjs(:,ms) = nouv_nd(jjs_lect(:,ancien_els(ms)))
       mesh%neighs(ms) = nouv_el(neighs_lect(ancien_els(ms)))
       mesh%sides(ms)  = sides_lect(ancien_els(ms))
    END DO

    !===Check number of edges
    edge = 0
    DO m = 1, mesh%me
       DO n = 1, 3
          mop = mesh%neigh(n,m)
          IF (mop==0) CYCLE !===Edge on boundary
          edge = edge + 1 !===New edge
       END DO
    END DO
    edge = edge/2 !===Number of internal edges
    IF (edge/=(3*mesh%me - mesh%mes)/2) THEN
       WRITE(*,*) ' BUG, edge/=(3*mesh%me + mesh%mes)/2'
       WRITE(*,*) ' edge ', edge, (3*mesh%me - mesh%mes)/2
       WRITE(*,*) ' mesh%mes ', mesh%mes, ' mesh%me ', mesh%me
       STOP
    END IF
    !===End check number of edges
    !---Fin renumerotation--------------------------------------------------------

    ALLOCATE(rr_lect(kd,np))

    IF (mesh_formatted) THEN
       DO n = 1, np
          READ(30,*) rr_lect(:,n)
       END DO
    ELSE
       READ(30) rr_lect
    END IF

    ALLOCATE(mesh%rr(kd,mesh%np))

    mesh%rr = rr_lect(:,ancien_nd(1:mesh%np))

    DEALLOCATE(jj_lect, neigh_lect, i_d_lect)
    DEALLOCATE(jjs_lect, neighs_lect, sides_lect)
    DEALLOCATE(rr_lect, virgin_el, virgin_nd)
    DEALLOCATE(nouv_nd,nouv_el,nouv_els,ancien_nd,ancien_el,ancien_els)
    !  END OF GRID READING --------------------------------------------------------

    ALLOCATE(mesh%iis(nws,mesh%mes))
    CALL dirichlet_nodes(mesh%jjs, SPREAD(1,1,mesh%mes), SPREAD(.TRUE.,1,1), mesh%j_s)
    CALL surf_nodes_i(mesh%jjs, mesh%j_s,  mesh%iis)
    mesh%nps = SIZE(mesh%j_s)

    WRITE (20,*)  'np_lect',np,'nw_lect',nw,'nws_lect',nws,'me_lect',me,&
         'mes_lect',mes
    WRITE (20,*)  'np ', mesh%np, 'me ',  mesh%me, &
         'mes ',mesh%mes,'nps ', mesh%nps 

    !   CALL Gauss_gen(mesh%np, mesh%me, mesh%nps, mesh%mes, &
    !                    mesh%jj, mesh%jjs, mesh%rr)
    IF (PRESENT(edge_stab)) THEN
       mesh%edge_stab = edge_stab
       IF (mesh%edge_stab) THEN
          CALL prep_interfaces(mesh) ! JLG April 2009
       ELSE
          mesh%mi = 0
       END IF
    ELSE
       mesh%edge_stab = .FALSE.
       mesh%mi = 0
    END IF
   !===Create Gauss points
    IF (kd==2) THEN
       CALL gauss_points_2d(mesh,type_fe)
    ELSE 
       WRITE(*,*) ' Finite element not yet programmed '
       STOP
    END IF
    !===End Create Gauss points

    CLOSE(20)
    CLOSE(30)

  END SUBROUTINE load_mesh_free_format

  !------------------------------------------------------------------------------

  SUBROUTINE load_mesh_formatted(dir, fil, list_dom, type_fe, mesh)

    USE def_type_mesh
    USE chaine_caractere
    USE mod_gauss_points_2d
    USE Dir_nodes

    IMPLICIT NONE
    CHARACTER(len=200),     INTENT(IN) :: dir, fil
    INTEGER, DIMENSION(:), INTENT(IN) :: list_dom
    INTEGER,               INTENT(IN) :: type_fe
    TYPE(mesh_type)                   :: mesh

    INTEGER, ALLOCATABLE, DIMENSION(:)   :: nouv_nd, nouv_el, nouv_els, &
         ancien_nd, ancien_el, ancien_els
    INTEGER, ALLOCATABLE, DIMENSION(:,:) :: jj_lect, neigh_lect, jjs_lect
    INTEGER, ALLOCATABLE, DIMENSION(:)   :: i_d_lect, sides_lect, neighs_lect
    LOGICAL, ALLOCATABLE, DIMENSION(:)   :: virgin_nd, virgin_el
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:) :: rr_lect

    LOGICAL :: test
    INTEGER :: mnouv, nnouv, i
    INTEGER :: n, m, ms, msnouv, neighs1, neighs2
    INTEGER :: np, nw, me, nws, mes, kd, nwneigh

    WRITE (*,*) 'Loading mesh-file ...'

    OPEN(30,FILE=TRIM(ADJUSTL(dir))//'/'//TRIM(ADJUSTL(fil)),FORM='formatted')
    !OPEN(30,FILE=dir(1:d_end)//'/'//fil(1:f_end),FORM='unformatted')
    OPEN(UNIT=20,FILE='error_mesh', FORM='formatted',STATUS='unknown')

    !  READ GRID DATA AND ARRAY ALLOCATION ----------------------------------------

    IF (type_fe == 2) THEN ! Skip P1 data if needed
       READ(30,*) np, nw, me, nws, mes
       !READ(30) np, nw, me, nws, mes
       DO m = 1, me
          READ(30,*)
       END DO
       !READ(30)
       DO ms = 1, mes
          READ(30,*)
       END DO
       !READ(30)
       DO n = 1, np
          READ(30,*)
       END DO
       !READ(30)
    END IF

    READ  (30, *)  np,  nw,  me,  nws,  mes
    !READ  (30)  np,  nw,  me,  nws,  mes

    IF (nw==3 .AND. nws==2) THEN ! Decide about space dimension
       kd = 2; nwneigh = 3
    ELSE IF (nw==6 .AND. nws==3) THEN
       kd = 2; nwneigh = 3
    ELSE IF (nw==4 .AND. nws==3) THEN
       kd = 3; nwneigh = 4
    ELSE IF (nw==10 .AND. nws==6) THEN
       kd = 3; nwneigh = 4
    ELSE 
       WRITE(*,*) ' Finite element not yet programmed '
       STOP
    END IF

    ALLOCATE (jj_lect(nw,me),neigh_lect(nwneigh,me),i_d_lect(me))
    ALLOCATE (nouv_nd(np),   ancien_nd(np), virgin_nd(np), &
         nouv_el(0:me), ancien_el(me), virgin_el(me))

    nouv_el = 0

    DO m = 1, me
       READ(30,*) jj_lect(:,m), neigh_lect(:,m), i_d_lect(m)
    END DO

    !READ(30) jj_lect, neigh_lect, i_d_lect

    !---Renumerotation------------------------------------------------------------
    virgin_nd = .TRUE.
    virgin_el = .TRUE.
    mnouv = 0
    nnouv = 0

    DO m = 1, me
       IF (MINVAL(ABS(list_dom-i_d_lect(m))) /= 0)  CYCLE ! i_d(m) dans la liste
       virgin_el(m) = .FALSE.
       mnouv = mnouv + 1  ! Nouvel element
       nouv_el(m) = mnouv
       ancien_el(mnouv) = m
       DO n = 1, nw; i = jj_lect(n,m)
          IF (virgin_nd(i)) THEN ! Nouveau point
             virgin_nd(i) = .FALSE.
             nnouv = nnouv + 1
             nouv_nd(i) = nnouv
             ancien_nd(nnouv) = i
          END IF
       END DO
    END DO

    mesh%me = mnouv
    mesh%np = nnouv

    ALLOCATE(mesh%jj(nw,mesh%me), mesh%neigh(nwneigh,mesh%me), mesh%i_d(mesh%me))

    DO m = 1, mesh%me
       mesh%jj(:,m) = nouv_nd(jj_lect(:,ancien_el(m)))
       mesh%neigh(:,m) = nouv_el(neigh_lect(:,ancien_el(m)))
       mesh%i_d(m) = i_d_lect(ancien_el(m))
    END DO

    !---Fin renumerotation--------------------------------------------------------
    ALLOCATE (jjs_lect(nws,mes), neighs_lect(mes), sides_lect(mes))

    DO ms = 1, mes
       READ(30,*) jjs_lect(:,ms), neighs_lect(ms), sides_lect(ms)
    END DO

    !READ(30) jjs_lect, neighs_lect, sides_lect

    !---Renumerotation------------------------------------------------------------
    ALLOCATE(nouv_els(mes), ancien_els(mes))
    DEALLOCATE(virgin_el)
    ALLOCATE(virgin_el(mes))
    virgin_el = .TRUE.
    msnouv = 0
    DO ms = 1, mes
       neighs1 = neighs_lect(ms)
       DO n = 1, nw
          IF (MINVAL(ABS(jjs_lect(:,ms)-jj_lect(n,neighs1))) /= 0) EXIT
       END DO
       neighs2 = neigh_lect(n,neighs1)
       test = .FALSE.
       IF (MINVAL(ABS(list_dom-i_d_lect(neighs1))) == 0) THEN
          test=.TRUE.
       ELSE IF (neighs2 /= 0) THEN
          IF (MINVAL(ABS(list_dom-i_d_lect(neighs2))) == 0) THEN
             test=.TRUE.
             neighs_lect(ms) = neighs2 ! On change de cote
          END IF
       END IF
       IF (.NOT.test) CYCLE
       virgin_el(ms) = .FALSE.
       msnouv = msnouv + 1
       nouv_els(ms) = msnouv
       ancien_els(msnouv) = ms
    END DO

    mesh%mes = msnouv

    ALLOCATE (mesh%jjs(nws,mesh%mes), mesh%neighs(mesh%mes), &
         mesh%sides(mesh%mes))

    DO ms = 1, mesh%mes
       mesh%jjs(:,ms) = nouv_nd(jjs_lect(:,ancien_els(ms)))
       mesh%neighs(ms) = nouv_el(neighs_lect(ancien_els(ms)))
       mesh%sides(ms)  = sides_lect(ancien_els(ms))
    END DO

    !---Fin renumerotation--------------------------------------------------------

    ALLOCATE(rr_lect(kd,np))

    DO n = 1, np
       READ(30,*) rr_lect(:,n)
    END DO
    !READ(30) rr_lect

    ALLOCATE(mesh%rr(kd,mesh%np))

    mesh%rr = rr_lect(:,ancien_nd(1:mesh%np))

    DEALLOCATE(jj_lect, neigh_lect, i_d_lect)
    DEALLOCATE(jjs_lect, neighs_lect, sides_lect)
    DEALLOCATE(rr_lect, virgin_el, virgin_nd)   
    DEALLOCATE(nouv_nd,nouv_el,nouv_els,ancien_nd,ancien_el,ancien_els)
    !  END OF GRID READING --------------------------------------------------------

    ALLOCATE(mesh%iis(nws,mesh%mes))
    CALL dirichlet_nodes(mesh%jjs, SPREAD(1,1,mesh%mes), SPREAD(.TRUE.,1,1), mesh%j_s)
    CALL surf_nodes_i(mesh%jjs, mesh%j_s,  mesh%iis)
    mesh%nps = SIZE(mesh%j_s)

    WRITE (20,*)  'np_lect',np,'nw_lect',nw,'nws_lect',nws,'me_lect',me,&
         'mes_lect',mes
    WRITE (20,*)  'np ', mesh%np, 'me ',  mesh%me, &
         'mes ',mesh%mes,'nps ', mesh%nps 

    !   CALL Gauss_gen(mesh%np, mesh%me, mesh%nps, mesh%mes, &
    !                    mesh%jj, mesh%jjs, mesh%rr)
    !===Create Gauss points
    IF (kd==2) THEN
       CALL gauss_points_2d(mesh,type_fe)
    ELSE 
       WRITE(*,*) ' Finite element not yet programmed '
       STOP
    END IF
    !===End Create Gauss points

    CLOSE(20)
    CLOSE(30)

  END SUBROUTINE load_mesh_formatted

  SUBROUTINE load_mesh(dir, fil, list_dom, type_fe, mesh)

    USE def_type_mesh
    USE chaine_caractere
    USE mod_gauss_points_2d
    USE Dir_nodes

    IMPLICIT NONE
    CHARACTER(len=200),     INTENT(IN) :: dir, fil
    INTEGER, DIMENSION(:), INTENT(IN) :: list_dom
    INTEGER,               INTENT(IN) :: type_fe
    TYPE(mesh_type)                   :: mesh

    INTEGER, ALLOCATABLE, DIMENSION(:)   :: nouv_nd, nouv_el, nouv_els, &
         ancien_nd, ancien_el, ancien_els
    INTEGER, ALLOCATABLE, DIMENSION(:,:) :: jj_lect, neigh_lect, jjs_lect
    INTEGER, ALLOCATABLE, DIMENSION(:)   :: i_d_lect, sides_lect, neighs_lect
    LOGICAL, ALLOCATABLE, DIMENSION(:)   :: virgin_nd, virgin_el
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:) :: rr_lect

    LOGICAL :: test
    INTEGER :: mnouv, nnouv, i
    INTEGER :: n, m, ms, msnouv, neighs1, neighs2
    INTEGER :: d_end, f_end
    INTEGER :: np, nw, me, nws, mes, kd, nwneigh
    CHARACTER(len=20)                 :: text
    CHARACTER(len=2)                  :: truc

    text = 'Mesh'   
    d_end = last_c_leng (20, text)
    DO n = 1, SIZE(list_dom)   
       d_end = last_c_leng (20, text)
       WRITE(truc,'(i2)') list_dom(n)
       f_end = start_of_string (truc)
       text = text(1:d_end)//'_'//truc(f_end:)
    END DO

    d_end = last_c_leng (20, text)
    IF (type_fe==1) THEN
       text = text(1:d_end)//'_FE_1'
    ELSE
       text = text(1:d_end)//'_FE_2'
    END IF

    WRITE (*,*) 'Loading mesh-file ...'
    !OPEN(30,FILE=dir(1:d_end)//'/'//fil(1:f_end),FORM='formatted')
    OPEN(30,FILE=TRIM(ADJUSTL(dir))//'/'//TRIM(ADJUSTL(fil)),FORM='unformatted')
    OPEN(UNIT=20,FILE=text, FORM='formatted', STATUS='unknown')

    !  READ GRID DATA AND ARRAY ALLOCATION ----------------------------------------

    IF (type_fe == 2) THEN ! Skip P1 data if needed
       !READ(30,*) np, nw, me, nws, mes
       READ(30) np, nw, me, nws, mes

       !DO m = 1, me
       !   READ(30,*)
       !END DO
       READ(30)

       !DO ms = 1, mes
       !   READ(30,*)
       !END DO

       ! DO m = 1, mes 
       !     READ(30)
       ! END DO
       READ(30)

       !DO n = 1, np
       !   READ(30,*)
       !END DO
       READ(30)
    END IF

    !READ  (30, *)  np,  nw,  me,  nws,  mes
    READ(30)  np,  nw,  me,  nws,  mes

    IF (nw==3 .AND. nws==2) THEN ! Decide about space dimension
       kd = 2; nwneigh = 3
    ELSE IF (nw==6 .AND. nws==3) THEN
       kd = 2; nwneigh = 3
    ELSE IF (nw==4 .AND. nws==3) THEN
       kd = 3; nwneigh = 4
    ELSE IF (nw==10 .AND. nws==6) THEN
       kd = 3; nwneigh = 4
    ELSE 
       WRITE(*,*) ' Finite element not yet programmed '
       STOP
    END IF

    ALLOCATE (jj_lect(nw,me),neigh_lect(nwneigh,me),i_d_lect(me))
    ALLOCATE (nouv_nd(np),   ancien_nd(np), virgin_nd(np), &
         nouv_el(0:me), ancien_el(me), virgin_el(me))

    nouv_el = 0

    !DO m = 1, me
    !   READ(30,*) jj_lect(:,m), neigh_lect(:,m), i_d_lect(m)
    !END DO

    READ(30) jj_lect, neigh_lect, i_d_lect

    !---Renumerotation------------------------------------------------------------
    virgin_nd = .TRUE.
    virgin_el = .TRUE.
    mnouv = 0
    nnouv = 0

    DO m = 1, me
       IF (MINVAL(ABS(list_dom-i_d_lect(m))) /= 0)  CYCLE ! i_d(m) dans la liste
       virgin_el(m) = .FALSE.
       mnouv = mnouv + 1  ! Nouvel element
       nouv_el(m) = mnouv
       ancien_el(mnouv) = m
       DO n = 1, nw; i = jj_lect(n,m)
          IF (virgin_nd(i)) THEN ! Nouveau point
             virgin_nd(i) = .FALSE.
             nnouv = nnouv + 1
             nouv_nd(i) = nnouv
             ancien_nd(nnouv) = i
          END IF
       END DO
    END DO

    mesh%me = mnouv
    mesh%np = nnouv

    ALLOCATE(mesh%jj(nw,mesh%me), mesh%neigh(nwneigh,mesh%me), mesh%i_d(mesh%me))

    DO m = 1, mesh%me
       mesh%jj(:,m) = nouv_nd(jj_lect(:,ancien_el(m)))
       mesh%neigh(:,m) = nouv_el(neigh_lect(:,ancien_el(m)))
       mesh%i_d(m) = i_d_lect(ancien_el(m))
    END DO

    !---Fin renumerotation--------------------------------------------------------
    ALLOCATE (jjs_lect(nws,mes), neighs_lect(mes), sides_lect(mes))

    !DO ms = 1, mes
    !   READ(30,*) jjs_lect(:,ms), neighs_lect(ms), sides_lect(ms)
    !END DO


    !TEST
    !DO ms = 1, mes
    !READ(30) jjs_lect(:,ms), neighs_lect(ms), sides_lect(ms)
    !END DO
    !TEST
    READ(30) jjs_lect, neighs_lect, sides_lect

    !---Renumerotation------------------------------------------------------------
    ALLOCATE(nouv_els(mes), ancien_els(mes))
    DEALLOCATE(virgin_el)
    ALLOCATE(virgin_el(mes))
    virgin_el = .TRUE.
    msnouv = 0
    DO ms = 1, mes

       neighs1 = neighs_lect(ms)
       DO n = 1, nw
          IF (MINVAL(ABS(jjs_lect(:,ms)-jj_lect(n,neighs1))) /= 0) EXIT
       END DO
       neighs2 = neigh_lect(n,neighs1)
       test = .FALSE.
       IF (MINVAL(ABS(list_dom-i_d_lect(neighs1))) == 0) THEN
          test=.TRUE.
       ELSE IF (neighs2 /= 0) THEN
          IF (MINVAL(ABS(list_dom-i_d_lect(neighs2))) == 0) THEN
             test=.TRUE.
             neighs_lect(ms) = neighs2 ! On change de cote
          END IF
       END IF

       IF (.NOT.test) CYCLE
       virgin_el(ms) = .FALSE.
       msnouv = msnouv + 1
       nouv_els(ms) = msnouv
       ancien_els(msnouv) = ms
    END DO

    mesh%mes = msnouv

    ALLOCATE (mesh%jjs(nws,mesh%mes), mesh%neighs(mesh%mes), &
         mesh%sides(mesh%mes))

    DO ms = 1, mesh%mes
       mesh%jjs(:,ms) = nouv_nd(jjs_lect(:,ancien_els(ms)))
       mesh%neighs(ms) = nouv_el(neighs_lect(ancien_els(ms)))
       mesh%sides(ms)  = sides_lect(ancien_els(ms))
    END DO

    !---Fin renumerotation--------------------------------------------------------

    ALLOCATE(rr_lect(kd,np))

    !DO n = 1, np
    !   READ(30,*) rr_lect(:,n)
    !END DO

    READ(30) rr_lect

    ALLOCATE(mesh%rr(kd,mesh%np))

    mesh%rr = rr_lect(:,ancien_nd(1:mesh%np))

    DEALLOCATE(jj_lect, neigh_lect, i_d_lect)
    DEALLOCATE(jjs_lect, neighs_lect, sides_lect)
    DEALLOCATE(rr_lect, virgin_el, virgin_nd)
    DEALLOCATE(nouv_nd,nouv_el,nouv_els,ancien_nd,ancien_el,ancien_els)
    !  END OF GRID READING --------------------------------------------------------

    ALLOCATE(mesh%iis(nws,mesh%mes))
    CALL dirichlet_nodes(mesh%jjs, SPREAD(1,1,mesh%mes), SPREAD(.TRUE.,1,1), mesh%j_s)
    CALL surf_nodes_i(mesh%jjs, mesh%j_s,  mesh%iis)
    mesh%nps = SIZE(mesh%j_s)

    WRITE (20,*)  'np_lect',np,'nw_lect',nw,'nws_lect',nws,'me_lect',me,&
         'mes_lect',mes
    WRITE (20,*)  'np ', mesh%np, 'me ',  mesh%me, &
         'mes ',mesh%mes,'nps ', mesh%nps 

    !===Create Gauss points
    IF (kd==2) THEN
       CALL gauss_points_2d(mesh,type_fe)
    ELSE 
       WRITE(*,*) ' Finite element not yet programmed '
       STOP
    END IF
    !===End Create Gauss points

    CLOSE(20)
    CLOSE(30)

  END SUBROUTINE load_mesh


  SUBROUTINE surf_nodes_i(jjs, j_s,  iis)

    !  generation of the surface element connectivity matrix  iis
    !  based on the surface node numbering, starting from the
    !  connectivity matrix  jjs  of the surface elements according
    !  to the volume node numbering, and from the array  j_s  of
    !  the boundary nodes according to the volume node numbering

    IMPLICIT NONE

    INTEGER, DIMENSION(:,:), INTENT(IN)  :: jjs
    INTEGER, DIMENSION(:),   INTENT(IN)  :: j_s
    INTEGER, DIMENSION(:,:), INTENT(OUT) :: iis

    INTEGER :: ms, ls, j, i

    DO ms = 1, SIZE(jjs,2)
       DO ls = 1, SIZE(jjs,1)
          j = jjs(ls,ms)
          DO i = 1, SIZE(j_s)
             IF ( j == j_s(i) )  iis(ls,ms) = i
          ENDDO
       ENDDO
    ENDDO

  END SUBROUTINE surf_nodes_i

  SUBROUTINE prep_interfaces(mesh)
    USE def_type_mesh
    IMPLICIT NONE
    TYPE(mesh_type)                   :: mesh
    LOGICAL, DIMENSION(mesh%me) :: virgin
    INTEGER :: m, mop, nw, me, l, lop, n, n1, n2, nmin, nmax, edge, nt,nws, f_dof, start, endf
    nw  = SIZE(mesh%jj,1)
    nws = SIZE(mesh%jjs,1)
    me  = mesh%me
    IF (SIZE(mesh%rr,1)==2) THEN
       nt=3
       f_dof= nws-2
    ELSE
       WRITE(*,*) ' BUG: prep_interfaces, 3D not programmed yet '
       STOP
       nt=4
    END IF

    virgin = .TRUE.
    edge = 0
    DO m = 1, me
       virgin(m) = .FALSE.
       DO n = 1, nt
          mop = mesh%neigh(n,m)
          IF (mop==0) CYCLE !Edge on boundary
          IF (.NOT.virgin(mop)) CYCLE !Edge already done
          edge = edge + 1 !New edge
       END DO
    END DO
    IF (SIZE(mesh%rr,1)==2) THEN
       IF (edge/=(3*mesh%me - mesh%mes)/2) THEN
          WRITE(*,*) ' BUG in prep_interfaces, edge/=(3*mesh%me + mesh%mes)/2'
          WRITE(*,*) ' edge ', edge, (3*mesh%me - mesh%mes)/2
          WRITE(*,*) ' mesh%mes ', mesh%mes, ' mesh%me ',mesh%me
          STOP
       END IF
    END IF

    mesh%mi = edge
    ALLOCATE(mesh%jji(nw,2,edge), mesh%jjsi(nws,edge),mesh%neighi(2,edge) )

    edge = 0
    virgin = .TRUE.
    DO m = 1, me
       virgin(m) = .FALSE.
       DO n = 1, nt
          mop = mesh%neigh(n,m)
          IF (mop==0) CYCLE !Edge on boundary
          IF (.NOT.virgin(mop)) CYCLE !Edge already done
          edge = edge + 1 !New edge
          n1 = MODULO(n,nt) + 1
          n2 = MODULO(n+1,nt) + 1 ! Works in 2D only
          nmin = MIN(n1,n2)
          nmax = MAX(n1,n2)
          mesh%jjsi(1,edge) = mesh%jj(nmin,m)
          mesh%jjsi(2,edge) = mesh%jj(nmax,m)
          IF (f_dof>0) THEN
             start = 3+(5-(n1+n2))*f_dof+1
             endf = start+f_dof-1
             mesh%jjsi(3:,edge) = mesh%jj(start:endf,m)
          END IF
          IF (m<mop) THEN 
             l = 1 ! Side 1 on smallest cell index
             lop = 2
             mesh%neighi(1,edge) = m
             mesh%neighi(2,edge) = mop
          ELSE
             l = 2 ! Side 2 on largest cell index
             lop = 1
             mesh%neighi(1,edge) = mop
             mesh%neighi(2,edge) = m
          END IF
          mesh%jji(:,l,edge) = mesh%jj(:,m)
          mesh%jji(:,lop,edge) = mesh%jj(:,mop)
       END DO
    END DO

  END SUBROUTINE prep_interfaces

END MODULE prep_maill




