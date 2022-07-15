MODULE sorting
  PUBLIC :: tri_jlg, sort_real
  PRIVATE
CONTAINS
    SUBROUTINE tri_jlg (a,  a_d, n_a_d)
    !=== sort in ascending order of the integer array  a  and generation
    !=== of the integer array  a_d  whose first  n_a_d  leading entries
    !=== contain different values in ascending order, while all the
    !=== remaining entries are set to zero
    !=== sorting by Shell's method.
    IMPLICIT NONE
    INTEGER, DIMENSION(:), INTENT(INOUT) :: a
    INTEGER, DIMENSION(:), INTENT(OUT)   :: a_d
    INTEGER,               INTENT(OUT)   :: n_a_d
    INTEGER :: n, na, inc, i, j, k, ia
    na = SIZE(a)
    !===sort phase
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
    !===compression phase
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

    SUBROUTINE sort_real (a, b, a_d, b_d, n_a_d)
    !=== sort in ascending order of the integer array  a  and generation
    !=== of the integer array  a_d  whose first  n_a_d  leading entries
    !=== contain different values in ascending order, while all the
    !=== remaining entries are set to zero
    !=== sorting by Shell's method.
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(:), INTENT(INOUT) :: a, b
    REAL(KIND=8), DIMENSION(:), INTENT(OUT)   :: a_d, b_d
    INTEGER,             INTENT(OUT)   :: n_a_d
    INTEGER :: n, na, inc, i, j, k
    REAL(KIND=8) :: aa, bb
    na = SIZE(a)
    !===sort phase
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
          aa = a(i)
          bb = b(i)
          j = i
          DO WHILE (a(j-inc) > aa)
             a(j) = a(j-inc)
             b(j) = b(j-inc)
             j = j - inc
             IF (j <= inc) EXIT
          ENDDO
          a(j) = aa
          b(j) = bb
       ENDDO
    ENDDO
    !===compression phase
    n = 1
    a_d(n) = a(1)
    b_d(n) = b(1)
    DO k = 2, na
       IF (a(k) > a(k-1)) THEN
          n = n + 1
          a_d(n) = a(k)
          b_d(n) = b(k)
       ELSE
          WRITE(*,*) 'We have a problem in the compression phase of tri_jlg', k, k-1
       ENDIF
    ENDDO
    n_a_d = n
    a_d(n_a_d + 1 : na) = 0
    b_d(n_a_d + 1 : na) = 0
  END SUBROUTINE sort_real
  
END MODULE sorting
