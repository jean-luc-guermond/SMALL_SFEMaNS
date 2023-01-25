MODULE CSR_transpose
  IMPLICIT NONE
CONTAINS
    SUBROUTINE transpose_op(mat,TYPE)
    IMPLICIT NONE
    TYPE(matrice_bloc), INTENT(INOUT):: mat
    CHARACTER(LEN=3),  INTENT(IN)   :: TYPE
    INTEGER, DIMENSION(SIZE(mat%ia)) :: iao
    INTEGER:: i, j, p, next
    IF  (TYPE/='min' .AND. TYPE/='max'.AND. TYPE/='nil') THEN
       WRITE(*,*) ' BUG in tanspose_op'
       STOP
    END IF
    iao = mat%ia
    DO i = 1, SIZE(mat%ia)-1
       DO p = mat%ia(i), mat%ia(i+1)-1 
          j = mat%ja(p)
          next = iao(j)
          iao(j) = next+1
          IF (j.LE.i) CYCLE
          IF (TYPE=='min') THEN
             mat%aa(next) = MIN(mat%aa(p),mat%aa(next))
             mat%aa(p) = mat%aa(next)
          ELSE IF (TYPE=='max') THEN
             mat%aa(next) = MAX(mat%aa(p),mat%aa(next))
             mat%aa(p) = mat%aa(next)
          ELSE IF (TYPE=='nil') THEN
             mat%aa(p) = mat%aa(next)
          END IF
       END DO
    END DO
  END SUBROUTINE transpose_op

  SUBROUTINE TRANSPOSE(mat,mat_tr)
    IMPLICIT NONE
    TYPE(matrice_bloc), INTENT(IN) :: mat
    TYPE(matrice_bloc), INTENT(OUT):: mat_tr
    INTEGER:: i, j, p, n, next
    INTEGER:: pp, tt
    !==Compute lengths of rows of A'.
    n = SIZE(mat%ia)-1
    mat_tr%ia = 0

    DO i = 1, n
       DO p = mat%ia(i), mat%ia(i+1)-1
          j = mat%ja(p) + 1
          mat_tr%ia(j) = mat_tr%ia(j) + 1
       END DO
    END DO
    !===Compute pointers from lengths.
    mat_tr%ia(1) = 1
    DO i = 1, n
       mat_tr%ia(i+1) = mat_tr%ia(i) + mat_tr%ia(i+1)
    END DO
    !===Do the actual copying.
    DO i = 1, n
       DO p = mat%ia(i), mat%ia(i+1)-1
          j = mat%ja(p)
          next = mat_tr%ia(j)
          mat_tr%aa(next) = mat%aa(p)
          mat_tr%ja(next) = i
          mat_tr%ia(j) = next + 1
       END DO
    END DO
    !===Reshift IAO and leave.
    DO i = n, 1, -1
       mat_tr%ia(i+1) = mat_tr%ia(i)
    END DO
    mat_tr%ia(1) = 1

!!$    DO i = 1, n
!!$       DO p = mat%ia(i), mat%ia(i+1)-1
!!$          j = mat%ja(p)
!!$          tt = 0
!!$          DO pp = mat_tr%ia(j), mat_tr%ia(j+1)-1
!!$             IF (mat_tr%ja(pp) == i) THEN
!!$                tt = 1
!!$                IF (ABS(mat_tr%aa(pp)-mat%aa(p))>1.d-15) THEN
!!$                   WRITE(*,*) 'BUG'
!!$                END IF
!!$                EXIT
!!$             END IF
!!$          END DO
!!$          IF (tt==0) THEN
!!$             WRITE(*,*) 'BUG ', tt
!!$          END IF
!!$       END DO
!!$    END DO
  END SUBROUTINE transpose
END MODULE CSR_transpose
