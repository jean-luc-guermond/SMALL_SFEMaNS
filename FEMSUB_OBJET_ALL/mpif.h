!
!Authors Jean-Luc Guermond, Raphael Laguerre, Copyrights 2005
!
EXTERNAL MPI_INIT, MPI_COMM_RANK,MPI_COMM_SIZE, MPI_SENDRECV, MPI_ALLREDUCE, MPI_ALLTOALL
 
INTEGER :: MPI_DOUBLE_PRECISION=1, MPI_SUM, MPI_COMM_WORLD
INTEGER, PARAMETER :: MPI_STATUS_SIZE=1

INTERFACE
   FUNCTION MPI_WTIME() RESULT(time)
   IMPLICIT NONE
   REAL(KIND=8) :: time
   END FUNCTION MPI_WTIME
END INTERFACE
