MODULE post_routines


  INTERFACE proc_loop
     SUBROUTINE proc_loop(vizdir,procIJK,stride,Gdata,FUNK)
       USE globals, ONLY: flen,varDIM,invprocmap,ax,ay,az,zero,one
       !USE post_routines
       IMPLICIT NONE
       CHARACTER(LEN=flen), INTENT(IN) :: vizdir
       INTEGER, DIMENSION(6) :: procIJK
       INTEGER, DIMENSION(2,2) :: stride
       DOUBLE PRECISION, DIMENSION(:,:,:,:), INTENT(INOUT) :: Gdata
       CHARACTER(LEN=flen), INTENT(IN) :: FUNK
     END SUBROUTINE proc_loop
  END INTERFACE


  INTERFACE block_ops
     SUBROUTINE block_ops(iodata,corner,stride,Gdata,FUNK)
       USE globals, ONLY: varDIM,flen
       !USE post_routines
       IMPLICIT NONE
       REAL(KIND=4), DIMENSION(:,:,:,:), INTENT(IN) :: iodata
       INTEGER, DIMENSION(3), INTENT(IN) :: corner
       INTEGER, DIMENSION(2,2), INTENT(IN) :: stride
       DOUBLE PRECISION, DIMENSION(:,:,:,:), INTENT(INOUT) :: Gdata
       CHARACTER(LEN=flen), INTENT(IN) :: FUNK
     END SUBROUTINE block_ops
  END INTERFACE
  


  INTERFACE SUBSUM3IK
     SUBROUTINE SUBSUM3IK(stride,vars,prof1D,file)
       USE globals, ONLY: flen
       IMPLICIT NONE
       DOUBLE PRECISION, DIMENSION(:,:) :: prof1D
       INTEGER, DIMENSION(:) :: vars
       INTEGER, DIMENSION(2,2) :: stride
       CHARACTER(LEN=flen) :: file
     END SUBROUTINE SUBSUM3IK
  END INTERFACE
  
  
!!$INTERFACE SUBSUM3IJ
!!$   SUBROUTINE SUBSUM3IJ(stride,vars,prof1D)
!!$     IMPLICIT NONE
!!$     DOUBLE PRECISION, DIMENSION(:,:) :: prof1D
!!$     INTEGER, DIMENSION(:) :: vars
!!$     INTEGER, DIMENSION(2,2) :: stride
!!$   END SUBROUTINE SUBSUM3IJ
!!$END INTERFACE
!!$
!!$
!!$INTERFACE SUBSUM3JK
!!$   SUBROUTINE SUBSUM3JK(stride,vars,prof1D)
!!$     IMPLICIT NONE
!!$     DOUBLE PRECISION, DIMENSION(:,:) :: prof1D
!!$     INTEGER, DIMENSION(:) :: vars
!!$     INTEGER, DIMENSION(2,2) :: stride
!!$   END SUBROUTINE SUBSUM3JK
!!$END INTERFACE
!!$
!!$
!!$INTERFACE PLAN_AVE
!!$   SUBROUTINE PLAN_AVE(direction,vars,plane2D)
!!$     IMPLICIT NONE
!!$     CHARACTER(LEN=10) :: direction
!!$     INTEGER, DIMENSION(:) :: vars
!!$     DOUBLE PRECISION, DIMENSION(:,:,:) :: plane2D
!!$   END SUBROUTINE PLAN_AVE
!!$END INTERFACE



END MODULE post_routines
