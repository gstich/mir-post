MODULE block

  INTERFACE block_ops
     SUBROUTINE block_ops(iodata,corner,stride,Gdata,FUNK)
       USE globals, ONLY: varDIM,flen
       IMPLICIT NONE
       REAL(KIND=4), DIMENSION(:,:,:,:), INTENT(IN) :: iodata
       INTEGER, DIMENSION(3), INTENT(IN) :: corner
       INTEGER, DIMENSION(2,2), INTENT(IN) :: stride
       DOUBLE PRECISION, DIMENSION(:,:,:,:), INTENT(INOUT) :: Gdata
       CHARACTER(LEN=flen), INTENT(IN) :: FUNK
     END SUBROUTINE block_ops
  END INTERFACE

END MODULE block


MODULE post_routines

  INTERFACE proc_loop
     SUBROUTINE proc_loop(vizdir,procIJK,stride,Gdata,FUNK)
       USE globals, ONLY: flen,varDIM,invprocmap,ax,ay,az,zero,one
       USE block
       IMPLICIT NONE
       CHARACTER(LEN=flen), INTENT(IN) :: vizdir
       INTEGER, DIMENSION(6) :: procIJK
       INTEGER, DIMENSION(2,2) :: stride
       DOUBLE PRECISION, DIMENSION(:,:,:,:), INTENT(INOUT) :: Gdata
       CHARACTER(LEN=flen), INTENT(IN) :: FUNK
     END SUBROUTINE proc_loop
  END INTERFACE

END MODULE post_routines




MODULE operators

  INTERFACE SUBSUM3IK
     SUBROUTINE SUBSUM3IK(stride,prof1D,file)
       USE globals, ONLY: flen
       USE post_routines
       IMPLICIT NONE
       DOUBLE PRECISION, DIMENSION(:,:) :: prof1D
       INTEGER, DIMENSION(2,2) :: stride
       CHARACTER(LEN=flen) :: file
     END SUBROUTINE SUBSUM3IK
  END INTERFACE


INTERFACE SUBSUM3JK
   SUBROUTINE SUBSUM3JK(stride,prof1D,file)
     USE globals, ONLY: flen
     USE post_routines
     IMPLICIT NONE
     DOUBLE PRECISION, DIMENSION(:,:) :: prof1D
     INTEGER, DIMENSION(2,2) :: stride
     CHARACTER(LEN=flen) :: file
   END SUBROUTINE SUBSUM3JK
END INTERFACE


INTERFACE SUBSUM3IJ
   SUBROUTINE SUBSUM3IJ(stride,prof1D,file)
     USE globals, ONLY: flen
     USE post_routines
     IMPLICIT NONE
     DOUBLE PRECISION, DIMENSION(:,:) :: prof1D
     INTEGER, DIMENSION(2,2) :: stride
     CHARACTER(LEN=flen) :: file
   END SUBROUTINE SUBSUM3IJ
END INTERFACE


  INTERFACE PLANE_IK
     SUBROUTINE PLANE_IK(stride,prof2D,file)
       USE globals, ONLY: flen
       USE post_routines
       IMPLICIT NONE
       DOUBLE PRECISION, DIMENSION(:,:,:) :: prof2D
       INTEGER, DIMENSION(2,2) :: stride
       CHARACTER(LEN=flen) :: file
     END SUBROUTINE PLANE_IK
  END INTERFACE

  INTERFACE PLANE_JK
     SUBROUTINE PLANE_JK(stride,prof2D,file)
       USE globals, ONLY: flen
       USE post_routines
       IMPLICIT NONE
       DOUBLE PRECISION, DIMENSION(:,:,:) :: prof2D
       INTEGER, DIMENSION(2,2) :: stride
       CHARACTER(LEN=flen) :: file
     END SUBROUTINE PLANE_JK
  END INTERFACE

  INTERFACE PLANE_IJ
     SUBROUTINE PLANE_IJ(stride,prof2D,file)
       USE globals, ONLY: flen
       USE post_routines
       IMPLICIT NONE
       DOUBLE PRECISION, DIMENSION(:,:,:) :: prof2D
       INTEGER, DIMENSION(2,2) :: stride
       CHARACTER(LEN=flen) :: file
     END SUBROUTINE PLANE_IJ
  END INTERFACE

  INTERFACE BLOCK
     SUBROUTINE BLOCK(dataB,file)
       USE globals
       USE post_routines
       IMPLICIT NONE
       DOUBLE PRECISION, DIMENSION(:,:,:,:), INTENT(INOUT) :: dataB
       CHARACTER(LEN=flen) :: file
     END SUBROUTINE BLOCK
  END INTERFACE




END MODULE operators
