MODULE globals
  IMPLICIT NONE
  INTEGER, PARAMETER :: flen=90          ! file length

  INTEGER :: nx,ny,nz
  INTEGER :: t1,tf,tloop
  CHARACTER(LEN=flen) :: jobdir

  INTEGER :: ax,ay,az
  INTEGER :: px,py,pz
  INTEGER :: varDIM
  INTEGER :: proc
  INTEGER, DIMENSION(:,:), ALLOCATABLE :: procmap
  INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: invprocmap
  
  DOUBLE PRECISION, PARAMETER :: pi=3.14159265358979323846d0,zero=0.0d0,one=1.0d0,two=2.0d0

END MODULE globals

MODULE sumdata
  IMPLICIT NONE

  INTEGER :: i1p,ifp,j1p,jfp,k1p,kfp
  INTEGER :: i1g,ifg,j1g,jfg,k1g,kfg
  LOGICAL :: fSS3IK = .FALSE.
  LOGICAL :: fSS3JK = .FALSE.
  LOGICAL :: fSS3IJ = .FALSE.
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: out
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: gout
  DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE :: pout,pgrid

  DOUBLE PRECISION :: norm


END MODULE sumdata

