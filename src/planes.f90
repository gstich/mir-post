SUBROUTINE PLANE_IK(stride,prof2D,file)
  USE globals
  USE post_routines
  IMPLICIT NONE
  DOUBLE PRECISION, DIMENSION(:,:,:) :: prof2D
  INTEGER, DIMENSION(3,2) :: stride
  CHARACTER(LEN=flen) :: file

  DOUBLE PRECISION, DIMENSION(:,:,:,:), ALLOCATABLE :: Gdata
  INTEGER, DIMENSION(6) :: procIJK
  CHARACTER(LEN=flen) :: FUNK
  DOUBLE PRECISION :: norm,tmp
  INTEGER :: i1p,ifp,j1p,jfp,k1p,kfp
  INTEGER :: i1g,ifg,j1g,jfg,k1g,kfg
  INTEGER :: i


  !! Set the index box size for the averaging kernal
  i1g = 1
  ifg = nx  

  j1g = stride(2,1)
  jfg = stride(2,2)

  k1g = 1
  kfg = nz

  !! Set the case for Summation rule
  FUNK = 'PLANE_IK'
  
  !! Allocate and initialize
  ALLOCATE(Gdata(nx,1,nz,DIM))
  Gdata = 0.0D0

  !! Set the processor bounds for the data
  tmp = dble(i1g)/dble(ax)
  i1p = CEILING(tmp) - 1
  tmp = dble(ifg)/dble(ax)
  ifp = CEILING(tmp) - 1

  tmp = dble(j1g)/dble(ay)
  j1p = CEILING(tmp) - 1
  tmp = dble(jfg)/dble(ay)
  jfp = CEILING(tmp) - 1

  tmp = dble(k1g)/dble(az)
  k1p = CEILING(tmp) - 1
  tmp = dble(kfg)/dble(az)
  kfp = CEILING(tmp) - 1

  procIJK(1) = i1p 
  procIJK(2) = ifp
  procIJK(3) = j1p
  procIJK(4) = jfp
  procIJK(5) = k1p
  procIJK(6) = kfp 

  !! Loop over all the processors and perform the operations
  CALL proc_loop(file,procIJK,stride,Gdata,FUNK)

  !! Normalize by the slice area
  norm = dble( jfg - j1g + 1)
  Gdata = Gdata / norm

  !! Send back to main in the proper format
  prof2D = Gdata(:,1,:,:)
  
     

END SUBROUTINE PLANE_IK


SUBROUTINE PLANE_JK(stride,prof2D,file)
  USE globals
  USE post_routines
  IMPLICIT NONE
  DOUBLE PRECISION, DIMENSION(:,:,:) :: prof2D
  INTEGER, DIMENSION(3,2), INTENT(IN) :: stride
  CHARACTER(LEN=flen) :: file

  DOUBLE PRECISION, DIMENSION(:,:,:,:), ALLOCATABLE :: Gdata
  INTEGER, DIMENSION(6) :: procIJK
  CHARACTER(LEN=flen) :: FUNK
  DOUBLE PRECISION :: norm,tmp
  INTEGER :: i1p,ifp,j1p,jfp,k1p,kfp
  INTEGER :: i1g,ifg,j1g,jfg,k1g,kfg
  INTEGER :: i


  !! Set the index box size for the averaging kernal
  i1g = stride(1,1)
  ifg = stride(1,2)
  
  j1g = 1
  jfg = ny

  k1g = 1
  kfg = nz

  !! Set the case for Summation rule
  FUNK = 'PLANE_JK'

  !! Allocate and initialize
  ALLOCATE(Gdata(1,ny,nz,DIM))
  Gdata = 0.0D0  

  !! Set the processor bounds for the data
  tmp = dble(i1g)/dble(ax)
  i1p = CEILING(tmp) - 1
  tmp = dble(ifg)/dble(ax)
  ifp = CEILING(tmp) - 1

  tmp = dble(j1g)/dble(ay)
  j1p = CEILING(tmp) - 1
  tmp = dble(jfg)/dble(ay)
  jfp = CEILING(tmp) - 1

  tmp = dble(k1g)/dble(az)
  k1p = CEILING(tmp) - 1
  tmp = dble(kfg)/dble(az)
  kfp = CEILING(tmp) - 1

  procIJK(1) = i1p 
  procIJK(2) = ifp
  procIJK(3) = j1p
  procIJK(4) = jfp
  procIJK(5) = k1p
  procIJK(6) = kfp 

  !! Loop over all the processors and perform the operations
  CALL proc_loop(file,procIJK,stride,Gdata,FUNK)  

  !! Normalize by the slice area
  norm = dble( ifg - i1g + 1) 
  Gdata = Gdata / norm

  !! Send back to main in the proper format
  prof2D = Gdata(1,:,:,:)
  


END SUBROUTINE PLANE_JK


SUBROUTINE PLANE_IJ(stride,prof2D,file)
  USE globals
  USE post_routines
  IMPLICIT NONE
  DOUBLE PRECISION, DIMENSION(:,:,:) :: prof2D
  INTEGER, DIMENSION(3,2), INTENT(IN) :: stride
  CHARACTER(LEN=flen) :: file

  DOUBLE PRECISION, DIMENSION(:,:,:,:), ALLOCATABLE :: Gdata
  INTEGER, DIMENSION(6) :: procIJK
  CHARACTER(LEN=flen) :: FUNK
  DOUBLE PRECISION :: norm,tmp
  INTEGER :: i1p,ifp,j1p,jfp,k1p,kfp
  INTEGER :: i1g,ifg,j1g,jfg,k1g,kfg
  INTEGER :: i

  !! Set the index box size for the averaging kernal
  i1g = 1
  ifg = nx
  j1g = 1
  jfg = ny
  
  k1g = stride(3,1)
  kfg = stride(3,2)

  !! Set the case for Summation rule
  FUNK = 'PLANE_IJ'

  !! Allocate and initialize
  ALLOCATE(Gdata(nx,ny,1,DIM))
  Gdata = 0.0D0  

  !! Set the processor bounds for the data
  tmp = dble(i1g)/dble(ax)
  i1p = CEILING(tmp) - 1
  tmp = dble(ifg)/dble(ax)
  ifp = CEILING(tmp) - 1

  tmp = dble(j1g)/dble(ay)
  j1p = CEILING(tmp) - 1
  tmp = dble(jfg)/dble(ay)
  jfp = CEILING(tmp) - 1

  tmp = dble(k1g)/dble(az)
  k1p = CEILING(tmp) - 1
  tmp = dble(kfg)/dble(az)
  kfp = CEILING(tmp) - 1

  procIJK(1) = i1p 
  procIJK(2) = ifp
  procIJK(3) = j1p
  procIJK(4) = jfp
  procIJK(5) = k1p
  procIJK(6) = kfp 

  !! Loop over all the processors and perform the operations
  CALL proc_loop(file,procIJK,stride,Gdata,FUNK)  

  !! Normalize by the slice area
  norm = dble( kfg - k1g + 1)
  Gdata = Gdata / norm

  !! Send back to main in the proper format
  prof2D = Gdata(:,:,1,:)
  
END SUBROUTINE PLANE_IJ

SUBROUTINE BLOCK(dataB,file)
  USE globals
  USE post_routines
  IMPLICIT NONE
  DOUBLE PRECISION, DIMENSION(:,:,:,:), INTENT(INOUT) :: dataB
  CHARACTER(LEN=flen) :: file

  !DOUBLE PRECISION, DIMENSION(:,:,:,:), ALLOCATABLE :: Gdata
  INTEGER, DIMENSION(3,2) :: stride
  INTEGER, DIMENSION(6) :: procIJK
  CHARACTER(LEN=flen) :: FUNK
  DOUBLE PRECISION :: norm,tmp
  INTEGER :: i1p,ifp,j1p,jfp,k1p,kfp
  INTEGER :: i1g,ifg,j1g,jfg,k1g,kfg
  INTEGER :: i

  !! Set the case for Summation rule
  FUNK = 'BLOCK'

  !! Set the index box size for the averaging kernal
  i1g = 1
  ifg = nx
  j1g = 1
  jfg = ny
  k1g = 1
  kfg = nz


  !! Set the processor bounds for the data
  tmp = dble(i1g)/dble(ax)
  i1p = CEILING(tmp) - 1
  tmp = dble(ifg)/dble(ax)
  ifp = CEILING(tmp) - 1

  tmp = dble(j1g)/dble(ay)
  j1p = CEILING(tmp) - 1
  tmp = dble(jfg)/dble(ay)
  jfp = CEILING(tmp) - 1

  tmp = dble(k1g)/dble(az)
  k1p = CEILING(tmp) - 1
  tmp = dble(kfg)/dble(az)
  kfp = CEILING(tmp) - 1

  procIJK(1) = i1p 
  procIJK(2) = ifp
  procIJK(3) = j1p
  procIJK(4) = jfp
  procIJK(5) = k1p
  procIJK(6) = kfp 

  !! Loop over all the processors and perform the operations
  CALL proc_loop(file,procIJK,stride,dataB,FUNK)  


  
END SUBROUTINE BLOCK


SUBROUTINE POINT(dataA,file,Ii,Ij,Ik)
  USE globals
  USE post_routines
  IMPLICIT NONE
  DOUBLE PRECISION, DIMENSION(DIM), INTENT(OUT) :: dataA
  CHARACTER(LEN=flen) :: file
  INTEGER, INTENT(IN) :: Ii,Ij,Ik

  DOUBLE PRECISION, DIMENSION(1,1,1,DIM) :: dataB
  INTEGER, DIMENSION(3,2) :: stride
  INTEGER, DIMENSION(6) :: procIJK
  CHARACTER(LEN=flen) :: FUNK
  DOUBLE PRECISION :: norm,tmp
  INTEGER :: i1p,ifp,j1p,jfp,k1p,kfp
  INTEGER :: i1g,ifg,j1g,jfg,k1g,kfg
  INTEGER :: i

  !! Set the case for Summation rule
  FUNK = 'POINT'

  !! Set the index box size for the averaging kernal
  i1g = Ii
  ifg = Ii
  j1g = Ij
  jfg = Ij
  k1g = Ik
  kfg = Ik
  
  stride(1,1) = Ii
  stride(1,2) = Ii
  stride(2,1) = Ij
  stride(2,2) = Ij
  stride(3,1) = Ik
  stride(3,2) = Ik



  !! Set the processor bounds for the data
  tmp = dble(i1g)/dble(ax)
  i1p = CEILING(tmp) - 1
  tmp = dble(ifg)/dble(ax)
  ifp = CEILING(tmp) - 1

  tmp = dble(j1g)/dble(ay)
  j1p = CEILING(tmp) - 1
  tmp = dble(jfg)/dble(ay)
  jfp = CEILING(tmp) - 1

  tmp = dble(k1g)/dble(az)
  k1p = CEILING(tmp) - 1
  tmp = dble(kfg)/dble(az)
  kfp = CEILING(tmp) - 1

  procIJK(1) = i1p 
  procIJK(2) = ifp
  procIJK(3) = j1p
  procIJK(4) = jfp
  procIJK(5) = k1p
  procIJK(6) = kfp 

  !! Loop over all the processors and perform the operations
  CALL proc_loop(file,procIJK,stride,dataB,FUNK)  

  dataA = dataB(1,1,1,:)

  
END SUBROUTINE POINT
