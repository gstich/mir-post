SUBROUTINE SUBSUM3IK(stride,prof1D,file)
  USE globals
  USE post_routines
  IMPLICIT NONE
  DOUBLE PRECISION, DIMENSION(:,:) :: prof1D
  INTEGER, DIMENSION(2,2) :: stride
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
  k1g = stride(2,1)
  kfg = stride(2,2)

  !! Set the case for Summation rule
  FUNK = 'SUBSUM3IK'
  
  !! Allocate and initialize
  ALLOCATE(Gdata(1,ny,1,DIM))
  Gdata = 0.0D0

  !! Set the processor bounds for the data
  tmp = dble(i1g)/dble(ax)
  i1p = CEILING(tmp) - 1
  tmp = dble(ifg)/dble(ax)
  ifp = CEILING(tmp) - 1

  j1p = 0       !j1g/ay
  jfp = py-1    !jfg/ay

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

  !CALL print_map
  
  !! Loop over all the processors and perform the operations
  CALL proc_loop(file,procIJK,stride,Gdata,FUNK)

  !! Normalize by the slice area
  norm = dble( ifg - i1g + 1) * dble( kfg - k1g + 1)
  Gdata = Gdata / norm

  !! Send back to main in the proper format
  DO i=1,ny
     prof1D(i,:) = Gdata(1,i,1,:)
  END DO
     

END SUBROUTINE SUBSUM3IK


SUBROUTINE SUBSUM3JK(stride,prof1D,file)
  USE globals
  USE post_routines
  IMPLICIT NONE
  DOUBLE PRECISION, DIMENSION(:,:) :: prof1D
  INTEGER, DIMENSION(2,2), INTENT(IN) :: stride
  CHARACTER(LEN=flen) :: file

  DOUBLE PRECISION, DIMENSION(:,:,:,:), ALLOCATABLE :: Gdata
  INTEGER, DIMENSION(6) :: procIJK
  CHARACTER(LEN=flen) :: FUNK
  DOUBLE PRECISION :: norm,tmp
  INTEGER :: i1p,ifp,j1p,jfp,k1p,kfp
  INTEGER :: i1g,ifg,j1g,jfg,k1g,kfg
  INTEGER :: i


  !! Set the index box size for the averaging kernal
  j1g = stride(1,1)
  jfg = stride(1,2)
  k1g = stride(2,1)
  kfg = stride(2,2)

  !! Set the case for Summation rule
  FUNK = 'SUBSUM3JK'

  !! Allocate and initialize
  ALLOCATE(Gdata(1,nx,1,DIM))
  Gdata = 0.0D0  

  !! Set the processor bounds for the data
  i1p = 0       
  ifp = px-1    

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
  norm = dble( jfg - j1g + 1) * dble( kfg - k1g + 1)
  Gdata = Gdata / norm

  !! Send back to main in the proper format
  DO i=1,nx
     prof1D(i,:) = Gdata(i,1,1,:)
  END DO  

END SUBROUTINE SUBSUM3JK


SUBROUTINE SUBSUM3IJ(stride,prof1D,file)
  USE globals
  USE post_routines
  IMPLICIT NONE
  DOUBLE PRECISION, DIMENSION(:,:) :: prof1D
  INTEGER, DIMENSION(2,2), INTENT(IN) :: stride
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
  j1g = stride(2,1)
  jfg = stride(2,2)

  !! Set the case for Summation rule
  FUNK = 'SUBSUM3IJ'

  !! Allocate and initialize
  ALLOCATE(Gdata(1,1,nz,DIM))
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

  k1p = 0       
  kfp = pz-1    

  procIJK(1) = i1p 
  procIJK(2) = ifp
  procIJK(3) = j1p
  procIJK(4) = jfp
  procIJK(5) = k1p
  procIJK(6) = kfp 

  !! Loop over all the processors and perform the operations
  CALL proc_loop(file,procIJK,stride,Gdata,FUNK)  

  !! Normalize by the slice area
  norm = dble( ifg - i1g + 1) * dble( jfg - j1g + 1)
  Gdata = Gdata / norm

  !! Send back to main in the proper format
  DO i=1,nz
     prof1D(i,:) = Gdata(1,1,i,:)
  END DO  

END SUBROUTINE SUBSUM3IJ
