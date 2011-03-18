SUBROUTINE SUBSUM3IK(stride,vars,prof1D,file)
  USE globals
  !USE sumdata
  IMPLICIT NONE
  DOUBLE PRECISION, DIMENSION(:,:) :: prof1D
  INTEGER, DIMENSION(:) :: vars
  INTEGER, DIMENSION(2,2) :: stride
  CHARACTER(LEN=flen) :: file

  DOUBLE PRECISION, DIMENSION(:,:,:,:), ALLOCATABLE :: Gdata
  INTEGER, DIMENSION(SIZE(vars)) :: voff
  INTEGER, DIMENSION(6) :: procIJK
  DOUBLE PRECISION :: tmp
  INTEGER :: i1p,ifp,j1p,jfp,k1p,kfp
  INTEGER :: i1g,ifg,j1g,jfg,k1g,kfg
  INTEGER :: i

  i1g = stride(1,1)
  ifg = stride(1,2)
  
  k1g = stride(2,1)
  kfg = stride(2,2)

  !fSS3IK = .TRUE.
  !ALLOCATE(out(ny,varDIM))
  !ALLOCATE(gout(ny,3))
  
  ALLOCATE(Gdata(1,ny,1,varDIM+3))
  Gdata = 0.0D0

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


  CALL print_map
   
  !norm = dble( ifg - i1g + 1) * dble( kfg - k1g + 1)


  CALL proc_loop(file,procIJK,stride,Gdata,'SUMSUM3IK')


  DO i=1,ny
     prof1D(i,1:3) = Gdata(1,i,1,varDIM+1:varDIM+3)
     prof1D(i,4:3+varDIM) = Gdata(1,i,1,1:varDIM)
  END DO
     

END SUBROUTINE SUBSUM3IK

!!$
!!$SUBROUTINE SUBSUM3JK(stride,vars,prof1D)
!!$  USE globals
!!$  USE sumdata
!!$  IMPLICIT NONE
!!$  DOUBLE PRECISION, DIMENSION(:,:) :: prof1D
!!$  INTEGER, DIMENSION(:) :: vars
!!$  INTEGER, DIMENSION(2,2), INTENT(IN) :: stride
!!$
!!$  INTEGER, DIMENSION(SIZE(vars)) :: voff
!!$  DOUBLE PRECISION :: tmp
!!$
!!$  j1g = stride(1,1)
!!$  jfg = stride(1,2)
!!$  
!!$  k1g = stride(2,1)
!!$  kfg = stride(2,2)
!!$
!!$  fSS3JK = .TRUE.
!!$  ALLOCATE(out(nx,varDIM))
!!$  ALLOCATE(gout(nx,3))
!!$  
!!$  i1p = 0       
!!$  ifp = px-1    
!!$
!!$  tmp = dble(j1g)/dble(ay)
!!$  j1p = CEILING(tmp) - 1
!!$  tmp = dble(jfg)/dble(ay)
!!$  jfp = CEILING(tmp) - 1
!!$
!!$  tmp = dble(k1g)/dble(az)
!!$  k1p = CEILING(tmp) - 1
!!$  tmp = dble(kfg)/dble(az)
!!$  kfp = CEILING(tmp) - 1
!!$
!!$  CALL print_map
!!$    
!!$  
!!$  norm = dble( jfg - j1g + 1) * dble( kfg - k1g + 1)
!!$
!!$  CALL viz_loop()
!!$
!!$
!!$  prof1D(:,1:3) = gout
!!$  voff = vars + 3
!!$  prof1D(:,voff) = out(:,vars)
!!$  
!!$
!!$  CALL clean_sum
!!$
!!$END SUBROUTINE SUBSUM3JK
!!$
!!$
!!$SUBROUTINE SUBSUM3IJ(stride,vars,prof1D)
!!$  USE globals
!!$  USE sumdata
!!$  IMPLICIT NONE
!!$  DOUBLE PRECISION, DIMENSION(:,:) :: prof1D
!!$  INTEGER, DIMENSION(2,2), INTENT(IN) :: stride
!!$  INTEGER, DIMENSION(:) :: vars
!!$
!!$  INTEGER, DIMENSION(SIZE(vars)) :: voff
!!$  DOUBLE PRECISION :: tmp
!!$
!!$  i1g = stride(1,1)
!!$  ifg = stride(1,2)
!!$  
!!$  j1g = stride(2,1)
!!$  jfg = stride(2,2)
!!$
!!$  fSS3IJ = .TRUE.
!!$  ALLOCATE(out(nz,varDIM))
!!$  ALLOCATE(gout(nz,3))
!!$  
!!$  tmp = dble(i1g)/dble(ax)
!!$  i1p = CEILING(tmp) - 1
!!$  tmp = dble(ifg)/dble(ax)
!!$  ifp = CEILING(tmp) - 1
!!$
!!$  tmp = dble(j1g)/dble(az)
!!$  j1p = CEILING(tmp) - 1
!!$  tmp = dble(jfg)/dble(az)
!!$  jfp = CEILING(tmp) - 1
!!$
!!$  k1p = 0       !j1g/ay
!!$  kfp = pz-1    !jfg/ay
!!$
!!$  CALL print_map
!!$
!!$  norm = dble( ifg - i1g + 1) * dble( jfg - j1g + 1)
!!$
!!$  CALL viz_loop()
!!$
!!$
!!$  prof1D(:,1:3) = gout
!!$  voff = vars + 3
!!$  prof1D(:,voff) = out(:,vars)
!!$  
!!$
!!$  CALL clean_sum
!!$
!!$END SUBROUTINE SUBSUM3IJ
