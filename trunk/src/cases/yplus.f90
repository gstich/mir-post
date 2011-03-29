!  Britton J. Olson
!  Stanford University
!  March. 4th, 2011

!  Description: Compute the mean velocity profiles for a turbulent boundary layer
!       and make the appropriate scaling in y+ and u_tau units
MODULE yplus_data
  SAVE
  
  INTEGER :: ix1 = 1        ! X-index left bound on averaging volume
  INTEGER :: ixn = 100      ! X-index right bound on averaging volume
  INTEGER :: iy1 = 1        ! Y-index left bound on averaging volume
  INTEGER :: iyn = 100      ! Y-index right bound on averaging volume
  INTEGER :: iz1 = 1        ! Z-index left bound on averaging volume
  INTEGER :: izn = 100      ! Z-index right bound on averaging volume

END MODULE yplus_data


PROGRAM yplus
  USE operators
  USE globals
  USE yplus_data
  IMPLICIT NONE
  INTEGER :: i
  INTEGER, DIMENSION(2,2) :: stride

  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: output
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: tave
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: tflc
  DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE :: tout
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: Uvd
  CHARACTER(LEN=flen) :: ofile
  CHARACTER(LEN=flen) :: vfile
  CHARACTER(LEN=flen) :: inputFile

  DOUBLE PRECISION :: dudy,tauw,mu_0,mu_w,rho_w,rho_0,del,utau,off,dUp
  INTEGER :: nviz,iviz,funit=34
  NAMELIST /yplus_vars/ ix1,ixn,iz1,izn



  !! General input file for t1,tf,nx,ny,nz,filename
  CALL read_input

  ! Read in namelist for specific case stuff  
  CALL GETARG(1,inputFile)
  OPEN(UNIT=funit,FILE=TRIM(inputFile),FORM='FORMATTED',STATUS='OLD')
  READ(UNIT=funit,NML=yplus_vars)
  CLOSE(funit)


  !! Set some variables from the input files
  stride(1,1) = ix1
  stride(1,2) = ixn
  stride(2,1) = iz1
  stride(2,2) = izn
  nviz=tf-t1+1

  !! Initialize the data arrays
  ALLOCATE(output( ny , DIM) )
  ALLOCATE(tout( ny , DIM, nviz) )
  ALLOCATE(tave( ny , DIM) )
  ALLOCATE(tflc( ny , DIM) )
  tflc = 0.0D0
  tave = 0.0D0

  !! Main loop over the viz-files
  DO i=1,nviz
     iviz = i + t1 - 1
     output = 0.0D0
     
     ! Get the next viz-directory and call kernal
     CALL viz_name(jobdir,iviz,vfile)
     CALL SUBSUM3IK(stride,output,vfile)

     ! Average data over time and save for fluctuating computation
     tout(:,:,i) = output
     tave = tave + output
  END DO
  tave = tave/DBLE(nviz)

  ! Get the temporal fluctuating part here
  DO i=1,nviz
     tflc = tflc + ( tave - tout(:,:,i))**2
  END DO
  tflc = tflc/dble(nviz)


  !! Get some y_plus units
  rho_w = output(1,rho)
  rho_0 = output(ny/2,rho)
  mu_0 = output(ny/2,mu)
  mu_w = output(1,mu)

  dudy = ( output(2,u) - output(1,u) )/ ( output(2,y_c) - output(1,y_c) )
  tauw = mu_w * dudy
  utau = sqrt( tauw / rho_w )
  print*,'U_tau',utau

  !! Scale y
  del = mu_w / ( rho_w * utau )
  output(:,y_c) = output(:,y_c) / del
  off = output(1,y_c)
  output(:,y_c) = output(:,y_c) - off
  print*,'Del',del


  !! Scale velocity
  output(:,u) = output(:,u) / utau


  !! Van Driest transformation
  ALLOCATE(Uvd(ny))
  Uvd(1) = 0.0D0
  DO i=2,ny/2
     dUp = output(i,u)-output(i-1,u)
     Uvd(i) = Uvd(i-1) + sqrt(output(i,rho)/rho_w) * dUp
  END DO

  Uvd(ny) = 0.0D0
  DO i=ny-1,ny/2+1,-1
     dUp = output(i,u)-output(i+1,u)
     Uvd(i) = Uvd(i+1) + sqrt(output(i,rho)/rho_w) * dUp
  END DO



  OPEN(UNIT=33,FILE='Uwall.dat',FORM='FORMATTED',STATUS='UNKNOWN')
  WRITE(33,*) '# y+, u+, rho, Uvd, <uu>'
  DO i=1,ny/2
     WRITE(33,'(6ES12.4)') output(i,y_c),output(i,u),output(i,rho),Uvd(i),tflc(i,u)
  END DO
  DO i=ny/2+1,ny
     WRITE(33,'(6ES12.4)') (output(ny,y_c)-output(i,y_c)),output(i,u),output(i,rho),Uvd(i),tflc(i,u)
  END DO

  CLOSE(33)





END PROGRAM yplus

