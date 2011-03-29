!  Britton J. Olson
!  Stanford University
!  March. 4th, 2011

!  Description: Compute the mean velocity profiles for a turbulent boundary layer
!       and make the appropriate scaling in y+ and u_tau units
MODULE length_data
  SAVE
  
  INTEGER :: ix1 = 1        ! X-index left bound on averaging volume
  INTEGER :: ixn = 100      ! X-index right bound on averaging volume
  INTEGER :: iy1 = 1        ! Y-index left bound on averaging volume
  INTEGER :: iyn = 100      ! Y-index right bound on averaging volume
  INTEGER :: iz1 = 1        ! Z-index left bound on averaging volume
  INTEGER :: izn = 100      ! Z-index right bound on averaging volume

END MODULE length_data


PROGRAM length
  USE operators
  USE globals
  USE length_data
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

  INTEGER :: nviz,iviz,funit=34
  NAMELIST /length_vars/ iy1,iyn,iz1,izn



  !! General input file for t1,tf,nx,ny,nz,filename
  CALL read_input

  ! Read in namelist for specific case stuff  
  CALL GETARG(1,inputFile)
  OPEN(UNIT=funit,FILE=TRIM(inputFile),FORM='FORMATTED',STATUS='OLD')
  READ(UNIT=funit,NML=length_vars)
  CLOSE(funit)


  !! Set some variables from the input files
  stride(1,1) = iy1
  stride(1,2) = iyn
  stride(2,1) = iz1
  stride(2,2) = izn
  nviz=tf-t1+1

  !! Initialize the data arrays
  ALLOCATE(output( nx , DIM) )
  ALLOCATE(tout( nx , DIM, nviz) )
  ALLOCATE(tave( nx , DIM) )
  ALLOCATE(tflc( nx , DIM) )
  tflc = 0.0D0
  tave = 0.0D0

  !! Main loop over the viz-files
  DO i=1,nviz
     iviz = i + t1 - 1
     output = 0.0D0
     
     ! Get the next viz-directory and call kernal
     CALL viz_name(jobdir,iviz,vfile)
     CALL SUBSUM3JK(stride,output,vfile)

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



  OPEN(UNIT=33,FILE='xprofile.dat',FORM='FORMATTED',STATUS='UNKNOWN')
  WRITE(33,*) '# X, U, rho, P, T'
  DO i=1,nx
     WRITE(33,'(5ES12.4)') output(i,x_c),output(i,u),output(i,rho),output(i,p),output(i,T)
  END DO

  CLOSE(33)





END PROGRAM length

