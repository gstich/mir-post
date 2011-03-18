!  Britton J. Olson
!  Stanford University
!  March. 4th, 2011

!  Description: Compute the mean velocity profiles for a turbulent boundary layer
!       and make the appropriate scaling in y+ and u_tau units
MODULE span_data
  SAVE
  
  INTEGER :: ix1 = 1        ! X-index left bound on averaging volume
  INTEGER :: ixn = 100      ! X-index right bound on averaging volume
  INTEGER :: iy1 = 1        ! Y-index left bound on averaging volume
  INTEGER :: iyn = 100      ! Y-index right bound on averaging volume
  INTEGER :: iz1 = 1        ! Z-index left bound on averaging volume
  INTEGER :: izn = 100      ! Z-index right bound on averaging volume

END MODULE span_data


PROGRAM span
  USE operators
  USE globals
  USE span_data
  IMPLICIT NONE
  INTEGER :: i
  INTEGER, DIMENSION(2,2) :: stride

  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: output
  CHARACTER(LEN=flen) :: ofile
  CHARACTER(LEN=flen) :: vfile
  CHARACTER(LEN=flen) :: inputFile

  !! Spatial correlation stuff
  DOUBLE PRECISION :: smean
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: sflc
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: corr
  INTEGER :: ivar
  INTEGER :: rr,ii

  INTEGER :: nviz,iviz,funit=34
  NAMELIST /span_vars/ ix1,ixn,iy1,iyn



  !! General input file for t1,tf,nx,ny,nz,filename
  CALL read_input

  ! Read in namelist for specific case stuff  
  CALL GETARG(1,inputFile)
  OPEN(UNIT=funit,FILE=TRIM(inputFile),FORM='FORMATTED',STATUS='OLD')
  READ(UNIT=funit,NML=span_vars)
  CLOSE(funit)


  !! Set some variables from the input files
  stride(1,1) = ix1
  stride(1,2) = ixn
  stride(2,1) = iy1
  stride(2,2) = iyn
  nviz=tf-t1+1

  !! Initialize the data arrays
  ALLOCATE(output( nz , DIM) )
  
  !! Main loop over the viz-files
  DO i=1,nviz
     iviz = i + t1 - 1
     output = 0.0D0
     
     ! Get the next viz-directory and call kernal
     CALL viz_name(jobdir,iviz,vfile)
     CALL SUBSUM3IJ(stride,output,vfile)

  END DO


  !! Use the profile, to make spatial correlation
  ivar = p
  ALLOCATE(sflc(nz*2))
  ALLOCATE(corr(nz))
  corr = 0.0D0
  smean = SUM( output(:,ivar) ) / dble(nz)
  sflc(1:nz) = output(:,ivar) - smean
  sflc(nz+1:2*nz) = sflc(1:nz)
  

  DO ii=1,nz
     DO rr=1,nz
        corr(rr) = corr(rr) + sflc(ii) * sflc( ii + rr - 1)
     END DO
  END DO
  corr = corr / corr(1)

  OPEN(UNIT=33,FILE='span.dat',FORM='FORMATTED',STATUS='UNKNOWN')
  WRITE(33,*) '# Z, U, rho, P, T'
  DO i=1,nz
     WRITE(33,'(5ES12.4)') output(i,z_c),output(i,u),output(i,rho),output(i,p),output(i,T)
  END DO

  CLOSE(33)

  OPEN(UNIT=33,FILE='corr.dat',FORM='FORMATTED',STATUS='UNKNOWN')
  WRITE(33,*) '# Z, uu'
  DO i=1,nz/2
     WRITE(33,'(5ES12.4)') output(i,z_c),corr(i)
  END DO

  CLOSE(33)




END PROGRAM span

