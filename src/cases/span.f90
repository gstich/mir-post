!  Britton J. Olson
!  Stanford University
!  March. 4th, 2011

!  Description: Compute the span averaged variables
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
  INTEGER :: i,j,k
  INTEGER, DIMENSION(2,2) :: stride

  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: output,plane2d
  CHARACTER(LEN=flen) :: ofile
  CHARACTER(LEN=flen) :: vfile
  CHARACTER(LEN=flen) :: inputFile

  !! Spatial correlation stuff
  DOUBLE PRECISION :: smean,count
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: corrM
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: corr
  INTEGER :: ivar
  INTEGER :: rr,ii
  

  INTEGER :: nviz,iviz,funit=34,tfreq
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
  
  tfreq = 1
  nviz=tf-t1+1
  nviz = nviz / tfreq

  !! Initialize the data arrays
  ALLOCATE(output( nz , DIM) )
  
  ALLOCATE(corr(nz/2))
  ALLOCATE(corrM(nz/2))
  
  
  corrM = zero
  count = zero
  !! Main loop over the viz-files
  DO i=1,nviz
     
     iviz = i*tfreq + t1 - 1
     
     
     !Loop over plane bounds
     DO j=1,(ixn-ix1+1)
        DO k=1,(iyn-iy1+1)
           
           ! Get the next viz-directory and call kernal
           CALL viz_name(jobdir,iviz,vfile)

           CALL get_coor(vfile,ix1+j-1,iy1+k-1,corr,p)
     
           corrM = corrM + corr
           count = count + one
        END DO
     END DO
  END DO
  
  ! Averaged in time and space
  corr = corrM / count
        
  OPEN(UNIT=33,FILE='corr.dat',FORM='FORMATTED',STATUS='UNKNOWN')
  WRITE(33,*) '# Z, uu'
  DO i=1,nz/2
     WRITE(33,'(1ES12.4)') corr(i)
  END DO
  CLOSE(33)


  !ALLOCATE(plane2d(nx,nz))
  CALL get_plane(vfile,1,1,4)



END PROGRAM span

SUBROUTINE get_plane(vfile,j1,jf,var)
  USE operators
  USE globals
  USE span_data
  IMPLICIT NONE
  INTEGER :: j1, jf, var
  CHARACTER(LEN=flen), INTENT(IN) :: vfile
  !DOUBLE PRECISION, DIMENSION(nx,nz), INTENT(OUT) :: plane2d

  INTEGER :: i, j,iv, ounit
  INTEGER, DIMENSION(2,2) :: stride
  DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE :: output
  CHARACTER(LEN=flen) :: ofile
  INTEGER, DIMENSION(:), ALLOCATABLE :: vars
  CHARACTER(LEN=flen) :: cvars
  
!!$  stride(1,1) = j1
!!$  stride(1,2) = jf
!!$
!!$  ALLOCATE(output(nx,nz,DIM))
!!$  CALL PLANE_IK(stride,output,vfile)
  
  stride(1,1) = j1
  stride(1,2) = jf

  ALLOCATE(output(nx,ny,DIM))
  CALL PLANE_IJ(stride,output,vfile)
  
  
  ofile = 'plane.tec'
  ALLOCATE(vars(6))
  vars = (/ u, v, w, rho, p, T /)
  cvars = ' "u","v","w","rho","p","T" '
  CALL write_tec(ofile,output,nx,ny,6,vars,cvars)

END SUBROUTINE get_plane

SUBROUTINE write_tec(ofile,data,n1,n2,nV,vars,cvars)
  USE globals
  IMPLICIT NONE
  
  CHARACTER(LEN=flen), INTENT(IN) :: ofile
  INTEGER, INTENT(IN) :: n1,n2,nV  
  DOUBLE PRECISION, DIMENSION(n1,n2,DIM) :: data
  INTEGER, DIMENSION(nV) :: vars
  CHARACTER(LEN=flen) :: cvars
  
  INTEGER :: i,j,iv,iiv
  INTEGER :: ounit = 25
  
  OPEN(UNIT=ounit,FILE=ofile,FORM='FORMATTED',STATUS='UNKNOWN')
  WRITE(ounit,*) ' VARIABLES = "X", "Y", "Z",',TRIM(cvars)
  WRITE(ounit,*) ' ZONE I=', n1, ', J=', n2, ', F=POINT'
  
  DO j=1,n2
     DO i=1,n1
        WRITE(ounit,*) data(i,j,x_c),data(i,j,y_c),data(i,j,z_c)
        DO iv=1,nV
           iiv = vars(iv)
           WRITE(ounit,*) data(i,j,iiv)
        END DO
     END DO
  END DO



END SUBROUTINE write_tec



SUBROUTINE get_coor(vfile,xslc,yslc,spc,ivar)
  USE operators
  USE globals
  USE span_data
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: xslc,yslc
  CHARACTER(LEN=flen), INTENT(IN) :: vfile
  DOUBLE PRECISION, DIMENSION(nz/2), INTENT(OUT) :: spc
  INTEGER,INTENT(IN) :: ivar

  INTEGER :: i
  INTEGER, DIMENSION(2,2) :: stride
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: output
  CHARACTER(LEN=flen) :: ofile

  CHARACTER(LEN=flen) :: inputFile

  !! Spatial correlation stuff
  DOUBLE PRECISION :: smean
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: sflc
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: corr

  INTEGER :: rr,ii

  INTEGER :: nviz,iviz,funit=34


  !! Set some variables from the input files
  stride(1,1) = xslc
  stride(1,2) = xslc
  stride(2,1) = yslc
  stride(2,2) = yslc

  !! Initialize the data arrays
  ALLOCATE(output( nz , DIM) )
  
  CALL SUBSUM3IJ(stride,output,vfile)

  
  !! Use the profile, to make spatial correlation
  !ivar = p
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

  !OPEN(UNIT=33,FILE='corr.dat',FORM='FORMATTED',STATUS='UNKNOWN')
  !WRITE(33,*) '# Z, uu'
  !DO i=1,nz/2
  !   WRITE(33,'(5ES12.4)') output(i,z_c),corr(i)
  !END DO

  !CLOSE(33)
  spc = corr(1:nz/2)




END SUBROUTINE get_coor
