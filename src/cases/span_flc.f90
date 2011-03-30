!  Britton J. Olson
!  Stanford University
!  March. 4th, 2011

!  Description: Compute the span averaged variables
MODULE span_flc_data
  SAVE
  
  INTEGER :: ix1 = 1        ! X-index left bound on averaging volume
  INTEGER :: ixn = 100      ! X-index right bound on averaging volume
  INTEGER :: iy1 = 1        ! Y-index left bound on averaging volume
  INTEGER :: iyn = 100      ! Y-index right bound on averaging volume
  INTEGER :: iz1 = 1        ! Z-index left bound on averaging volume
  INTEGER :: izn = 100      ! Z-index right bound on averaging volume

END MODULE span_flc_data


PROGRAM span_flc
  USE operators
  USE globals
  USE span_flc_data
  IMPLICIT NONE
  INTEGER :: i,j,k
  INTEGER, DIMENSION(2,2) :: stride

  DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE :: pmean,pave,pflc
  CHARACTER(LEN=flen) :: ofile
  CHARACTER(LEN=flen) :: vfile
  CHARACTER(LEN=flen) :: inputFile

  !! Spatial correlation stuff
  DOUBLE PRECISION :: smean,count
  INTEGER :: ivar
  INTEGER :: rr,ii

  INTEGER, DIMENSION(:), ALLOCATABLE :: vars
  CHARACTER(LEN=flen) :: cvars
  
  INTEGER :: nviz,iviz,funit=34,tfreq
  NAMELIST /span_flc_vars/ iz1,izn



  !! General input file for t1,tf,nx,ny,nz,filename
  CALL read_input

  ! Read in namelist for specific case stuff  
  CALL GETARG(1,inputFile)
  OPEN(UNIT=funit,FILE=TRIM(inputFile),FORM='FORMATTED',STATUS='OLD')
  READ(UNIT=funit,NML=span_flc_vars)
  CLOSE(funit)


  tfreq = 1
  nviz=tf-t1+1
  nviz = nviz / tfreq

  !! Initialize the data arrays
  ALLOCATE(pmean(nx,ny,DIM))
  ALLOCATE(pave(nx,ny,DIM))
  ALLOCATE(pflc(nx,ny,DIM))
  
  pave = zero
  pflc = zero
  count = zero
  
  !! Main loop over the viz-files
  DO i=1,nviz
     iviz = i*tfreq + t1 - 1
     ! Get the next viz-directory and call kernal
     CALL viz_name(jobdir,iviz,vfile)
     CALL get_mean(vfile,pmean)
     pave = pave + pmean
     count = count + one
     
     CALL get_flc(vfile,pmean)
     pflc = pflc + pmean

  END DO

  ! Averaged in time and space
  pmean = pave / count

!!$  pflc = zero
!!$  count = zero
!!$  !! Main loop over the viz-files
!!$  DO i=1,nviz
!!$     iviz = i*tfreq + t1 - 1
!!$     pave = pmean
!!$     ! Get the next viz-directory and call kernal
!!$     CALL viz_name(jobdir,iviz,vfile)
!!$     CALL get_flc(vfile,pave)
!!$     pflc = pflc + pave
!!$     count = count + one
!!$  END DO
  
  ! Averaged in time and space
  pflc = pflc / count

  !! Put the proper grid info here.
  pflc(:,:,x_c:y_c) = pmean(:,:,x_c:y_c)
  pflc(:,:,z_c) = zero


  ofile = 'autocorr.tec'
  ALLOCATE(vars(6))
  vars = (/ u, v, w, rho, p, T /)
  cvars = ' "u","v","w","rho","p","T" '
  CALL write_tec(ofile,pmean,nx,ny,6,vars,cvars)

 END PROGRAM span_flc

SUBROUTINE get_plane(vfile,j1,jf,plane2d)
  USE operators
  USE globals
  USE span_flc_data
  IMPLICIT NONE
  INTEGER :: j1, jf, var
  CHARACTER(LEN=flen), INTENT(IN) :: vfile
  DOUBLE PRECISION, DIMENSION(nx,ny,DIM), INTENT(OUT) :: plane2d

  INTEGER, DIMENSION(2,2) :: stride
  
  stride(1,1) = j1
  stride(1,2) = jf

  CALL PLANE_IJ(stride,plane2d,vfile)
  

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
  !WRITE(ounit,*) ' VARIABLES = "X", "Y", "Z",',TRIM(cvars)
  WRITE(ounit,*) ' VARIABLES = "X", "Y", ',TRIM(cvars)
  WRITE(ounit,*) ' ZONE I=', n1, ', J=', n2, ', F=POINT'
  
  DO j=1,n2
     DO i=1,n1
        WRITE(ounit,*) data(i,j,x_c),data(i,j,y_c)!,data(i,j,z_c)
        DO iv=1,nV
           iiv = vars(iv)
           WRITE(ounit,*) data(i,j,iiv)
        END DO
     END DO
  END DO



END SUBROUTINE write_tec

SUBROUTINE get_mean(vfile,spc)
  USE operators
  USE globals
  USE span_flc_data
  IMPLICIT NONE
  CHARACTER(LEN=flen), INTENT(IN) :: vfile
  DOUBLE PRECISION, DIMENSION(nx,ny,DIM), INTENT(OUT) :: spc

  ! Just get the span averaged plane
  CALL get_plane(vfile,1,nz,spc)


END SUBROUTINE get_mean


SUBROUTINE get_flc(vfile,spc)
  USE operators
  USE globals
  USE span_flc_data
  IMPLICIT NONE
  CHARACTER(LEN=flen), INTENT(IN) :: vfile
  DOUBLE PRECISION, DIMENSION(nx,ny,DIM), INTENT(INOUT) :: spc

  DOUBLE PRECISION, DIMENSION(nx,ny,DIM) :: mean
  DOUBLE PRECISION, DIMENSION(nx,ny,DIM) :: flc
  DOUBLE PRECISION, DIMENSION(nx,ny,DIM) :: corr  
  INTEGER :: i

  mean = spc

  
  !! Loop over the k indices and accumulate the flc terms
  corr = zero
  DO i=1,nz
     CALL get_plane(vfile,i,i,spc)
     flc = spc - mean 
     corr = corr + flc**2
  END DO
  corr = corr / dble(nz)
  
  spc = corr

END SUBROUTINE get_flc
