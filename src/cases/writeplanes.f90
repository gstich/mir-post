!  Britton J. Olson
!  Stanford University
!  March. 4th, 2011

!  Description: Compute the span averaged variables
MODULE writeplanes_data
  SAVE
  
  INTEGER :: ix1 = 1        ! X-index left bound on averaging volume
  INTEGER :: ixn = 100      ! X-index right bound on averaging volume
  INTEGER :: iy1 = 1        ! Y-index left bound on averaging volume
  INTEGER :: iyn = 100      ! Y-index right bound on averaging volume
  INTEGER :: iz1 = 1        ! Z-index left bound on averaging volume
  INTEGER :: izn = 100      ! Z-index right bound on averaging volume

END MODULE writeplanes_data


PROGRAM writeplanes
  USE operators
  USE globals
  USE writeplanes_data
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
  NAMELIST /writeplanes_vars/ iz1,izn

  CHARACTER(LEN=flen) :: DONEfile, command
  LOGICAL :: compute


  !! General input file for t1,tf,nx,ny,nz,filename
  CALL read_input

  ! Read in namelist for specific case stuff  
  CALL GETARG(1,inputFile)
  OPEN(UNIT=funit,FILE=TRIM(inputFile),FORM='FORMATTED',STATUS='OLD')
  READ(UNIT=funit,NML=writeplanes_vars)
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
     iviz = (i-1)*tfreq + t1

     ! Get the next viz-directory and call kernal
     CALL viz_name(jobdir,iviz,vfile)

     ! Write out the planes in the vizdir
     WRITE(DONEfile,'(2A)') TRIM(vfile),'/planes/DONE'
     INQUIRE(FILE=DONEfile, EXIST=compute)
     
     IF(.NOT. compute) THEN
        
        WRITE(command,'(3A)') 'mkdir ',TRIM(vfile),'/planes'
        CALL SYSTEM(command)
        

        PRINT*,'Computing planes for ',vfile
        CALL get_mean(vfile,pmean)
        pave = pmean
        
        CALL get_flc(vfile,pmean)
        pflc = pmean
        
        ! Fix the grid variables
        pflc(:,:,x_c:y_c) = pave(:,:,x_c:y_c)
        pflc(:,:,z_c) = zero
        
        ! Write files and flags for each viz_dump
        IF(.NOT. ALLOCATED(vars)) ALLOCATE(vars(6))
        vars = (/ u, v, w, rho, p, T /)
        cvars = ' "u","v","w","rho","p","T" '

        WRITE(ofile,'(2A)') TRIM(vfile),'/planes/mean.tec'
        CALL write_tec(ofile,pave,nx,ny,6,vars,cvars)
        
        WRITE(ofile,'(2A)') TRIM(vfile),'/planes/flc.tec'
        CALL write_tec(ofile,pflc,nx,ny,6,vars,cvars)

        WRITE(command,'(3A)') 'touch ',TRIM(vfile),'/planes/DONE'
        CALL SYSTEM(command)

     ELSE
        PRINT*,'Planes exist, skipping ',vfile
     END IF


  END DO




 END PROGRAM writeplanes

SUBROUTINE get_plane(vfile,j1,jf,plane2d)
  USE operators
  USE globals
  USE writeplanes_data
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
  USE writeplanes_data
  IMPLICIT NONE
  CHARACTER(LEN=flen), INTENT(IN) :: vfile
  DOUBLE PRECISION, DIMENSION(nx,ny,DIM), INTENT(OUT) :: spc

  ! Just get the span averaged plane
  CALL get_plane(vfile,1,nz,spc)


END SUBROUTINE get_mean


SUBROUTINE get_flc(vfile,spc)
  USE operators
  USE globals
  USE writeplanes_data
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
