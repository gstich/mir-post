


SUBROUTINE read_input
  USE globals, ONLY: t1,tf,nx,ny,nz,ax,ay,az,px,py,pz
  USE globals, ONLY: flen,jobdir,invprocmap,procmap,proc
  IMPLICIT NONE
  INTEGER :: funit,ioUnit
  CHARACTER(LEN=flen) ::  inputFile           ! input file tmp
  CHARACTER(LEN=flen) ::  iodir
  CHARACTER(LEN=flen) ::  comments
  INTEGER :: i,j,k,p
  INTEGER, DIMENSION(:,:), ALLOCATABLE :: procmap_3d
  
  NAMELIST /INPUT/ t1,tf,nx,ny,nz,jobdir

  funit  = 11
  ioUnit = 17 

  ! Check to make sure an input file was given
  IF (iargc() .EQ. 0) STOP 'Usage: ./exec inputfile'
  CALL GETARG(1,inputFile)

  ! Read in namelist file and use to read in grid and setup
  ! arrays for the solver

  OPEN(UNIT=funit,FILE=TRIM(inputFile),FORM='FORMATTED',STATUS='OLD')
  READ(UNIT=funit,NML=INPUT)
  CLOSE(funit)


  !  Get Overall Domain size and figure out the extruded direction
  PRINT*,'2D Domain size:  nx=',nx,'ny=',ny,'nz=',nz

  IF(nx == 1 .OR. ny == 1 .OR. nz == 1) THEN
     PRINT*,'ERROR: This is 2D set, enter 3d set'
     STOP
  END IF

  
  !  Read in the procmap 
  WRITE(iodir,'(2A)') TRIM(jobdir),'/procmap'
  OPEN(UNIT=ioUnit,FILE=TRIM(iodir),FORM='FORMATTED',STATUS='UNKNOWN')
  READ(ioUnit,*) pz, py, px
  CLOSE(ioUnit)
  !READ(ioUnit,*) comments
  print*,px,py,pz
  proc = px*py*pz
  ax = nx/px;ay = ny/py;az = nz/pz
  ALLOCATE(procmap(proc,4))
  ALLOCATE(procmap_3d(proc,4))
  ALLOCATE(invprocmap(px,py,pz))
  
  !  Get the proc IJK layout
  p = 1
  DO k=1,pz
     DO j=1,py
        DO i=1,px
           procmap_3d(p,1)=p-1
           procmap_3d(p,2)=i-1
           procmap_3d(p,3)=j-1
           procmap_3d(p,4)=k-1
           invprocmap(i,j,k) = p-1
           p = p + 1
        END DO
     END DO
  END DO
  

  procmap = procmap_3d
  

END SUBROUTINE read_input



SUBROUTINE viz_name(file,nviz,vizfile)
  USE globals, ONLY: flen
  IMPLICIT NONE  
  CHARACTER(LEN=flen) :: vizfile,file
  INTEGER :: nviz
  WRITE(vizfile,'(2A,I4.4)') TRIM(file),'/vis',nviz
END SUBROUTINE viz_name



SUBROUTINE proc_loop(vizdir,procIJK,stride,Gdata,FUNK)
  USE globals, ONLY: flen, varDIM, DIM, invprocmap,ax,ay,az,zero,one
  USE block
  IMPLICIT NONE
  CHARACTER(LEN=flen), INTENT(IN) :: vizdir
  INTEGER, DIMENSION(6) :: procIJK
  INTEGER, DIMENSION(2,2) :: stride
  DOUBLE PRECISION, DIMENSION(:,:,:,:), INTENT(INOUT) :: Gdata
  CHARACTER(LEN=flen), INTENT(IN) :: FUNK


  REAL(KIND=4), DIMENSION(ax,ay,az,varDIM+3) :: iodata

  CHARACTER(LEN=flen) :: procdir
  INTEGER :: i,j,k,p,v
  INTEGER :: punit=27
  INTEGER :: ig,jg,kg
  INTEGER :: i1p,ifp,j1p,jfp,k1p,kfp

  INTEGER, DIMENSION(3) :: corner

  !ALLOCATE(iodata(ax,ay,ax,varDIM))
  i1p = procIJK(1)
  ifp = procIJK(2)
  j1p = procIJK(3)
  jfp = procIJK(4)
  k1p = procIJK(5)
  kfp = procIJK(6)


  ! Loop only over PROCS in pertainate range
  DO i=i1p,ifp
     DO j=j1p,jfp
        DO k=k1p,kfp
           

           
           ! Global indices on corner of proc grid
           ig = i*ax + 1
           jg = j*ay + 1
           kg = k*az + 1
           corner(1)=ig;corner(2)=jg;corner(3)=kg;
           
           !! i,j,k are (0,px0-1).... p is (0,nproc-1)
           p = invprocmap(i+1,j+1,k+1)


           WRITE(procdir,'(2A,I6.6)') TRIM(vizdir),'/p',p
           print*,'proc read',p,ig,jg,kg
           OPEN(UNIT=punit,FILE=TRIM(procdir),FORM='UNFORMATTED',STATUS='OLD')
           DO v=1,varDIM
              READ(punit) iodata(:,:,:,v)
           END DO
           CLOSE(punit)

           !! Get the grid data as well
           WRITE(procdir,'(2A,I6.6)') TRIM(vizdir),'/../grid/p',p
           OPEN(UNIT=punit,FILE=TRIM(procdir),FORM='UNFORMATTED',STATUS='OLD')
           DO v=varDIM+1,DIM
              READ(punit) iodata(:,:,:,v)
           END DO
           CLOSE(punit)
           
           !! Call the main kernal on the chunk of data
           CALL block_ops(iodata,corner,stride,Gdata,FUNK)

        END DO
     END DO
  END DO

END SUBROUTINE proc_loop


