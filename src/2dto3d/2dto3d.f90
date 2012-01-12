!  Britton J. Olson
!  Stanford University
!  Jan. 5th, 2010

!  Program will generate 3-d restart file for miranda by
!  extruding in the 3rd direction and blending two 2-D restart 
!  files together for use in Miranda code


PROGRAM twod
IMPLICIT NONE
INTEGER, PARAMETER :: flen=90          ! file length
DOUBLE PRECISION, PARAMETER :: pi=3.14159265358979323846d0,one=1.0d0,two=2.0d0
INTEGER, DIMENSION(:), ALLOCATABLE :: ix,iy,iz,ix1,iy1,iz1
INTEGER, DIMENSION(:,:), ALLOCATABLE :: I_2DA,I_2DB
INTEGER, DIMENSION(:,:), ALLOCATABLE :: procmap_2d, procmap_3d
DOUBLE PRECISION :: blend,mode
DOUBLE PRECISION, DIMENSION(:,:,:,:), ALLOCATABLE :: iodata_2da,iodata_2db,iodata_3d,iodata_2dtmp

INTEGER :: i,j,k,p,ext,varDIM,funit


INTEGER :: nx, ny, nz
INTEGER :: nx3,ny3,nz3
INTEGER :: ax_2d, ay_2d, az_2d
INTEGER :: ax_3d, ay_3d, az_3d
INTEGER :: px_2d, py_2d, pz_2d
INTEGER :: px_3d, py_3d, pz_3d
INTEGER :: proc2d,proc3d,in3d

CHARACTER(LEN=flen)  comments
CHARACTER(LEN=flen)  iodir2d             ! 2d restart and viz directories
CHARACTER(LEN=flen)  jobdir2d            ! 2d job path
INTEGER, DIMENSION(2) :: resnum          ! Restart numbers to use for blending
CHARACTER(LEN=flen)  iodir3d             ! 3d restart and viz directories
CHARACTER(LEN=flen)  jobdir3d            ! 3d job path
CHARACTER(LEN=flen)  inputFile           ! input file tmp
CHARACTER(LEN=7)  mkdir                  ! mkdir
INTEGER, PARAMETER :: io2Unit=12,io3Unit=13
NAMELIST /INPUT/ nx,ny,nz,nx3,ny3,nz3,jobdir2d,resnum,jobdir3d,px_3d,py_3d,pz_3d


varDIM=8
mode = 8.0d0     ! Number of modes fluctuating between restart dumps

! Check to make sure an input file was given
IF (iargc() .EQ. 0) STOP 'Usage: ./exec inputfile'
CALL GETARG(1,inputFile)

! Read in namelist file and use to read in grid and setup
! arrays for the solver
funit=11
OPEN(UNIT=funit,FILE=TRIM(inputFile),FORM='FORMATTED',STATUS='OLD')
  READ(UNIT=funit,NML=INPUT)
CLOSE(funit)


!  Get Overall Domain size and figure out the extruded direction
PRINT*,'2D Domain size:  nx=',nx,'ny=',ny,'nz=',nz

IF(nx .ne. 1 .and. ny .ne. 1 .and. nz .ne. 1) THEN
   PRINT*,'ERROR: This is 3D set, enter 2d set'
   STOP
ELSE
   IF(nx==1) THEN
      PRINT*,'Enter 3d nx=';!READ*,nx3;
      ny3=ny;nz3=nz;ext=1
   ELSE IF(ny==1) THEN
      PRINT*,'Enter 3d ny=';!READ*,ny3;
      nx3=nx;nz3=nz;ext=2
   ELSE IF(nz==1) THEN
      PRINT*,'Enter 3d nz=';!READ*,nz3;
      ny3=ny;nx3=nx;ext=3
   END IF
END IF

PRINT*,'3D Domain size:  nx=',nx3,'ny=',ny3,'nz=',nz3

!  Get the 2d-file location & restart files
PRINT*,"2d Path: ",jobdir2d


!  Read in the procmap data for 2D case
WRITE(iodir2d,'(2A)') TRIM(jobdir2d),'/procmap'
OPEN(UNIT=io2Unit,FILE=TRIM(iodir2d),FORM='FORMATTED',STATUS='UNKNOWN')
READ(io2Unit,*) pz_2d, py_2d, px_2d
READ(io2Unit,*) comments
print*,px_2d,py_2d,pz_2d
proc2d = px_2d*py_2d*pz_2d
ax_2d = nx/px_2d;ay_2d = ny/py_2d;az_2d = nz/pz_2d
ALLOCATE(procmap_2d(proc2d,4))
DO i=1,proc2d
!**!! Error in procmap file.. ordering is reversed... go with it
  READ(io2Unit,*) procmap_2d(i,1),procmap_2d(i,4),procmap_2d(i,3),procmap_2d(i,2)
END DO
CLOSE(io2Unit)



!  Set the procmap for the 3d grid
!  Get the 3d-file location
!PRINT*,"Enter the path of the 3d data directory: '/full/path/here'"
!READ*, jobdir3d
PRINT*,"3d Path: ",jobdir3d

proc3d = px_3d*py_3d*pz_3d
ax_3d = nx3/px_3d;ay_3d = ny3/py_3d;az_3d = nz3/pz_3d
ALLOCATE(procmap_3d(proc3d,4))

p = 1
DO k=1,pz_3d
   DO j=1,py_3d
      DO i=1,px_3d
         procmap_3d(p,1)=p-1
         procmap_3d(p,2)=i-1
         procmap_3d(p,3)=j-1
         procmap_3d(p,4)=k-1
         p = p + 1
      END DO
   END DO
END DO

ALLOCATE(iodata_2dtmp(ax_2d,ay_2d,az_2d,varDIM))
ALLOCATE(iodata_2da(nx,ny,nz,varDIM))
ALLOCATE(iodata_2db(nx,ny,nz,varDIM))
ALLOCATE(ix(ax_2d))
ALLOCATE(iy(ay_2d))
ALLOCATE(iz(az_2d))
ALLOCATE(ix1(ax_2d))
ALLOCATE(iy1(ay_2d))
ALLOCATE(iz1(az_2d))

DO i=1,ax_2d;ix1(i)=i;END DO
DO j=1,ay_2d;iy1(j)=j;END DO
DO k=1,az_2d;iz1(k)=k;END DO
!  Read in both 2d data sets
!    set-a
DO p=1,proc2d
  WRITE(iodir2d,'(2A,I4.4)') TRIM(jobdir2d),'/res',resnum(1)
  WRITE(iodir2d,'(2A,I6.6)') TRIM(iodir2d),'/p',p-1
  OPEN(UNIT=io2Unit,FILE=TRIM(iodir2d),FORM='UNFORMATTED',STATUS='OLD',action='read')
    READ(io2Unit) iodata_2dtmp
  CLOSE(io2Unit) 
  ix = ix1 + procmap_2d(p,2)*ax_2d
  iy = iy1 + procmap_2d(p,3)*ay_2d
  iz = iz1 + procmap_2d(p,4)*az_2d
  iodata_2da(ix,iy,iz,:) = iodata_2dtmp  
END DO

!    set-b
DO p=1,proc2d
  WRITE(iodir2d,'(2A,I4.4)') TRIM(jobdir2d),'/res',resnum(2)
  WRITE(iodir2d,'(2A,I6.6)') TRIM(iodir2d),'/p',p-1
  OPEN(UNIT=io2Unit,FILE=TRIM(iodir2d),FORM='UNFORMATTED',STATUS='OLD',action='read')
    READ(io2Unit) iodata_2dtmp
  CLOSE(io2Unit) 
  ix = ix1 + procmap_2d(p,2)*ax_2d
  iy = iy1 + procmap_2d(p,3)*ay_2d
  iz = iz1 + procmap_2d(p,4)*az_2d
  iodata_2db(ix,iy,iz,:) = iodata_2dtmp  
END DO

!  Blend and send to the individual 3D procs
DEALLOCATE(ix,iy,iz,ix1,iy1,iz1)
ALLOCATE(ix(ax_3d))
ALLOCATE(iy(ay_3d))
ALLOCATE(iz(az_3d))
ALLOCATE(ix1(ax_3d))
ALLOCATE(iy1(ay_3d))
ALLOCATE(iz1(az_3d))
ALLOCATE(iodata_3d(ax_3d,ay_3d,az_3d,varDIM))

DO i=1,ax_3d;ix1(i)=i;END DO
DO j=1,ay_3d;iy1(j)=j;END DO
DO k=1,az_3d;iz1(k)=k;END DO



!  For each proc, blend solutions together and write the restart file
mkdir = "mkdir '"
WRITE(iodir3d,'(2A,I4.4)') mkdir, TRIM(jobdir3d)
WRITE(iodir3d,'(2A,I4.4)') TRIM(iodir3d),"'"
CALL system(iodir3d)

WRITE(iodir3d,'(2A,I4.4)') mkdir, TRIM(jobdir3d)
WRITE(iodir3d,'(2A,I4.4)') TRIM(iodir3d),"/res_init'"
CALL system(iodir3d)

DO p=1,proc3d
   WRITE(iodir3d,'(2A,I4.4)') TRIM(jobdir3d),'/res_init'
   WRITE(iodir3d,'(2A,I6.6)') TRIM(iodir3d),'/p',p-1
   ix = ix1 + procmap_3d(p,2)*ax_3d
   iy = iy1 + procmap_3d(p,3)*ay_3d
   iz = iz1 + procmap_3d(p,4)*az_3d

   !  Get blending function 
   DO k=1,az_3d
      blend = (one + cos( mode*two*pi*dble(iz(k)-1)/dble(nz3-1) ))/two
      iodata_3d(:,:,k,:) = iodata_2da(ix,iy,1,:)*blend + iodata_2db(ix,iy,1,:)*(one-blend)
   END DO

   OPEN(UNIT=io3Unit,FILE=TRIM(iodir3d),FORM='UNFORMATTED',status='unknown',action='write')
   WRITE(io3Unit) iodata_3d
   CLOSE(io3Unit)
   PRINT*,p

END DO
PRINT*, 'DONE!  Files written to:', jobdir3d


END PROGRAM twod
