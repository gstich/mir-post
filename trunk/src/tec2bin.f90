! Convert the tecplot format to binary
PROGRAM tec2bin
  IMPLICIT NONE

  INTEGER, PARAMETER :: flen = 100
  INTEGER, PARAMETER :: funit = 21
  CHARACTER(LEN=flen) :: ftec,fbin,comments
  INTEGER :: nx,ny,nvar
  INTEGER :: i,j,v
  DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE :: data

  ftec = '/p/lscratchd/olson45/nozzle/nozzlecoarse3d/vis0020/planes/mean.tec'
  fbin = '/p/lscratchd/olson45/nozzle/nozzlecoarse3d/vis0020/planes/mean.bin'
  nx = 512
  ny = 128
  nvar = 8

  ALLOCATE(data(nx,ny,nvar))

  OPEN(UNIT=funit,FILE=ftec,FORM='FORMATTED',STATUS='UNKNOWN')
  READ(funit,*) comments
  READ(funit,*) comments
  
  DO j=1,ny
     DO i=1,nx
        READ(funit,*) data(i,j,1),data(i,j,2)
        DO v=3,nvar
           READ(funit,*) data(i,j,v)
        END DO
     END DO
  END DO
  CLOSE(funit)


  ! Dump the bin here
  OPEN(UNIT=funit,FILE=fbin,FORM='FORMATTED',STATUS='UNKNOWN')
  WRITE(funit,*) data
  CLOSE(funit)

  DEALLOCATE(data)
  

END PROGRAM tec2bin
