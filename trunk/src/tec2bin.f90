! Convert the tecplot format to binary
PROGRAM tec2bin
  IMPLICIT NONE

  INTEGER, PARAMETER :: flen = 100
  INTEGER, PARAMETER :: funit = 21
  CHARACTER(LEN=flen) :: ftec,comments
  INTEGER :: nx,ny
  

  ftec = 

  OPEN(UNIT=funit,FILE=ftec,FORM='FORMATTED',STATUS='UNKNOWN')
  READ(funit,*) comments
  READ(funit,*) comments,comments,nx,comments,comments,ny

  CLOSE(funit)
