MODULE globals
  IMPLICIT NONE
  INTEGER, PARAMETER :: flen=90          ! Default string length

  !! Variable aliases and definitions
  INTEGER, PARAMETER :: varDIM=12        ! Number of variables in viz dumps
  INTEGER, PARAMETER :: u=1              ! U-velocity index
  INTEGER, PARAMETER :: v=2              ! V-velocity index
  INTEGER, PARAMETER :: w=3              ! W-velocity index
  INTEGER, PARAMETER :: rho=4            ! density index
  INTEGER, PARAMETER :: e=5              ! internal energy index
  INTEGER, PARAMETER :: p=6              ! pressure index  
  INTEGER, PARAMETER :: T=7              ! temperature index
  INTEGER, PARAMETER :: c=8              ! sound speed index
  INTEGER, PARAMETER :: mu=9             ! shear viscosity index
  INTEGER, PARAMETER :: bulk=10          ! bulk viscosity index
  INTEGER, PARAMETER :: ktc=11           ! thermal conductivity index
  INTEGER, PARAMETER :: Y=12             ! mass fraction index
  
  INTEGER, PARAMETER :: coordDIM=3
  INTEGER, PARAMETER :: x_c=varDIM+1     ! x-coord index
  INTEGER, PARAMETER :: y_c=varDIM+2     ! y-coord index
  INTEGER, PARAMETER :: z_c=varDIM+3     ! z-coord index

  INTEGER, PARAMETER :: DIM=varDIM+coordDIM

  INTEGER :: nx,ny,nz
  INTEGER :: t1,tf,tloop
  CHARACTER(LEN=flen) :: jobdir

  INTEGER :: ax,ay,az
  INTEGER :: px,py,pz
  
  INTEGER :: proc
  INTEGER, DIMENSION(:,:), ALLOCATABLE :: procmap
  INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: invprocmap
  
  DOUBLE PRECISION, PARAMETER :: pi=3.14159265358979323846d0,zero=0.0d0,one=1.0d0,two=2.0d0

END MODULE globals

