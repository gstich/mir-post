!! Britton J. Olson
!! Pre-processor for Delery-bump grid generation

MODULE inputs
  IMPLICIT NONE
  DOUBLE PRECISION, PARAMETER :: Lx = 3.0
  DOUBLE PRECISION, PARAMETER :: Ly = 1.5
  DOUBLE PRECISION, PARAMETER :: Lb = 0.2
  DOUBLE PRECISION, PARAMETER :: Hb = 0.02
  DOUBLE PRECISION, PARAMETER :: Db = 1.5

  INTEGER, PARAMETER :: nx = 512
  INTEGER, PARAMETER :: ny = 128

  INTEGER, PARAMETER :: segB = nx

END MODULE inputs



PROGRAM delery_bump
  USE inputs
  IMPLICIT NONE

  DOUBLE PRECISION, DIMENSION(segB) :: x,y
  INTEGER :: funit,i

  !! Get the shape of the bump.  Use a gauss for now...
  x = (/( Lx*dble(i-1)/dble(nx-1), i=1, nx)/)
  y = Hb * exp( - (x-Db)**2 / Lb**2*10  )
  
  funit=17
  OPEN(UNIT=funit,FILE='seg.dat',STATUS='UNKNOWN')


  
  
  WRITE(funit,*) 0.000,Ly,1
  WRITE(funit,*) Lx,Ly,1
  WRITE(funit,*) Lx,0.000,1
  DO i=segB-1,2,-1
     WRITE(funit,*) x(i),y(i)
  END DO
  WRITE(funit,*) 0.000,0.000,1

  CLOSE(funit)

  CALL write_prm()
  CALL stretch()
  

END PROGRAM delery_bump


SUBROUTINE write_prm()
  USE inputs, ONLY: nx,ny
  IMPLICIT NONE
  INTEGER :: funit

  funit = 76
  !  Write the file here... one less thing to cart around
  OPEN(UNIT=funit,FILE='bump.prm',STATUS='UNKNOWN')
  WRITE(funit,*) 'input seg.dat'
  WRITE(funit,*) 'output bump.grid'
  WRITE(funit,*) 'grid stretch.bump'
  !WRITE(funit,*) 'nx',nx
  !WRITE(funit,*) 'ny',ny
  WRITE(funit,*) 'ppe 0'
  WRITE(funit,*) 'nnnodes 10'
  WRITE(funit,*) 'precision 1.0e-11'
  WRITE(funit,*) 'newton 1'
  WRITE(funit,*) 'rectangle rect.bump'
  WRITE(funit,*) 'sigmas sigmas.bump'
  CLOSE(funit)

END SUBROUTINE write_prm




SUBROUTINE stretch()
USE inputs
IMPLICIT NONE
DOUBLE PRECISION, DIMENSION (nx,ny) :: xb,yb,xbt,blend
DOUBLE PRECISION, DIMENSION (ny-1) :: dyV,dyV2
DOUBLE PRECISION :: dyU,dymin,dymax,dymax2,ytmp
DOUBLE PRECISION :: C1,C3,L0,nn,jj,sig
DOUBLE PRECISION :: C1p,C3p,ybc,ybw,aa,bb, sech,widef
DOUBLE PRECISION :: wide,shift,dy,sumy,sum1,sum2,wide2,dwide,ds_dw,eps
DOUBLE PRECISION :: wall,xwall,A,B,dx,Lb_t,thick
INTEGER :: i,j,funit,ny_tmp,i1,i2,iw,locj

L0 = 1.0d0    ! Make mesh on the unit square
wall = .02d0
xwall = .05


!! Newton Rhapson Parameters
dwide = 1.0d-4
eps = 1.0d-12


ny_tmp = ny*2
sumy = 1.0d0
aa = dble(ny_tmp+1)/2.0d0
bb = dble(ny_tmp) - aa
dymin = wall/dble(ny_tmp-1)*2.0d0 
wide = ny_tmp

dwide = 1.0d-4
eps = 1.0d-12

DO WHILE ( abs(sumy) > eps)
   wide2 = wide + dwide
   sum1 = wide - (sech(bb/wide))**2 / (dymin*tanh(bb/wide))
   sum2 = wide2 - (sech(bb/wide2))**2 / (dymin*tanh(bb/wide2))
   ds_dw = (sum2-sum1)/dwide   
   wide = wide - sum1/ds_dw
   sumy = sum1
END DO
widef = wide

!  Stretch in Y
DO j=1,ny
   yb(:,j) = (1.0d0+ tanh( (dble(j)-aa)/widef) / abs(tanh(bb/widef)) ) / 2.0d0
END DO
yb = yb / yb(1,ny)


OPEN(UNIT=funit,FILE='stretch.bump',STATUS='UNKNOWN')

!ybw = 0.0
!ybc = 0.0
yb(:,1) = 0.0
xb(1,:) = 0.0d0
xbt(1,:) = 0.0d0

!! Bump cluster parameters
iw = 20
Lb_t = Lb*1.3d0
i1 = (Lb_t-2.0*Db)*nx*xwall/ ( (Lb_t-2.0*Db)*xwall - Lb_t) / 2.0d0
i2 = nx - i1 

!! Blend parameters
locj = ny*.8
thick = 10.0d0


DO j=1,ny
   DO i=2,nx
      blend(i,j) = tanh( (dble(j-locj) )/thick)
      blend(i,j) = 0.5d0 * (blend(i,j) + 1.0d0)
      ! Uniform in X
      A = .5d0*(1.0d0-tanh( dble(i-i1)/dble(iw)))
      B = .5d0*(1.0d0+tanh( dble(i-i2)/dble(iw)))
      dx = (A+B)*(1.0d0-xwall) + xwall
      ! Focused core grid
      xbt(i,j) = xbt(i-1,j) + dx
      ! Uniform grid
      xb(i,j) = dble(i-1)*L0/dble(nx-1)

   END DO

END DO

xbt = xbt / xbt(nx,1)
xb = xb / xb(nx,1)

xb = xbt*(1.0d0-blend) + xb*blend

DO j=1,ny
   DO i=1,nx
      WRITE(funit,*) xb(i,j), yb(i,j)
   END DO
END DO




!xb(:,1) = xb(:,2)

CLOSE(funit)

END SUBROUTINE stretch

FUNCTION sech(x)
DOUBLE PRECISION :: x,sech
sech = 1.0d0 / cosh(x)
END FUNCTION 
