!----------------------------------------
!----------------------------------------
!  Britton J. Olson
!  Stanford University
!  Department of Aero/Astro
!----------------------------------------
!----------------------------------------
!  This program will generate the structured grid
!  boundaries for the experimental nozzle setup
!  as published by Papamoschou AIAA 2004.
MODULE inputs
DOUBLE PRECISION , PARAMETER :: zero=0.0d0,one=1.0d0,two=2.0d0,three=3.0d0,half=.5d0,four=4.0d0
DOUBLE PRECISION , PARAMETER :: pi=3.141592653589793238462643383279 

INTEGER, PARAMETER :: seg1=10
INTEGER, PARAMETER :: seg2=100
INTEGER, PARAMETER :: seg3=30
INTEGER, PARAMETER :: seg4=25
INTEGER, PARAMETER :: seg5=100
INTEGER, PARAMETER :: seg6=2
INTEGER, PARAMETER :: seg7=10
INTEGER, PARAMETER :: seg=seg1+seg2+seg3+seg4+seg5+seg6+seg7

!  Specify Units in [mm]
DOUBLE PRECISION , PARAMETER :: throat = 17.8d0            ! Throat height
DOUBLE PRECISION , PARAMETER :: Aratio = 1.6               ! Area-ratio of the nozzle
DOUBLE PRECISION , PARAMETER :: NPR = 1.70d0                ! Nozzle Pressure Ratio (P_res/P_a) [1.2,1.8]
DOUBLE PRECISION , PARAMETER :: gam = 1.4d0                ! Ratio of specific heats
DOUBLE PRECISION , PARAMETER :: len = 117.0d0              ! Length of nozzle (throat to exit)
DOUBLE PRECISION , PARAMETER :: noz_x = 0.0d0              ! X-location of throat
DOUBLE PRECISION , PARAMETER :: noz_y = throat/two         ! Y-location of throat
DOUBLE PRECISION , PARAMETER :: inlet_x = 25.0d0          ! X-location of inlet (relative-up from noz_x)
DOUBLE PRECISION , PARAMETER :: inlet_y = 0.0d0           ! Y-location of inlet (relative-left from noz_y)

!  Segment parameters - Outer Boundaries
LOGICAL                      :: full = .TRUE.              ! Make full nozzle or half
DOUBLE PRECISION , PARAMETER :: R3 = 0.25d0*noz_y            ! Radius of turn for segment #3
DOUBLE PRECISION , PARAMETER :: R5=1.5d1*Aratio*two*noz_y  ! Radius of background mesh

!  Grid Parameters
INTEGER , PARAMETER          :: nx = 768                    ! Number of points in /xi direction
INTEGER , PARAMETER          :: ny = 256                    ! Number of points in /eta direction
DOUBLE PRECISION             :: wall  = 0.08d0              ! Zoom factor at wall in uniform spacing (overwritten by wall_raw)
DOUBLE PRECISION , PARAMETER :: wall_raw  = 7.0d-3          ! Raw spacing in [mm] .. this value will set wall
LOGICAL, PARAMETER           :: rwall_on = .TRUE.           ! Use raw to set wall? or not?
DOUBLE PRECISION , PARAMETER :: wall2 = .2d0 !0.03d0        ! Thickness of refined region near wall  (This does nothing 3/2011 BJO)
DOUBLE PRECISION , PARAMETER :: wall3 = 10.0d0              ! Width of interface from fine to course in grid points  (ALSO does nothing)
DOUBLE PRECISION , PARAMETER :: core  = 0.9d0               ! Zoom factor of the plume in uniform spacing units
DOUBLE PRECISION , PARAMETER :: core2 = pi/4.0d0            ! Thickness of refined region in the plume in [radians]
DOUBLE PRECISION , PARAMETER :: core3 = 10.0d0              ! Width of interface from fine to course in grid points

DOUBLE PRECISION , PARAMETER :: loci = dble(nx)*.9d0  ! Location in grid number of stretch-> iso interface
DOUBLE PRECISION , PARAMETER :: thick = 32.0d0              ! Thickness of transition region above

!  Quasi-1D initialization 
INTEGER , PARAMETER          :: np = 200                   ! Number of point in 1-D direction
DOUBLE PRECISION , PARAMETER :: offset = 5.2d0
DOUBLE PRECISION , PARAMETER :: Xstt = noz_x-inlet_x       ! X-location to start 1D analysis
DOUBLE PRECISION , PARAMETER :: Xend = noz_x+len           ! X-location to stop 1D analysis

!  Flags
INTEGER :: init

END MODULE


!___________________________
!!!!!!!!!!!!!!!!!!!!!!!!!!!|
!!  MAIN PROGRAM  !!!!!!!!!|
!!!!!!!!!!!!!!!!!!!!!!!!!!!|
!---------------------------
PROGRAM noz_gridgen
USE inputs
CHARACTER * 100 :: BUFFER

  funit=1
  CALL GETARG(1,BUFFER)
  READ(BUFFER,*) init
  
  IF (init==0) THEN
     CALL bounds()       ! Get the outer boundary points
     CALL stretch()      ! Get the mapped grid coordinates
     CALL write_prm()    ! Write the GG input file
     CALL oneD_init()    ! Solve the [M,P,rho,A] = F( x ), quasi-1D flow for given Pressure ratio for use flow initialization
     PRINT*,'Done with grid setup files...'
  ELSE
     CALL init_field()
  END IF

END PROGRAM noz_gridgen
!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!
!--------------------------





SUBROUTINE bounds()
USE inputs
IMPLICIT NONE
DOUBLE PRECISION, DIMENSION(seg1) :: x1,y1
DOUBLE PRECISION, DIMENSION(seg2) :: x2,y2
DOUBLE PRECISION, DIMENSION(seg3) :: x3,y3,theta3
DOUBLE PRECISION, DIMENSION(seg4) :: x4,y4
DOUBLE PRECISION, DIMENSION(seg5) :: x5,y5,theta5
DOUBLE PRECISION, DIMENSION(seg6) :: x6,y6
DOUBLE PRECISION, DIMENSION(seg7) :: x7,y7
DOUBLE PRECISION, DIMENSION(seg1+seg2+seg3+seg4+seg5+seg6+seg7) :: x,y

INTEGER :: i,funit,count
DOUBLE PRECISION :: x10,x1f,y10,y1f
DOUBLE PRECISION :: x20,x2f,y20,y2f
DOUBLE PRECISION :: x30,x3c,y30,y3c,dydx3,theta30,theta3f
DOUBLE PRECISION :: x40,x4f,y40,y4f
DOUBLE PRECISION :: x5c,y5c,theta50,theta5f
DOUBLE PRECISION :: x60,x6f,y60
DOUBLE PRECISION :: x70,y70,y7f
DOUBLE PRECISION :: c5,c4,c3,c2,xrel,d2T


!!!!!!!!!!!!!!!!!!!!!!! --- Diagram of the Nozzle segments --- !!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                         |--------______              !!
!!                                                         |               \            !!
!!                                                 --#4--  |                \  --#5--   !!
!!                                                         |                 \          !!
!!-------__  --#1--                                        |                  \         !!
!!         --_                                            -                    \        !!
!!            -_                  --#2-- _______________--  --#3--              \       !!
!!              - _         ______-------                                        \      !!
!!                 -_  _ ---                                                      |     !!
!! --#7--            --                                                           |     !!
!!                                                                                |     !!
!!________________________________________________________________________________|     !!
!!                            --#6--                                                    !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!---  Segment #1 ---!!!
!! 5th-Order Polynomial with 2 dirichleit & 2 nuemann BC and 2 2nd derivatives...
d2T = three*throat/two*(Aratio-one )/len**2
x10 = noz_x - inlet_x
x1f = noz_x
y10 = noz_y + inlet_y
y1f = noz_y

!  These coefficients match 1st der. only and set 2nd der to zero at BC
c5 = - 6.0d0*inlet_y/inlet_x**5
c4 = -15.0d0*inlet_y/inlet_x**4
c3 = -10.0d0*inlet_y/inlet_x**3
c2 = zero

!  Same but match 2nd derivative with nozzle and other 2nd der. is zero
!c5 = (d2T * inlet_x**2 - 12.0d0*inlet_y )/ (two*inlet_x**5)
!c4 = (3.0d0*d2T*inlet_x**2-30.0d0*inlet_y)/(two*inlet_x**4)
!c3 = (3.0d0*d2T*inlet_x**2-20.0d0*inlet_y)/(two*inlet_x**3)
!c2 =  d2T / two

!  Match 1,2,3 derivatives and relax 2nd der. at other BC
!c5 = -(d2T * inlet_x**2-four*inlet_y)/inlet_x**5
!c4 = -(three*d2T*inlet_x**2-10.0d0*inlet_y)/(two*inlet_x**4)
!c3 = 0.0d0
!c2 = d2T / two


funit=1
OPEN(UNIT=funit,FILE='seg.dat',STATUS='UNKNOWN')
WRITE(funit,*) x10,y10,1
DO i=1,seg1
   x1(i) = (x1f-x10)/dble(seg1)*dble(i) + x10
   xrel = x1(i) - noz_x
   y1(i) = c5*xrel**5 + c4*xrel**4 + c3*xrel**3 +c2*xrel**2
   y1(i) = y1(i) + y1f
   WRITE(funit,*) x1(i),y1(i)
END DO

!!!---  Segment #2 ---!!!
!! 3rd-Order Polynomial with 2 dirichleit & 1 nuemann BC & 1 2nd derivative
!! This is the solution to the cantalever beam problem
x20 = noz_x
x2f = noz_x + len
y20 = noz_y 
y2f = noz_y*Aratio 

DO i=1,seg2
   x2(i) = (x2f-x20)/dble(seg2)*dble(i) + x20
   y2(i) = (-(x2(i) - x2f)* (x2(i)**2 - three* x20**2 - two* x2(i)* x2f + 6.0d0* x20* x2f - & 
        & two* x2f**2)* y20 + (x2(i) - x20)**2 *(x2(i) + two* x20 - three* x2f)* y2f)/(two* (x20 - x2f)**3)
   WRITE(funit,*) x2(i),y2(i)
END DO

!!!---  Segment #3 ---!!!
!! Simple Radial-arc of radius -R- and matching derivative and location from cantalever beam
dydx3 = three/two * (y2f-y20)/(x2f-x20)
x30 = noz_x + len
y30 = y2(seg2)
x3c = -(dydx3**2* R3**2)/sqrt(dydx3**2 *(one + dydx3**2)* R3**2) + x30
y3c = sqrt(R3**2/(one + dydx3**2)) + y30 
theta30 = -pi/two + abs(atan(dydx3))
theta3f = zero
  
DO i=1,seg3
   theta3(i) = (theta3f-theta30)/dble(seg3)*dble(i)+theta30
   x3(i) = R3*cos(theta3(i)) + x3c
   y3(i) = R3*sin(theta3(i)) + y3c
   WRITE(funit,*) x3(i),y3(i)
END DO

!!!---  Segment #4 ---!!!
!! Simple verticle Line
x40 = x3(seg3)
x4f = x40
y40 = y3(seg3) 
y4f = R5

DO i=1,seg4
   x4(i) = x40
   y4(i) = (R5-y40)/dble(seg4) * dble(i) + y40
   IF (i==seg4) THEN    
      WRITE(funit,*) x4(i),y4(i),1
   ELSE
      WRITE(funit,*) x4(i),y4(i)
   END IF
END DO

!!!---  Segment #5 ---!!!
!! Radial Sweep for back ground outflow
x5c = x40
y5c = zero
theta50 = pi/two 
theta5f = zero

DO i=1,seg5
   theta5(i) = (theta5f - theta50)/dble(seg5) * dble(i) + theta50
   x5(i) = R5*cos(theta5(i)) + x5c
   y5(i) = R5*sin(theta5(i)) + y5c
   IF (i==seg5 .and. full .eq. .FALSE.) THEN
      WRITE(funit,*) x5(i),y5(i),1
   ELSE
      WRITE(funit,*) x5(i),y5(i)
   END IF
END DO


!!!---  Reflect the segments and run backwards  ---!!! 
IF(full) THEN
   DO i=seg5-1,1,-1
      WRITE(funit,*) x5(i),-y5(i)
   END DO
   
   DO i=seg4,1,-1
      IF (i==seg4) THEN    
         WRITE(funit,*) x4(i),-y4(i),1
      ELSE
         WRITE(funit,*) x4(i),-y4(i)
      END IF
   END DO
 
   DO i=seg3,1,-1
      WRITE(funit,*) x3(i),-y3(i)
   END DO

   DO i=seg2,1,-1
      WRITE(funit,*) x2(i),-y2(i)
   END DO

   DO i=seg1,1,-1
      WRITE(funit,*) x1(i),-y1(i)
   END DO
END IF

!!!---  Segment #6 ---!!!
!! Straight line back to  inlet, only call if half-nozzle
x60 = x5(seg5)
x6f = x10
y60 = zero 

DO i=1,seg6
   x6(i) = (x6f-x60)/dble(seg6)*dble(i) + x60
   y6(i) = y60
   IF(full .NEQV. .TRUE.) THEN
      IF (i==seg6) THEN
         WRITE(funit,*) x6(i),y6(i),1
      ELSE
         WRITE(funit,*) x6(i),y6(i)
      END IF
   END IF
END DO

!!!---  Segment #7 ---!!!
!! Straight line for inlet
x70 = x6(seg6)
y70 = zero 
y7f = y10

DO i=1,seg7-1
   x7(i) = x70
   y7(i) = (y7f - y70)/dble(seg7)* dble(i) + y70
   IF (full .neqv. .TRUE.) THEN
      WRITE(funit,*) x7(i),y7(i)
   END IF
END DO


IF(full) THEN
   WRITE(funit,*) x70,-y7f,1
   DO i=seg7-1,1,-1
      WRITE(funit,*) x7(i),-y7(i)
   END DO
   WRITE(funit,*) x70, zero
   DO i=1,seg7-1
      WRITE(funit,*) x7(i),y7(i)
   END DO
END IF

CLOSE(UNIT=funit)

END SUBROUTINE bounds

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


SUBROUTINE stretch()
USE inputs
IMPLICIT NONE
DOUBLE PRECISION, DIMENSION (nx,ny) :: xb,yb
DOUBLE PRECISION, DIMENSION (ny-1) :: dyV,dyV2
DOUBLE PRECISION :: dyU,dymin,dymax,dymax2,ytmp
DOUBLE PRECISION :: C1,C3,L0,nn,jj,sig
DOUBLE PRECISION :: C1p,C3p,ybc,ybw,aa,bb, sech,widef
DOUBLE PRECISION :: blend,wide,shift,dy,sumy,sum1,sum2,wide2,dwide,ds_dw,eps
INTEGER :: i,j,funit

IF(rwall_on) THEN
   wall = wall_raw / ( throat / dble(ny) ) 
END IF


L0 = one    ! Make mesh on the unit square

!! Newton Rhapson Parameters
dwide = 1.0d-4
eps = 1.0d-12

!! Hyperbolic Tangent parameters

! Boundary Layer stretch dy = tanh(f1(j)) + tanh(f2(j))
!!$dymin = wall*L0/dble(ny-1)
!!$dymax = 2.0*L0/dble(ny-1)
!!$wide = wall3/2.0d0
!!$shift = (wall2)*dble(ny-1)/wall/2.0d0

!!$sumy = one
!!$DO WHILE ( abs(sumy) > eps)
!!$   dymax2 = dymax + dwide
!!$   CALL sumDY(wide,dymin,dymax,ny,dyV,sum1,shift)
!!$   CALL sumDY(wide,dymin,dymax2,ny,dyV,sum2,shift)
!!$   sum1 = sum1 - L0 + eps
!!$   sum2 = sum2 - L0 + eps
!!$   ds_dw = (sum2-sum1)/dwide   
!!$   dymax = dymax - sum1/ds_dw
!!$   sumy = sum1
!!$END DO
!!$CALL sumDY(wide,dymin,dymax,ny,dyV,sumy,shift)

! BL stretch...y = tanh(f(j))
sumy = one
aa = dble(ny+1)/two
bb = dble(ny) - aa
dymin = wall*L0/dble(ny-1)*two 
wide = ny
DO WHILE ( abs(sumy) > eps)
   wide2 = wide + dwide
   sum1 = wide - (sech(bb/wide))**2 / (dymin*tanh(bb/wide))
   sum2 = wide2 - (sech(bb/wide2))**2 / (dymin*tanh(bb/wide2))
   ds_dw = (sum2-sum1)/dwide   
   wide = wide - sum1/ds_dw
   sumy = sum1
END DO
widef = wide


! Exit Plume...focus grid at the core
dymin = core*L0/dble(ny-1)
dymax = 2.0*L0/dble(ny-1)
wide = core3/2.0d0
shift =  dble(ny/4) - (core2/pi)*dble(ny-1)/core/4.0d0

sumy = one
DO WHILE ( abs(sumy) > eps)
   dymax2 = dymax + dwide
   CALL sumDY(wide, dymax,dymin,ny,dyV2,sum1,shift)
   CALL sumDY(wide,dymax2,dymin,ny,dyV2,sum2,shift)
   sum1 = sum1 - L0 + eps
   sum2 = sum2 - L0 + eps
   ds_dw = (sum2-sum1)/dwide
   dymax = dymax - sum1/ds_dw     
   sumy = sum1
END DO
CALL sumDY(wide,dymax,dymin,ny,dyV2,sumy,shift)


OPEN(UNIT=funit,FILE='stretch.dat',STATUS='UNKNOWN')
ybw = zero
ybc = zero
yb(:,1) = zero
DO j=1,ny
   DO i=1,nx
      blend = tanh( (dble(i)-loci)/thick)
      blend = half * (blend + one)
      ! Uniform in X
      xb(i,j) = dble(i-1)*L0/dble(nx-1)
      !  Hyperbolic-Tangent
      if (j>1) then
         !ybw = ytmp + dyV(j-1)
         ybw =  (one+ tanh( (dble(j)-aa)/widef) / abs(tanh(bb/widef)) ) / two
         ybc = ytmp + dyV2(j-1)
      else
         ybw = zero
         ybc = zero
      end if
      yb(i,j) = ybw *(one-blend) + blend*ybc 
      WRITE(funit,*) xb(i,j), yb(i,j)
   END DO
   ytmp = ybc
END DO

xb(:,1) = xb(:,2)

CLOSE(funit)

END SUBROUTINE stretch


SUBROUTINE write_prm()

!  Write the file here... one less thing to cart around
OPEN(UNIT=1,FILE='noz.prm',STATUS='UNKNOWN')
WRITE(1,*) 'input seg.dat'
WRITE(1,*) 'output nozzle.grid'
WRITE(1,*) 'grid stretch.dat'
WRITE(1,*) 'nppe 0'
WRITE(1,*) 'nnnodes 10'
WRITE(1,*) 'precision 1.0e-11'
WRITE(1,*) 'newton 1'
WRITE(1,*) 'rectangle rect.noz'
WRITE(1,*) 'sigmas sigmas.noz'
CLOSE(1)

END SUBROUTINE write_prm



SUBROUTINE sumDY(wide,dymin,dymax,ny,dy,sumy,shift)
IMPLICIT NONE
DOUBLE PRECISION :: shift,sumy,dymin,dymax,wide
DOUBLE PRECISION, PARAMETER :: one=1.0d0, two=2.0d0
INTEGER :: j,ny
DOUBLE PRECISION, DIMENSION (ny-1) :: dy

!!  Hyperbolic Tangent
wide = wide
sumy = 0.0d0
DO j=1,ny-1
   dy(j) = dymin +   (dymax-dymin) * (one + tanh( dble(j-1)/wide - shift ))/two
   dy(j) = dymin +   (dy(j)-dymin) * (one - tanh( (dble(j - (ny-1) ) )/wide + shift ))/two
   sumy = sumy + dy(j)
END DO

END SUBROUTINE sumDY

!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE oneD_init()
USE inputs
IMPLICIT NONE
DOUBLE PRECISION, DIMENSION(np) :: Mach,P,rho,A,x,h
DOUBLE PRECISION, DIMENSION(seg1+seg2+seg3+seg4) :: xb,yb
DOUBLE PRECISION :: Xst,Ast,fMach,gp1,gm1
DOUBLE PRECISION :: M,M1,M2,dM,res,dv_dm,fun1,fun2
DOUBLE PRECISION :: A_jmp,P_jmp,rho_jmp,Pt_2
INTEGER :: funit,i,j,topL,jm,s_loc,max

funit=1
topL = seg1 + seg2 + seg3 + seg4
OPEN(UNIT=funit,FILE='seg.dat',STATUS='UNKNOWN')
DO i=1,topL
   READ(funit,*) xb(i),yb(i)
END DO
CLOSE(funit)

! Guess shock location 
s_loc = np*3/4

! Check Xstt to make sure its not too small
Xst = Xstt
IF (Xstt < xb(1) ) Xst = xb(1)

! Set x to be uniform
DO i=1,np
   x(i) = (Xend - Xstt) / ( dble(np-1) )*dble(i-1) + Xstt
END DO

! Interpolate height using y_b (simple linear)
CALL interp(xb,x,yb,h,size(xb),np)

gp1 = gam + one
gm1 = gam - one

! Get Area ratios
Ast = throat
!A = two * h / Ast 
A = two * h / Ast + (two*h - Ast)/Ast*(offset-one)
! Solve for Profiles
! Pressure is in units of P_amb
max = 1
DO WHILE(max < np)
   max = max + 1

   DO i = 1 , np
      
      dM = 1.0d-5
      IF ( x(i) < noz_x ) M = dM
      IF ( x(i) < noz_x .and. i > 1) M = Mach(i-1)
      IF ( x(i) > noz_x ) M = 1.1d0
      Mach(i) = fMach(A(i),gam,M)
      
      ! Find Pressure
      P(i) = (one+(gam-one)*Mach(i)**2/two)**(-gam/gm1) ! In ambient units
      rho(i) = (one + gm1/two * Mach(i)**2)**(-one/gm1)     ! In res units

   END DO

   ! Apply Jump Relations
   M1 = Mach(s_loc)
   M2 = sqrt(  ( gm1*M1**2 + two) / (two*gam*M1**2 - gm1)  )  !  Post-shock Mach Number
   Pt_2 = (gp1*M1**2/(gm1*M1**2+two))**(gam/gm1) * (gp1/(two*gam*M1**2-gm1))**(one/gm1)  !  Post-shock total pressure
   
   A_jmp = Pt_2
   !A = two*h / Ast
   A = two * h / Ast + (two*h - Ast)/Ast*(offset-one)
   A(s_loc:np) = A(s_loc:np) * A_jmp
   
   ! Post-shock region, M<1
   M = dM
   DO i=s_loc,np
      Mach(i) = fMach(A(i),gam,M)
      
      ! Find Pressure
      P(i) = Pt_2*(one+(gam-one)*Mach(i)**2/two)**(-gam/(gam-one))
      rho(i) = Pt_2*(one + (gam-one)/two * Mach(i)**2)**(-one/(gam-one))
   END DO

   ! Check ambient pressure and move shock location accordingly
   IF (P(np) < 1.*one/NPR) s_loc = s_loc - 1
   IF (P(np) > 1.*one/NPR) s_loc = s_loc + 1
   
END DO

OPEN(UNIT=funit,FILE='profile.dat',STATUS='UNKNOWN')
A  = h 
DO i=1,np
   WRITE(funit,'(5ES12.4,I5)') x(i),mach(i),p(i),rho(i),A(i)
END DO
CLOSE(funit)

PRINT*,'Done with 1d profile setup...'

END SUBROUTINE oneD_init

FUNCTION fMach(A,gamma,sup)
IMPLICIT NONE
DOUBLE PRECISION :: fMach,A,gamma,M,M2,dM,res,fun1,fun2,dv_dm,sup,Astar
DOUBLE PRECISION, PARAMETER :: acur=1.0d-3
!sup = .1d0    ! Subsonic
!sup = 1.1d0   ! Supersonic
M = sup
res = 1.0d0
dM = 1.0d-5

DO WHILE ( abs(res) > 1.0d-8)
  M2 = M + dM
  fun1 = Astar(M ,gamma) - A
  fun2 = Astar(M2,gamma) - A
  dv_dm = (fun2-fun1)/dM
  !IF ( fun2 == fun1 ) THEN
  !  M = M + dM/2.0d0
  !  res = 0.0d0
  !ELSE
    M = M - fun1/dv_dm
    res = abs(fun1)
  !  dM = dM/2.0d0
  !END IF      
  !pr
END DO

fMach = M

END FUNCTION


FUNCTION Astar(M,gamma)
IMPLICIT NONE
DOUBLE PRECISION :: M,gamma,gm1,gp1,Astar

gm1 = gamma - 1.0d0
gp1 = gamma + 1.0d0
Astar = 1.0d0/M * ( 2.0d0/gp1 + gm1/gp1*M**2 )**(gp1/(2.0d0*gm1))

END FUNCTION



SUBROUTINE init_field
USE inputs
IMPLICIT NONE
DOUBLE PRECISION, DIMENSION(nx,ny) :: x_c,y_c,u,v,p,rho,M,dx,dy,plume
DOUBLE PRECISION, DIMENSION(nx) :: tmp,tmp1d
DOUBLE PRECISION, DIMENSION(ny) :: tmpy
DOUBLE PRECISION, DIMENSION(np) :: x_1d,Mach,p_1d,rho_1d
DOUBLE PRECISION :: bL
INTEGER :: funit,i,j

x_c = zero;y_c = zero;
OPEN(UNIT=funit,FILE='nozzle.grid',STATUS='OLD')
DO j=1,ny
   DO i=1,nx
      READ(funit,*) x_c(i,j), y_c(i,j)
   END DO
END DO
CLOSE(funit)
print*,'Done reading grid...making initialization file'

OPEN(UNIT=funit,FILE='profile.dat',STATUS='OLD')
DO i=1,np
   READ(funit,*) x_1d(i),Mach(i),p_1d(i),rho_1d(i)
END DO
CLOSE(funit)

! Solve for initial field, pressure, density
P = one/NPR
rho = rho_1d(np)
u = zero
v = zero
M = zero
DO i=1,nx
   DO j=1,ny
      plume(i,j) = (one+tanh( (x_c(i,j)-1.d0*Xend)/(noz_y/10.d0)  ))/two - (one-(one+tanh( (abs(y_c(i,j))-Aratio*noz_y)/(noz_y/5.d0)  )  )/two)
      plume(i,j) = plume(i,j) + (one+tanh((x_c(i,j)-Xend*1.3d0)/noz_y))/two
   END DO
END DO

WHERE (plume>one) plume=one
WHERE (plume<zero) plume=zero

DO j=1,ny
   tmp = (one + tanh( (x_c(:,j)-Xend)/(throat/10.d0)))/two   
   
   CALL interp(x_1d,x_c(:,j),p_1d,tmp1d,np,nx)
   P(:,j) = tmp1d + plume(:,j)*(P(:,j)-tmp1d)
   
   CALL interp(x_1d,x_c(:,j),rho_1d,tmp1d,np,nx)
   rho(:,j) = tmp1d + plume(:,j)*(rho(:,j)-tmp1d)
   
   CALL interp(x_1d,x_c(:,j),mach,tmp1d,np,nx)
   M(:,j) = tmp1d + plume(:,j)*(M(:,j)-tmp1d)
   
END DO

! Split up Mach number into vectors aligned with the grid
! take some simple derivatives to get orientations
DO i=2,nx-1
   dx(i,:) = ( x_c(i+1,:) - x_c(i-1,:) )/ two
   dy(i,:) = ( y_c(i+1,:) - y_c(i-1,:) )/ two
END DO
dx(1,:) = x_c(2,:)-x_c(1,:)
dy(1,:) = y_c(2,:)-y_c(1,:)
dx(nx,:) = x_c(nx,:)-x_c(nx-1,:)
dy(nx,:) = y_c(nx,:)-y_c(nx-1,:)

u = dx/ sqrt(dx**2 + dy**2) * M
v = dy/ sqrt(dx**2 + dy**2) * M

! Boundary Layer tanh approximation over a few points 
bL = 3.0d-1
DO j=1,ny
   tmpy(j) = (one-(one+tanh((dble(j)-1)*bL))/two) +  (one+tanh((dble(j)-ny)*bL ))/two
   tmpy(j) = one - tmpy(j)
   u(:,j) = u(:,j) * tmpy(j)
   v(:,j) = v(:,j) * tmpy(j)
END DO

OPEN(1,file='nozzle_init.tec',status='unknown') 
WRITE (1,*) ' VARIABLES = "X", "Y", "P","rho","U","V" '
WRITE(1,*) "ZONE I=", nx, ", J=", ny, ", F=POINT"
DO  J=1,ny
   DO  I=1, nx
      WRITE (1,*) x_c(i,j), y_c(i,j), p(I,J), rho(I,J),u(I,J),v(I,J)
   END DO
END DO
CLOSE(1)


OPEN(1,file='nozzle.inflow',status='unknown')
WRITE (1,*) NPR,P_1d(1),rho_1d(1),Mach(1)
CLOSE(1)

END SUBROUTINE init_field


SUBROUTINE interp(X1,X2,Y1,Y2,N1,N2)
IMPLICIT NONE
INTEGER :: N1,N2
DOUBLE PRECISION, DIMENSION(N1) :: X1,Y1
DOUBLE PRECISION, DIMENSION(N2) :: X2,Y2
INTEGER :: j,jm,i

j = 1
jm= 0 
DO i=1,n2
  DO WHILE(1==1)
    IF( x2(i) == x1(j) ) THEN  ! Are we dead on ?
      y2(i) = y1(j)
      EXIT
    ELSEIF( j == n1) THEN  ! Were at the end and still havent found a bounds.... just give it last point
      y2(i) = y1(n1)
      EXIT
    ELSEIF( x1(j) < x2(i) ) THEN  ! Are we too small
      j = j + 1
      jm= j - 1
    ELSEIF( x1(j) > x2(i) .and. x1(jm) < x2(i)  ) THEN ! In range... interpolate
      jm = j - 1
      y2(i) = y1(jm) + ( y1(j) - y1(jm))/(x1(j)-x1(jm))*(x2(i)-x1(jm))
      EXIT
    ELSE  ! Too far right... move left one
      j = j - 1
      jm= j - 1
    END IF
  END DO
END DO

END SUBROUTINE


FUNCTION sech(x)
DOUBLE PRECISION :: x,sech
sech = 1.0d0 / cosh(x)
END FUNCTION 
