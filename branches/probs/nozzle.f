! you must re-define this module name to be unique across all problems
MODULE nozzle_data
  SAVE
  ! ---------------------------------------------------------------------------
  ! Place Problem dependent state data here, similar to inputs.f and globals.f
  ! ---------------------------------------------------------------------------  
  CHARACTER(LEN=30) 	:: init_type    = '2D'
  CHARACTER(LEN=90) 	:: restart_dir  = '/path/to/directory'
  CHARACTER(LEN=30) 	:: gridfile    = 'nozzle.grid'
  DOUBLE PRECISION 	:: Pinitial     = 1.01d6 
  DOUBLE PRECISION 	:: Tinitial     = 300.0d0
  ! Not set in Namelist, these variables are read in and used to set the inflow BC
  CHARACTER(LEN=30) 	:: prob_jobname = 'nozzle'
  DOUBLE PRECISION 	:: Mach  = 1.0d0 
  DOUBLE PRECISION      :: P_in   = 1.01d6
  DOUBLE PRECISION      :: rho_in = 1.0d0
  DOUBLE PRECISION      :: U_in = 1000.0d0
  DOUBLE PRECISION      :: T_in = 300.0d0
  DOUBLE PRECISION      :: e_in = 1.0d0      
  DOUBLE PRECISION      :: NPR   = 1.5d0
  REAL,ALLOCATABLE, DIMENSION(:,:)  :: randBS
  !! For the Recycling rescaling routines
  DOUBLE PRECISION 	:: del_BL       = 1.0d-1 
  DOUBLE PRECISION 	:: Re_BL        = 1.0d4
  DOUBLE PRECISION 	:: BLalpha      = 1.0D0      
  DOUBLE PRECISION 	:: mu_0         = 1.0D0      
  DOUBLE PRECISION 	:: Pr           = 0.7D0    
  INTEGER        	:: rcy_pt_g     = 100  
  DOUBLE PRECISION 	:: Len_BL       = 4.0d0
  INTEGER               :: nvar         = 4
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: Qave
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: yp_r,yp_i,yp_o
  DOUBLE PRECISION :: simtime_old
  INTEGER :: res_dump
  INTEGER :: rand_seed = 12
  LOGICAL :: init_stats = .FALSE.     ! Initialize the temporal averages? Will do this if simtime < dt (at startup)
  LOGICAL :: BL_flag = .TRUE.              ! A flag to call routine once in BC

  !! Stuff to save for RR communication
  INTEGER, DIMENSION(:,:,:,:,:), ALLOCATABLE :: MAP
  INTEGER, DIMENSION(:), ALLOCATABLE :: pmapID

  LOGICAL :: pdonor,precv
  INTEGER :: idonor,irecv
  INTEGER :: yproc
  LOGICAL :: flippy
  DOUBLE PRECISION :: tflip,tflow,tiid
  INTEGER :: step_old

  INTEGER :: check = 0
  INTEGER :: filmax = 30
  LOGICAL :: Ifil = .FALSE.

END MODULE nozzle_data

! -----------------------------------------------------------------------------
! prob_inputs
! Read problem-dependent inputs and save in module variables
! must also define 'jobname'
! Called by the Miranda main program.
! -----------------------------------------------------------------------------
SUBROUTINE prob_inputs(fileName)
  USE mpi
  USE prob_interface, ONLY : jobname
  USE nozzle_data
  USE globals, ONLY : inputUnit
  USE inputs, ONLY : nx,ny,nz
  IMPLICIT NONE
  CHARACTER(LEN=*), INTENT(IN) :: fileName
  INTEGER :: nxp,nyp
  INTEGER :: in_unit 
  INTEGER, DIMENSION(1) :: seed

  ! Uncomment to define a namelist
  NAMELIST /nozzle_vars/ Tinitial,Pinitial,del_BL,Re_BL,Pr,rcy_pt_g,init_type,restart_dir,init_stats,BL_flag,&
      & gridfile, Ifil,filmax,BLalpha

  ! Uncomment to open and read a namelist file
  OPEN(UNIT=inputUnit,FILE=TRIM(fileName),FORM='FORMATTED',STATUS='OLD')
  READ(UNIT=inputUnit,NML=nozzle_vars)
  CLOSE(inputUnit)
  
  ! give the right name
  jobname = TRIM(prob_jobname)
  
  ! Read the inflow boundary conditions                                     
  in_unit = 26
  OPEN(UNIT=in_unit, file='nozzle.inflow')
  READ(in_unit,*) NPR,P_in,rho_in,Mach
  CLOSE(in_unit)

  ! Set some time independent random numbers
  ALLOCATE(randBS(2,nz))
  seed = 12
  CALL random_seed(put=seed)
  CALL random_number(randBS)
  
END SUBROUTINE prob_inputs



! -------------------------------------------------------------------------------------
! Set up problem geometry. (This is called on startup and restart before timestepping.)
! -------------------------------------------------------------------------------------
 SUBROUTINE prob_setup()
  USE mpi
  USE constants
  USE inputs
  USE globals
  USE extend
  USE nozzle_data
  USE metrics, ONLY: x_c
  USE interfaces, ONLY : spherecoord, ERF, ghost, SUMproc, ran1,set_time_seed
  IMPLICIT NONE
  DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE :: Gplane,Lplane
  DOUBLE PRECISION :: Rgas,delX
  DOUBLE PRECISION, DIMENSION(1) :: randT
  INTEGER :: i

  ! Some prob constants which are input dependent
  Rgas = Runiv/molwts(1)                        ! Specific Gas constant
  rho_in = rho_in*NPR*Pinitial/(Tinitial*Rgas)  ! Density of inflow- Convert to physical units 
  P_in = P_in*NPR*Pinitial                      ! Pressure of inflow- Convert to physical units
  U_in = Mach*sqrt(P_in*gamma/rho_in)           ! Get the velocity on the inflow
  e_in = (P_in/(gamma-one))/rho_in              ! Energy of inflow- Convert to physical units 
  T_in = P_in / (rho_in * Rgas)                 ! Temp. of inflow- Convert to physical units 
  mu_0 = U_in * rho_in * del_BL / Re_BL         ! Physical viscosity based on inlet parameters
  

  ! Set-up the grid and metric terms for the calculation
  CALL prob_geom()
  check = 0

  !! Get Boundary Layer length from lower right hand corner
  delX = 0.0D0
  IF(ix(1)==1 .AND. iy(1)==1 .AND. iz(1)==1) delX = ABS( x_c(2,1,1)-x_c(1,1,1) )
  delX = SUMproc(delX)
  Len_BL = dble(rcy_pt_g) * delX


  !! Set up the RR mapping routines here.
  CALL yprocBL()
  ALLOCATE(MAP(2,2,pz,yproc,3))
  ALLOCATE(pmapID(4))
  ALLOCATE(Lplane(nvar,ay,az))
  ALLOCATE(Gplane(nvar,ay*yproc,az))

  idonor = rcy_pt_g
  irecv  = 1
  CALL setup_maps(MAP,pmapID,pdonor,idonor,precv,irecv,yproc)

  !! RR boundary layer
  !flippy = .TRUE.                      ! Set the flip flag
  flippy = .FALSE.                      ! Set the flip flag
  tflow = Len_BL/U_in                  ! Flow through time on the BL
  tiid = .25d0                         ! Range in Tflow units of IID RN
  randT = 0.0D0                        ! Init the RN
  CALL set_time_seed(1)                ! Clock based RN seed
  IF (xyzcom_id==0) CALL ran1(randT)    ! Get the RN on master
  randT(1) = SUMproc(randT(1))               ! Send procs the RN
  tflip = tflow*(one-tiid*(randT(1)-one)) ! Get next random time period
  tflip = tflip + simtime              ! Get the absolute value of flip time
  
  IF(xyzcom_id==0) print*,tflow,Len_BL

  !! Print out the map
  IF(world_id == master ) THEN
      !print*,'Yproc',yproc

      print*,'---Inflow---'
      DO i=1,yproc
         print*,MAP(1,2,:,i,1)
      END DO

      DO i=yproc,1,-1
         print*,MAP(1,1,:,i,1)
      END DO
      print*,'-----------'

      print*,'---Donor---'
      DO i=1,yproc
         print*,MAP(2,2,:,i,1)
      END DO

      DO i=yproc,1,-1
         print*,MAP(2,1,:,i,1)
      END DO
      print*,'-----------'

      !print*,xyzcom_id,pmapID
  END IF




 END SUBROUTINE prob_setup




! -----------------------------------------------------------------------------
! prob_init
! Initialize state variables at start of run
! Called by the Miranda main program.
! -----------------------------------------------------------------------------
SUBROUTINE prob_init(rho,u,v,w,e,Y,p,T)
  USE mpi
  USE constants
  USE globals, ONLY: simtime,molwts,ax,ay,az,ix,iy,iz,xloc,yloc,zloc,x1proc,xnproc,y1proc
  USE globals, ONLY: iodata,nres,iodir
  USE interfaces, ONLY: restart,filter,boundary
  USE inputs
  USE prob_interface, ONLY : jobname
  USE nozzle_data
  USE metrics, ONLY : x_c,y_c
  IMPLICIT NONE
  DOUBLE PRECISION, DIMENSION(:,:,:), INTENT(OUT) :: rho,u,v,w,e
  DOUBLE PRECISION, DIMENSION(:,:,:,:), INTENT(OUT) :: Y
  DOUBLE PRECISION, DIMENSION(:,:,:), INTENT(OUT)   :: p,T
  DOUBLE PRECISION, DIMENSION(ax,ay,az) :: dum,tmp0
  DOUBLE PRECISION, DIMENSION(nx,ny) :: tmpP,tmpRho,tmpU,tmpV
  LOGICAL :: noz
  INTEGER :: in_unit,i,j,k,re_unit  
  DOUBLE PRECISION :: Mflow,Pflow,Gam,mwt,thick,x0,y0,jump,dumy
  INTEGER :: nxp,nyp
  CHARACTER(LEN=30) :: comments,filename

  SELECT CASE(init_type)

      CASE('2D')
         !  Read serial ASCII file in .tec format to initialize the flow
         !  This file is generated solving the Quasi-1D Nozzle Area Equations
         !  All procs read in file then keep their chunks (init.tec is 2d)
         filename = 'nozzle.init'
         in_unit = 20
         OPEN(UNIT=in_unit,FILE=filename,STATUS='OLD')
         READ(in_unit,*) comments
         READ(in_unit,*) comments
         DO J=1,ny
            DO I = 1,nx
               READ(in_unit,*) dumy,dumy,tmpP(i,j),tmpRho(i,j),tmpU(i,j),tmpV(i,j)
            END DO
         END DO
         CLOSE(in_unit)

         ! Copy proc data through third direction... perfect extrusion
         DO k=1,nz
            P(:,:,k) = tmpP(ix,iy)*NPR*Pinitial ! Gives P_exit=Pa
            rho(:,:,k) = tmpRho(ix,iy)*NPR*Pinitial/(Tinitial*Runiv/molwts(1))
            u(:,:,k) = tmpU(ix,iy) ! x-dir Mach number
            v(:,:,k) = tmpV(ix,iy) ! y-dir Mach number
         END DO
         !P = Pinitial
         !T = Tinitial
         !rho = Pinitial/(Tinitial*Runiv/molwts(1))


         !! Smooth out the initial field for stability
         IF(Ifil) THEN
            DO i=1,filmax
               CALL filter(gfilter,P,tmp0)
               CALL filter(gfilter,tmp0,P)
               CALL filter(gfilter,rho,tmp0)
               CALL filter(gfilter,tmp0,rho)
               CALL filter(gfilter,u,tmp0)
               CALL filter(gfilter,tmp0,u)
               CALL filter(gfilter,v,tmp0)
               CALL filter(gfilter,tmp0,v)
            END DO
         END IF


         e = (p/(gamma-one))/rho
         T = P / (rho*Runiv/molwts(1))
         u = u*sqrt(P*gamma/rho) ! convert to actual velocity
         v = v*sqrt(P*gamma/rho) ! convert to actual velocity
         w = zero
         Y = one

         !! Add noise in the BL to trip the BL
         IF(BL_flag) THEN
            CALL random_BL(u,v,w)
         END IF

      CASE('3D')
         !  Another option here to read in a restart file as the initialization of the flow field.
         !  The global variable 'iodata' is a pointer to the global flow field variables.  We only need
         !  to set iodata for a full initialization.  A new restart file will be written with the identical 
         !  data.
         WRITE(iodir,'(2A,I4.4)') TRIM(restart_dir) !,'/res',dump
         CALL restart('r',TRIM(iodir),iodata(:,:,:,1:nres))

         DO i=1,nres
            CALL filter(gfilter,iodata(:,:,:,i),iodata(:,:,:,i))
         END DO

     END SELECT



END SUBROUTINE prob_init

! -----------------------------------------------------------------------------
! prob_eos
! Problem specific equation of state call
! Called by routine "eos()" in eos.f in the case when eostype == 'CASES'
! -----------------------------------------------------------------------------
 SUBROUTINE  prob_eos(rho,e,Y,p,T,c,dTdE,kR,kP,dkPdE)
  USE mpi
  USE constants, only : zero
  USE prob_interface, ONLY : jobname
  IMPLICIT NONE
  DOUBLE PRECISION, DIMENSION(:,:,:), INTENT(IN) :: rho,e
  DOUBLE PRECISION, DIMENSION(:,:,:,:), INTENT(IN) :: Y
  DOUBLE PRECISION, DIMENSION(:,:,:), INTENT(INOUT) :: p,T
  DOUBLE PRECISION, DIMENSION(:,:,:), INTENT(OUT), OPTIONAL :: c,dTdE,kR,kP,dkPdE
  ! should never be called for Incompressible problem
  CALL fatal(-1,"prob_eos called for Incompressible problem ",jobname)
  ! give INTENT(OUT) args a value to get rid of annoying compiler warnings
  if (present(c)) c = zero
  if (present(dTde)) dTde = zero
  if (present(kR)) kR = zero
  if (present(kP)) kP = zero
  if (present(dkPdE)) dkPdE = zero
 END SUBROUTINE prob_eos
 
! ----------------------------------------------------------------------------------
! prob_stats
! Problem specific function to compute and write statistics data every "stats"
! timesteps
! Called by routine "monitor()".
! -----------------------------------------------------------------------------------
SUBROUTINE prob_stats(simtime,rho,u,v,w,e,Y,p,T,c)
  use mpi
  USE globals, ONLY : flen

  IMPLICIT NONE
  DOUBLE PRECISION, INTENT(IN) :: simtime
  DOUBLE PRECISION, DIMENSION(:,:,:), INTENT(IN) :: rho,u,v,w,e

  DOUBLE PRECISION, DIMENSION(:,:,:,:), INTENT(IN) :: Y
  DOUBLE PRECISION, DIMENSION(:,:,:), INTENT(IN) :: p,T,c

  INTEGER, PARAMETER :: statsUnit=20
  CHARACTER(LEN=flen) :: statsFile

  ! construct name of stats file
  !WRITE(statsFile,'(2A)') TRIM(jobdir),'/statistics'


END SUBROUTINE prob_stats


! -----------------------------------------------------------------------------------
! prob_plots
! Problem specific function to write x-y plotfile data.  Called just after
! of graphics dump in "viz()" of file "viz.f".
! -----------------------------------------------------------------------------------
SUBROUTINE prob_plots(plotdir)
  USE mpi
  USE globals, ONLY: flen,ax,ay,az,ix,iy,iz,rho,p,u,v,w,Y
  USE inputs, ONLY: nx,ny,nz,z1,dz,gfilter
  USE prob_interface, ONLY : jobname
  USE interfaces, ONLY : SUMproc,subsum3xy,subsum3yz,subsum3xz,filter
  USE nozzle_data
  USE metrics, ONLY: x_c,y_c,z_c
  IMPLICIT NONE
  CHARACTER(LEN=*), INTENT(IN) :: plotdir
  CHARACTER(LEN=flen) :: plotFile
  INTEGER, PARAMETER :: plotUnit=21
  INTEGER :: k
  DOUBLE PRECISION :: x_a,x_b,y_a,y_b,z_a,z_b,core,x0
  DOUBLE PRECISION, DIMENSION(ny) :: ubar,vbar,wbar,rhobar,pbar,ybar
  DOUBLE PRECISION, DIMENSION(nx) :: p_dwn,p_up,p_cen,xbar
  DOUBLE PRECISION, DIMENSION(ax,ay,az) :: rhop,up,vp,wp,pp
  DOUBLE PRECISION, DIMENSION(ax,ay,az) :: tmp0,tmp1,box

  ! Get the x0 value
  x0 = 0.0D0
  IF(ix(1)==1 .AND. iy(1)==1 .AND. iz(1)==1) x0 = x_c(1,1,1)
  x0 = SUMproc(x0)

  ! Mean statistics for incoming BL
  x_a = x0 + 0.9D0 * Len_BL
  x_b = x0 + Len_BL
  z_a = 0.0d0 
  z_b = z_a + DBLE(nz-1)*dz
  
  ! Smooth box will give smoother profiles
  WHERE (x_c >= x_a .and. x_c <= x_b .and. z_c >= z_a .and. z_c <= z_a)
      box = 1.0d0
  ELSEWHERE
      box = 0.0d0
  END WHERE
  !CALL filter(gfilter,box,tmp0)
  !CALL filter(gfilter,tmp0,box)
  !where (box > 1.0d0) box = 1.0d0
  !where (box < 0.0d0) box = 0.0d0

  ubar = SUBSUM3XZ(u*box)/SUBSUM3XZ(box)
  vbar = SUBSUM3XZ(v*box)/SUBSUM3XZ(box)
  wbar = SUBSUM3XZ(w*box)/SUBSUM3XZ(box)
  rhobar = SUBSUM3XZ(rho*box)/SUBSUM3XZ(box)
  pbar = SUBSUM3XZ(p*box)/SUBSUM3XZ(box)
  ybar = SUBSUM3XZ(y_c*box)/SUBSUM3XZ(box)


! master cpu writes plot file--------------------------------------------------------------------
  SELECT CASE(xyzcom_id)
  CASE(master)
    WRITE(plotFile,'(2A)') TRIM(plotdir),'/yprofiles.dat'
    OPEN(UNIT=plotUnit,FILE=TRIM(plotFile),FORM='FORMATTED',STATUS='REPLACE')
    WRITE(plotUnit,*) "%# <1-6> y-loc,ubar,vbar,wbar,rhobar,Pbar"
    DO k=1,ny
      WRITE(plotUnit,'(6ES12.4)') ybar(k),ubar(k),vbar(k),wbar(k),rhobar(k),pbar(k)
    END DO
    CLOSE(plotUnit)
  END SELECT

  ! Pressure down length of nozzle
  z_a = 0.0d0 !! dble(nz-1)*dz / 2.0d0 
  z_b = z_a + dble(nz-1)*dz

  where (z_c >= z_a .and. z_c <= z_a)
      box = 1.0d0
  elsewhere
      box = 0.0d0
  end where
  !CALL filter(gfilter,box,tmp0)
  !CALL filter(gfilter,tmp0,box)
  !where (box > 1.0d0) box = 1.0d0
  !where (box < 0.0d0) box = 0.0d0

  tmp0 = 0.0d0
  if (iy(1)==1) tmp0(:,1,:) = 1.0d0
  tmp1 = box*tmp0

  p_dwn = SUBSUM3YZ(p*tmp1)/SUBSUM3YZ(tmp1)
  xbar = SUBSUM3YZ(x_c*tmp1)/SUBSUM3YZ(tmp1)

  tmp0 = 0.0d0
  if (iy(ay)==ny) tmp0(:,ay,:) = 1.0d0
  tmp1 = box*tmp0
  p_up = SUBSUM3YZ(p*tmp1)/SUBSUM3YZ(tmp1)

  core = .1d0
  where (y_c<=-core .and. y_c>=core) box=0.0d0
  p_cen = SUBSUM3YZ(p*box)/SUBSUM3YZ(box)


  ! master cpu writes plot file--------------------------------------------------------------------
  SELECT CASE(xyzcom_id)
  CASE(master)
    WRITE(plotFile,'(2A)') TRIM(plotdir),'/pressure.dat'
    OPEN(UNIT=plotUnit,FILE=TRIM(plotFile),FORM='FORMATTED',STATUS='REPLACE')
    WRITE(plotUnit,*) "%# <1-4> x-loc,p_lower,p_upper,p_cen"
    DO k=1,nx
      WRITE(plotUnit,'(4ES12.4)') xbar(k),p_dwn(k),p_up(k),p_cen(k)
    END DO
    CLOSE(plotUnit)
  END SELECT


END SUBROUTINE prob_plots

! -----------------------------------------------------------------------------------
! prob_bc
! Set the problem-specific boundary values for the state data
! Note that only problem specific BCs need be set.  If the type of the boundary
! is general, such as WALL, SYMM, SLIP, PFLO, PINF, it has already been set.
! This routine is called from "boundary()" in boundary.f
! -----------------------------------------------------------------------------------
SUBROUTINE prob_bc(rho,u,v,w,e,Y)
  USE mpi
  USE globals, ONLY : p,x1proc,y1proc,z1proc,xnproc,ynproc,znproc,molwts,ax,ay,az,iy,ix,iz
  USE globals, ONLY : simtime,step
  USE inputs, ONLY : nx,ny,nz,gamma,gfilter,lfilter
  USE constants
  USE nozzle_data
  USE metrics
  USE interfaces, ONLY : filter, ran1, SUMproc
  IMPLICIT NONE
  DOUBLE PRECISION, DIMENSION(:,:,:), INTENT(INOUT) :: rho,u,v,w
  DOUBLE PRECISION, DIMENSION(:,:,:), INTENT(INOUT), OPTIONAL :: e
  DOUBLE PRECISION, DIMENSION(:,:,:,:), INTENT(INOUT), OPTIONAL :: Y
  DOUBLE PRECISION, DIMENSION(ax,ay,az) :: dumF,dumT
  DOUBLE PRECISION, DIMENSION(size(rho,1),size(rho,3))  :: tmpw,blowx,blowz,blow
  DOUBLE PRECISION, DIMENSION(SIZE(u,1),SIZE(u,3)) :: dxdA,dydA,dzdA,t_mag,u_mag	!  A-column of Jacobian Tensor
  DOUBLE PRECISION, DIMENSION(1) :: randT
  DOUBLE PRECISION :: del
  DOUBLE PRECISION :: filpt,thick
  INTEGER :: i,k

                                                           
   
      !!!  INFLOW  !!!
      IF (x1proc) THEN
         !  Set the Inflow BC to be constant
         u(1,:,:)    = U_in
         v(1,:,:)    = zero
         w(1,:,:)    = zero
         rho(1,:,:)  = rho_in
         e(1,:,:)    = e_in
      END IF
      

      !! RR boundary layer
      IF(step .NE. step_old) THEN ! Only if we are on 1st step of RK-scheme
         IF( simtime .GT. tflip) THEN          ! Are we beyond the last period?
            flippy = .NOT. flippy                   ! Set the flip flag
            IF (xyzcom_id==0) print*,"Flipping BL ",flippy       ! echos stuff
            randT = 0.0D0                           ! Init the RN
            IF (xyzcom_id==0) CALL ran1(randT)       ! Get the RN on master
            randT(1) = SUMproc(randT(1))            ! Send procs the RN
            tflip = tflow*(one-tiid*(randT(1)-one)) ! Get next random time period
            tflip = tflip + simtime                 ! Get the absolute value of flip time
            CALL MPI_BARRIER(xyzcom,mpierr)
         END IF
      END IF
      step_old = step
      
      CALL inflow_rcy_v2(rho,u,v,w,e)


      !!!  TOP WALL  !!!
      IF (ynproc) THEN

         
         u(:,ay,:) = zero
         v(:,ay,:) = zero
         w(:,ay,:) = zero
         ! Adiabatic Wall
         !rho(:,ay,:) = (54.0d0*rho(:,ay-1,:)-27.0d0*rho(:,ay-2,:)+6.0d0*rho(:,ay-3,:) )/33.0d0
         !e(:,ay,:) = (54.0d0*e(:,ay-1,:)-27.0d0*e(:,ay-2,:)+6.0d0*e(:,ay-3,:) )/33.0d0 
         rho(:,ay,:) = rho(:,ay-1,:) 
         e(:,ay,:)   = e(:,ay-1,:) 



      END IF

      !!!   BOTTOM WALL   !!!
      IF (y1proc) THEN

         u(:,1,:) = zero
         v(:,1,:) = zero
         w(:,1,:) = zero
         
         ! Adiabatic Wall 
         !rho(:,1,:) = (54.0d0*rho(:,2,:)-27.0d0*rho(:,3,:)+6.0d0*rho(:,4,:) )/33.0d0
         !e(:,1,:) = (54.0d0*e(:,2,:)-27.0d0*e(:,3,:)+6.0d0*e(:,4,:) )/33.0d0 
         rho(:,1,:) = rho(:,2,:) 
         e(:,1,:)   = e(:,2,:) 
         

      END IF


      ! Gradually apply the exit boundary conditions  
      filpt = dble(nx-5)
      thick = 3.0d0
      DO i=1,ax
         dumT(i,:,:)=(one+tanh((dble(ix(i))-filpt)/thick))/two
      END DO
      
      !!!  OUT-FLOW  !!!
      IF (xnproc) THEN 
         ! Force back pressure to remain ambient and let density float from NSCBC
         e = e + dumT * ( Pinitial/(gamma-one)/rho - e )                                           
      END IF  
      
      ! Gussian Filter for last N points in x-direction (A)
      CALL filter(gfilter,u,dumF)
      u = u + dumT*(dumF-u) 
      CALL filter(gfilter,v,dumF)
      v = v + dumT*(dumF-v)
      CALL filter(gfilter,w,dumF)
      w = w + dumT*(dumF-w)
      CALL filter(gfilter,e,dumF)
      e = e + dumT*(dumF-e)
      CALL filter(gfilter,rho,dumF)
      rho = rho + dumT*(dumF-rho)
      
      IF(nz==1) w = zero
  
END SUBROUTINE prob_bc



! -----------------------------------------------------------------------------------
! prob_acceleration
! Set the body force source term for a compressible calculation
! called by "acceleration" in "acceleration.f"
! -----------------------------------------------------------------------------------
SUBROUTINE prob_acceleration(simtime,gx,gy,gz)
  USE mpi
  USE prob_interface, ONLY : jobname
  USE inputs, ONLY : compressible,accelx,accely,accelz
  USE constants, ONLY : zero
  IMPLICIT NONE
  DOUBLE PRECISION, INTENT(IN) :: simtime
  DOUBLE PRECISION, DIMENSION(:,:,:), INTENT(OUT) :: gx,gy,gz

  gx = accelx
  gy = accely
  gz = accelz


END SUBROUTINE prob_acceleration


SUBROUTINE prob_esource()

! Do something

END SUBROUTINE prob_esource


! ------------------------------------------------------------------------------
! Set the body force, due to frame acceleration, for incompressible cases.
! ------------------------------------------------------------------------------
 SUBROUTINE prob_Iacceleration(simtime,g)
  USE inputs, ONLY : accelz
  IMPLICIT NONE
  DOUBLE PRECISION, INTENT(IN) :: simtime
  DOUBLE PRECISION, INTENT(OUT) :: g
  g = accelz
 END SUBROUTINE prob_Iacceleration
 
 
 ! ------------------------------------------------------------------------------
! Assign species enthalpies compatible with prob_eos. These are necessary in
! order to compute the enthalpy diffusion term in the energy equation.
! ------------------------------------------------------------------------------
 FUNCTION prob_enthalpies()
  USE constants
  USE inputs, ONLY : ns
  USE globals, ONLY : ax,ay,az,p
  USE leosvars, ONLY : part_r,part_e
  IMPLICIT NONE
  DOUBLE PRECISION, DIMENSION(ax,ay,az,ns) :: prob_enthalpies
  INTEGER :: i
     DO i=1,ns
       prob_enthalpies(:,:,:,i) = part_e(:,:,:,i) + p/part_r(:,:,:,i) ! assumes all species are at the same pressure
     END DO
 END FUNCTION prob_enthalpies

! ------------------------------------------------------------------------------
! Assign problem-specific viscosity, thermal conductivity,
! species diffusivities and magnetic diffusivity.
! ------------------------------------------------------------------------------ 
 SUBROUTINE prob_properties()
  USE constants, ONLY : zero,one,half,three,Runiv
  USE inputs, ONLY : ns,viscous,conductive,diffusive,resistive,diffusive,gamma
  USE globals, ONLY : materials,T
  USE globals, ONLY : molwts
    USE globals, ONLY : mu,bulk,ktc,Diff,eta ! *** OUTPUTS ***
  USE leosvars, ONLY : part_r
  USE nozzle_data
  IMPLICIT NONE
  DOUBLE PRECISION :: T_0,ST,R,cp,cv

  ! Specific Gas Constant
  R = Runiv/molwts(1)
      
  ! Sutherland's Law for Viscosity
  T_0 = 273.15D0  ! Kelvin reference temperature
  ST = 110.4D0    ! Sutherland temperature
  mu = mu_0 * (T/T_0)**(three*half) * (T_0 + ST) / (T + ST)

  ! Bulk viscosity
  bulk = zero
  
  ! Thermal diffusion
  cv = R/(gamma-1.0D0)
  cp = gamma*cv
  ktc = cp * mu / Pr

  IF (ASSOCIATED(Diff)) THEN
    Diff = zero
  END IF
  IF (ASSOCIATED(eta)) THEN
    eta = zero
  END IF
 END SUBROUTINE prob_properties

!-------------------------------------------------------------------------------
! Assign problem-specific opacities. This routine is only called if
! eostype='TABLE' and opacities=2.
! ------------------------------------------------------------------------------ 
 SUBROUTINE prob_opacities()
  USE constants, ONLY : zero,one
  USE inputs, ONLY : ns,radiative
  USE globals, ONLY : materials,rho,p,T,Y
  USE globals, ONLY : atomic_weight,atomic_number,atomic_number23,atomic_number2
  USE globals, ONLY : kP,kR ! *** OUTPUTS ***
  USE leosvars, ONLY : part_r
  IMPLICIT NONE
  DOUBLE PRECISION, DIMENSION(SIZE(Y,1),SIZE(Y,2),SIZE(Y,3)) :: tmp,dum
  DOUBLE PRECISION, DIMENSION(SIZE(Y,1),SIZE(Y,2),SIZE(Y,3),SIZE(Y,4)) :: zeff
  INTEGER :: i,j,k,n
  IF (ASSOCIATED(kP)) THEN
    kP = zero
  END IF
  IF (ASSOCIATED(kR)) THEN
    kR = zero
  END IF
 END SUBROUTINE prob_opacities


!-------------------------------------------------------------------------------
! Assign problem-specific geometry. This routine is only called if
! upon initialization only
! ------------------------------------------------------------------------------ 
SUBROUTINE prob_geom()
 USE constants
 USE metrics, ONLY: x_c,y_c,z_c
 USE nozzle_data
 USE mpi
 USE constants
 USE inputs
 USE globals
 USE interfaces, ONLY : filter
 IMPLICIT NONE
 DOUBLE PRECISION, DIMENSION(ax,ay,az) :: tophat,tmp
 DOUBLE PRECISION, DIMENSION(nx,ny) :: tmpx,tmpy
 INTEGER :: i,j,k,nz_tmp
 INTEGER :: re_unit


  re_unit = 19
  
  
  !  Put file read routine here 
  re_unit = 32
  OPEN(UNIT=re_unit,FILE=gridfile,STATUS='OLD')
  DO J=1,ny
      DO I = 1,nx
         READ(re_unit,*) tmpx(I,J),tmpy(I,J)
      END DO
  END DO
   
  x_c(:,:,1) = tmpx(ix,iy)
  y_c(:,:,1) = tmpy(ix,iy)
  CLOSE(re_unit)

  !  Z-direction extrusion
  DO k=1,az
      z_c(:,:,k) = zcom_id*az*dz + dz*dble(k-1)
      x_c(:,:,k) = x_c(:,:,1)
      y_c(:,:,k) = y_c(:,:,1)
  END DO

  !  For Nozzle-grid convert mm -> cm
  !  Z is left out because it is specified through dz
  x_c = x_c*1.0d-1
  y_c = y_c*1.0d-1


  !! Make sure everyone waits till script if finished 
  CALL MPI_BARRIER(MPI_COMM_WORLD,mpierr)

  ! Make sure Jacobian is 2d
  nz_tmp = nz
  nz = 1
  CALL get_jacobian()
  nz = nz_tmp

END SUBROUTINE




FUNCTION sech(x)
 DOUBLE PRECISION :: x,sech
      sech = 1.0d0 / cosh(x)
END FUNCTION 


SUBROUTINE interpZv(X1,X2,Y1,Y2,n1,n2,nZ,var)
  !! Assumes there is no variation in x1,x2,y1,y2 int the z (2nd) index direction
  IMPLICIT NONE
  INTEGER :: N1,N2,nz,var
  DOUBLE PRECISION, DIMENSION(N1),    INTENT(IN)  :: X1
  DOUBLE PRECISION, DIMENSION(var,N1,nz), INTENT(IN)  :: Y1
  DOUBLE PRECISION, DIMENSION(N2),    INTENT(IN)  :: X2
  DOUBLE PRECISION, DIMENSION(var,N2,nz), INTENT(OUT) :: Y2
  INTEGER :: j,jm,i
  
      j = 1
      jm= 0 
      DO i=1,n2
         DO WHILE(1==1)
            IF( x2(i) == x1(j) ) THEN ! Are we dead on ?
               y2(:,i,:) = y1(:,j,:)
               EXIT
            ELSEIF( j == n1) THEN ! Were at the end and still havent found a bounds.... just give it last point
               y2(:,i,:) = y1(:,n1,:)
               EXIT
            ELSEIF( x1(j) < x2(i) ) THEN ! Are we too small
               j = j + 1
               jm= j - 1
            ELSEIF( x1(j) > x2(i) .and. x1(jm) < x2(i)  ) THEN ! In range... interpolate
               jm = j - 1
               y2(:,i,:) = y1(:,jm,:) + ( y1(:,j,:) - y1(:,jm,:))/(x1(j)-x1(jm))*(x2(i)-x1(jm))
               EXIT
            ELSE                ! Too far right... move left one
               j = j - 1
               jm= j - 1
            END IF
         END DO
      END DO


END SUBROUTINE



!!  Rescale and recycle the flow at a given plane and use as the inlet boundary condition
!!  See Urbin and Knight 2001
SUBROUTINE inflow_rcy_v2(rho,u,v,w,e) !!,rcy_pt_g,del_BL,Len_BL,U_in)
      USE mpi, ONLY: xyzcom_id,master
      USE constants, ONLY: zero,half,one,two,four,five,six,ten,Runiv
      USE globals, ONLY: ax,ay,az,molwts,ix,iy,iz,simtime,x1proc,dump
      USE globals, ONLY: jobdir,step,dt
      USE nozzle_data
      !USE interfaces, ONLY: curl_c
      USE inputs, ONLY: nx,ny,nz,gamma,dx,dy,dz
      USE interfaces, ONLY : SUBSUM3XZ,newdir
      USE metrics, ONLY: y_c

      IMPLICIT NONE
      DOUBLE PRECISION, DIMENSION(ax,ay,az), INTENT(INOUT) :: rho,u,v,w,e
      DOUBLE PRECISION, DIMENSION(1,ay,1) :: tmp
      DOUBLE PRECISION, DIMENSION(ax,ay,az) :: T   ! or maybe use globals
      DOUBLE PRECISION, DIMENSION(nvar,yproc*ay,az) :: Riave,Roave,Riflc,Roflc
      DOUBLE PRECISION, DIMENSION(nvar,yproc*ay,az) :: Qflow,Qflc
      DOUBLE PRECISION, DIMENSION(nvar,ay,az) :: Lplane
      DOUBLE PRECISION, DIMENSION(yproc*ay,az) :: wL,wwL,tmpP 
      DOUBLE PRECISION, DIMENSION(ay,az) :: wL2 
      DOUBLE PRECISION, DIMENSION(nvar*2,yproc*ay,az) :: InTmp1,InTmp2 

      DOUBLE PRECISION, DIMENSION(ny) :: y_G 
      DOUBLE PRECISION, DIMENSION(yproc*ay) :: y_r

      INTEGER :: side,top,bottom
      INTEGER :: plane,recv,donor
      INTEGER :: yprocID,zprocID
      INTEGER :: inlet,rcy_pt_l
      INTEGER :: uvar,vvar,wvar,Tvar
      INTEGER :: iy1,iyn

      DOUBLE PRECISION :: in_plane,line
      INTEGER :: i,j,k,k_loc,isync,rrunit

      DOUBLE PRECISION :: alpha,beta
      DOUBLE PRECISION :: Talpha,Tbeta
      DOUBLE PRECISION :: t_nm1,t_n,dt_n,Rgas,gm1
      DOUBLE PRECISION :: Uvd_inf,AA,iBB,BB,U_inf
      CHARACTER(len=80) :: uaveFile,uaveDir
      LOGICAL :: onMAP
      

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!         Section 1- Variable Definitions            !!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !! BLplane refers to the (ay*yproc,az) plane on the inlet procs
      !! Pplane refers to the (ay,az) plane which contains some BL data
      !! Qflow - the BLplane containing (u,v,w,T) from the donor station
      !! Qave  - the BLplane of time averaged Qflow data
      !! Qflc  - the BLplane of fluctuating Qflow data
      !! Lplane- the Pplane containing (u,v,w,T) from the donor station
      !! Roave - the interpolated outer layer, averaged variables
      !! Riave - the interpolated inner layer, averaged variables
      !! Roflc - the interpolated outer layer, fluctuating variables
      !! Riflc - the interpolated inner layer, fluctuating variables

      !! wL - weight function to apply lower RR BL
      !! wwL - weight function to apply inner/outer layer for lower BL

      !! y_G - raw global 1d array for y location
      !! yp_r- rr station global wall distance
      !! yp_i- inlet station inner normalized wall distance
      !! yp_o- inlet station outer normalized wall distance

      !! InTmp1- Temp array for interpZv input
      !! InTmp2- Temp array for interpZv output


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!          Section 2- Constants and Setup            !!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !! Get some local problem variables based on input parameters
      Rgas = Runiv/molwts(1)

      U_inf = U_in
      gm1 = gamma - one
      T = e*gm1/Rgas
      AA = dsqrt( (gm1/two*Mach**2*Pr )/(one+gm1/two*Mach**2*Pr)  )
      Uvd_inf = U_inf/AA * dasin ( AA )
            
      !! Variable Map stuff
      uvar = 1  ! U-velocity map
      vvar = 2  ! V-velocity map
      wvar = 3  ! W-velocity map
      Tvar = 4  ! Temperature map

      !! Constants for inner/outer BL and for RR region in BL units
      iBB = 0.2d0
      BB = 3.0d0

      !! Map-comm stuff
      onMAP = .FALSE.
      IF(pmapID(1)>=1) onMAP = .TRUE.

      !! Get the proc's BL ID information
      plane = pmapID(1)
      side = pmapID(2)
      zprocID = pmapID(3)
      yprocID = pmapID(4)

      !! Map constants to avoid confusion
      recv   = 1    ! Int Tag for procs on the receiving plane (inflow)
      donor  = 2    ! Int Tag for procs on the donor plane 
      bottom = 1    ! Int Tag for procs on the bottom BL
      top    = 2    ! Int Tag for procs on the top BL


      ! Get recycled plane indices, global(given) and local
      rcy_pt_l = MOD(rcy_pt_g,ax) ! This is the local i index of the recycling plane
      IF (rcy_pt_l .EQ. 0) rcy_pt_l = ax 
      i = rcy_pt_l

      inlet = irecv   ! This is always the case
      
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!        Section 3- Read/Write Temporal Averages     !!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ! Time avergaging coefficients and quantities/fluctuations
      ! Need to write a restart file with averages stored
      t_nm1 = simtime_old 
      t_n = simtime
      simtime_old = simtime
      dt_n = t_n - t_nm1
      Talpha =  t_nm1 / ( t_n )
      Tbeta = one - Talpha
      WRITE (uaveDir,'(2A)') TRIM(jobdir),'/RR'
      

      !! If we've already allocated the time averages, skip, else allocate and initialize
      IF( .NOT. ALLOCATED(Qave) ) THEN
         ALLOCATE(Qave(nvar,ay*yproc,az))
         Qave(uvar,:,:) = U_inf
         Qave(vvar,:,:) = zero
         Qave(wvar,:,:) = zero
         Qave(Tvar,:,:) = T_in
         simtime_old = 0.0d0

         ! Set up directory
         CALL newdir(LEN_TRIM(uaveDir),TRIM(uaveDir),isync)
         
         ! Read File for time averages when NOT at t0=0
         IF (simtime .GT. dt .and. init_stats .ne. .TRUE.) THEN
            IF(xyzcom_id == master) PRINT*,'Reading old RR files'

            ! READ FILE only if on recv. map            
            IF (onMAP .AND. precv) THEN
               WRITE (uaveFile,'(2A,I6.6)') TRIM(uaveDir),'/p',xyzcom_id
               rrunit = 37
               OPEN(UNIT=rrunit,FILE=TRIM(uaveFile),FORM='UNFORMATTED', &
               & STATUS='UNKNOWN')
               READ(rrunit) simtime_old
               READ(rrunit) Qave
               CLOSE(rrunit)
            END IF
         END IF
      END IF
      
    
      ! Check to see if dump has incremented
      ! If it has, write new Uave file before they change for this next time step
      IF (res_dump .NE. dump) THEN
         ! WRITE FILE
         IF (onMAP .AND. precv) THEN
            WRITE (uaveFile,'(2A,I6.6)') TRIM(uaveDir),'/p',xyzcom_id
            rrunit = 17
            OPEN(UNIT=rrunit,FILE=TRIM(uaveFile),FORM='UNFORMATTED', &
            & STATUS='REPLACE')
            WRITE(rrunit) simtime_old
            WRITE(rrunit) Qave
            CLOSE(rrunit)
         END IF
      END IF
      res_dump = dump


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!          Section 4- Setup and send Data            !!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       
      ! Get the planes of data from the donor proc, to the recv procs.... shifting and flipping done here
      Lplane(uvar,:,:) = u(rcy_pt_l,:,:)  ! here i is local proc x index for donor point.  Is wrong on non-donor procs but doesnt matter as data not used
      Lplane(vvar,:,:) = v(rcy_pt_l,:,:)
      Lplane(wvar,:,:) = w(rcy_pt_l,:,:)
      Lplane(Tvar,:,:) = T(rcy_pt_l,:,:)
      
      CALL send_plane(Lplane,Qflow,pdonor,precv,nvar,yproc)

      ! Temporal averages... this is a running average that can be "flushed" with the init_stats flag
      Qave = Talpha * Qave + Tbeta * Qflow

      ! Fluctuating quantities calculated
      Qflc = Qflow - Qave
          
      ! Convert to Van Driest velocity (only done for U velocity)
      tmpP = Qave(uvar,:,:)
      Qave(uvar,:,:) = U_inf/AA * dasin( AA*tmpP / U_inf )
      Uvd_inf = U_inf/AA * dasin ( AA )

      
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!          Section 5- Setup boundary layer vars      !!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   

      ! Get the global line of u_tau and del_nu (for y+)
      ! alpha = (del_rec/del_inl) where del is the BL thickness.. > 1
      ! beta = u_tau_inl / u_tau_rec where u_tau is the friction velocity
      ! Equation from Smits and Dussauge for TBL growth
      IF (BLalpha == 0) THEN   !! Use input alpha unless its zero, then use this relation
         alpha = (( one + (Len_BL)/del_BL * .27d0**(six/five)*Re_BL**(-one/five) )**(five/six)) 
      ELSE
         alpha = BLalpha
      END IF
      beta =  alpha**(one/ten)


      IF(.NOT. ALLOCATED(yp_r)) THEN

         ALLOCATE(yp_r(yproc*ay))
         ALLOCATE(yp_i(yproc*ay))
         ALLOCATE(yp_o(yproc*ay))

         in_plane = zero
         line = zero
         IF (ix(1) == 1) in_plane = one
         IF (iz(1) == 1) line = one
         tmp(1,:,1) = y_c(1,:,1)*in_plane*line
         y_G = SUBSUM3XZ(tmp)
         y_r = y_G(1:ay*yproc)  ! Only take the needed height out... not whole profile... also, assume symmtry of BLs...
         

         ! bottom boundary layer interpolate
         yp_r = y_r - y_r(1) 
         yp_r = yp_r / del_BL ! We non-dimensionalize the height here by del_BL, yp_in = "/eta"
         yp_i = yp_r * beta
         yp_o = yp_r * alpha
      
      END IF

      ! Lower wall weights
      DO k=1,az
         wwL(:,k) = half*( one + ( dtanh( (four*(yp_r-iBB))/((one-two*iBB)*yp_r + iBB))/dtanh(four)))
         wL(:,k) = half*(one - dtanh( (yp_r-BB)/(.1d0)))
      END DO
      

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!          Section 6- Interpolation & Rescaling      !!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      IF (precv .AND. onMAP) THEN  ! Only the actual RR procs do this (garbage data will break interpZv)
      
         InTmp1(1:nvar,:,:)        = Qave
         InTmp1(nvar+1:2*nvar,:,:) = Qflc
         
         ! Inner layer mean & flc
         CALL interpZv(yp_r,yp_i,InTmp1,InTmp2,ay*yproc,ay*yproc,az,nvar*2)
         Riave = InTmp2(1:nvar,:,:)
         Riflc = InTmp2(nvar+1:2*nvar,:,:)
         
         ! Outer layer mean & flc
         CALL interpZv(yp_r,yp_o,InTmp1,InTmp2,ay*yproc,ay*yproc,az,nvar*2)
         Roave = InTmp2(1:nvar,:,:)
         Roflc = InTmp2(nvar+1:2*nvar,:,:)

      END IF

      !! Scaling on flow variables here

      ! Scale the Fluctuations
      Roflc = Roflc * beta  ! Outer layer
      Riflc = Riflc * beta  ! Inner layer

      ! Scale Mean stream-wise velocity, inner layer
      tmpP = Riave(uvar,:,:) * beta
      Riave(uvar,:,:) = U_inf/AA *  dsin( AA*tmpP / U_inf ) ! Transform back from Van Driest

      ! Scale Mean stream-wise velocity, outer layer
      tmpP = Roave(uvar,:,:)* beta + (one-beta)*Uvd_inf
      Roave(uvar,:,:) = U_inf/AA * dsin( tmpP/U_inf * AA) ! Transform back from Van Driest

      ! Convert Mean back from Van Driest velocity (only done for U velocity)
      tmpP = Qave(uvar,:,:)
      Qave(uvar,:,:) = U_inf/AA * dsin( AA*tmpP / U_inf )


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!          Section 7- Reassemble and Reassign        !!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ! Form the weighted sum (inner and outer regions)
      DO i=1,nvar
         Qflow(i,:,:) = (Riave(i,:,:) + Riflc(i,:,:))*(one - wwL) + (Roave(i,:,:) + Roflc(i,:,:))*wwL
      END DO

      ! Get the data in the right place for the inlet procs
      IF (precv .AND. onMAP ) THEN
         
         IF(side==bottom) THEN     !! Bottom
            iy1 = (yprocID-1)*ay + 1
            iyn = iy1 + ay - 1
         END IF
         
         IF(side==top) THEN !! Top
            iy1 = yproc*ay - yprocID*ay + 1
            iyn = iy1 + ay - 1
            
            !! Now, re-order top BL data to have BL start at j=ay*yproc, not j=1
            Qflow = Qflow(:,ay*yproc:1:-1,:)
            wL = wL(yproc*ay:1:-1,:)
         END IF

         ! Blend the calc. RR inflow only in pertainant region and add to inlet value
         wL2 = wL(iy1:iyn,:)
         u(inlet,:,:) = Qflow(uvar,iy1:iyn,:)*wL2 + u(inlet,:,:)*(one-wL2)
         v(inlet,:,:) = Qflow(vvar,iy1:iyn,:)*wL2 + v(inlet,:,:)*(one-wL2)
         w(inlet,:,:) = Qflow(wvar,iy1:iyn,:)*wL2 + w(inlet,:,:)*(one-wL2)
         T(inlet,:,:) = Qflow(Tvar,iy1:iyn,:)*wL2 + T(inlet,:,:)*(one-wL2)

         e(inlet,:,:) = T(inlet,:,:)*Rgas/gm1
         rho(inlet,:,:) = (P_in/gm1)/e(inlet,:,:)
         
         ! Ensure 2d... otherwise things drift
         IF(nz==1) w = zero

      END IF


       
END SUBROUTINE



SUBROUTINE random_BL(u,v,w)
      USE mpi
      USE inputs, ONLY: nx,ny,nz,gfilter,lfilter
      USE globals, ONLY: ax,ay,az,ix,iy,iz,y1proc,ynproc,x1proc,xnproc
      USE metrics, ONLY : x_c,y_c
      USE constants, ONLY: zero,half,one,two
      USE interfaces, ONLY: restart,filter,filtery,set_time_seed,ran1,SUBSUM3YZ
      USE nozzle_data
      IMPLICIT NONE
      DOUBLE PRECISION, DIMENSION(ax,ay,az),INTENT(INOUT) :: u,v,w
      DOUBLE PRECISION, DIMENSION(ax,ay,az) :: wU,wL,wN,ywall,weight
      DOUBLE PRECISION, DIMENSION(ax,1,1) :: tmp
      DOUBLE PRECISION, DIMENSION(ay,az) :: up,vp,wp
      DOUBLE PRECISION, DIMENSION(ny) :: y_U,y_L
      DOUBLE PRECISION, DIMENSION(size(u,1),size(u,2),size(u,3))  :: tmpy,tmpx,tmpu
      INTEGER :: in_unit,i,j,k,re_unit,pt
      DOUBLE PRECISION :: dumy,bL,dum1,dum2,shift,umag,BB,f1w,f2w,Uinf,dt_step,thick
      DOUBLE PRECISION :: plane
      INTEGER :: nxp,nyp,wall
      CHARACTER(LEN=30) :: comments,filename

      !! Use this to get smooth profile at wall
      IF(y1proc) u(:,1,:)  = zero
      IF(ynproc) u(:,ay,:) = zero
      IF(ynproc .and. x1proc) print*,'yeppers'
      DO i=1,10
         CALL filtery(gfilter,u,tmpx)
         CALL filtery(gfilter,tmpx,u)
      END DO

      CALL set_time_seed(1)
      CALL MPI_BARRIER(MPI_COMM_WORLD,mpierr)
      
      DO i=1,ax
         DO j=1,ay
            DO k=1,az
               wL(i,j,k) = half*( one - tanh( dble(iy(j)-30) / 4))
               wU(i,j,k) = half*( one + tanh( dble(iy(j)-( ny-30) ) / 4))
            END DO
         END DO
      END DO

      weight = wL + wU

      weight = weight !* .1d0

      f1w = .3d-1    ! High frquency
      f2w = .4d0    ! Low  frequency
      
      tmpu = u
      
      CALL ran1(tmpy)
      tmpy = two*(tmpy-half)           ! Make RN from (-1,1)
      CALL filter(lfilter,tmpy,tmpx)   ! High Wave Number filter
      u = u * (one + weight * f1w*tmpx) ! /sqrt(tmpy))
      CALL ran1(tmpy)
      tmpy = two*(tmpy-half)           ! Make RN from (-1,1)
      CALL filter(gfilter,tmpy,tmpx)   ! Low Pass filter
      u = u * (one + weight * f2w*tmpx) ! /sqrt(tmpy))
      
      
      
      f1w = .1d-1                ! High frquency
      f2w = .2d0                ! Low  frequency

      CALL ran1(tmpy)
      tmpy = two*(tmpy-half)           ! Make RN from (-1,1)
      CALL filter(lfilter,tmpy,tmpx)   ! High Wave Number filter
      v = v +  tmpu* (weight * f1w*tmpx) ! /sqrt(tmpy))
      CALL ran1(tmpy)
      tmpy = two*(tmpy-half)           ! Make RN from (-1,1)
      CALL filter(gfilter,tmpy,tmpx)   ! Low Pass filter
      v = v + tmpu* (weight * f2w*tmpx) ! /sqrt(tmpy))


      f1w = .1d-1                ! High frquency
      f2w = .2d0                ! Low  frequency

      CALL ran1(tmpy)
      tmpy = two*(tmpy-half)           ! Make RN from (-1,1)
      CALL filter(lfilter,tmpy,tmpx)   ! High Wave Number filter
      w = w + tmpu* ( weight * f1w*tmpx) ! /sqrt(tmpy))
      CALL ran1(tmpy)
      tmpy = two*(tmpy-half)           ! Make RN from (-1,1)
      CALL filter(gfilter,tmpy,tmpx)   ! Low Pass filter
      w = w +  tmpu* ( weight * f2w*tmpx) ! /sqrt(tmpy))

      IF(y1proc) THEN
         u(:,1,:) = zero
         v(:,1,:) = zero
         w(:,1,:) = zero
      END IF

      IF(ynproc) THEN
         u(:,ay,:) = zero
         v(:,ay,:) = zero
         w(:,ay,:) = zero
      END IF

      
END SUBROUTINE random_BL

!! This routine will send a plane of yz data that is local on a domain-yz plane
!! to another domain-yz plane.  The procs on the recv. plane will have the entire 
!! plane of data, not just the proc data.  Useful for shifting data in the z-direction
!! (spanwise) and for interpolation in the y-direction (wall normal).
SUBROUTINE send_plane(Lplane,Gplane,pdonor,precv,nvar,yproc)
USE inputs, ONLY: nx,ny,nz,px,py,pz
USE globals, ONLY: ax,ay,az
USE mpi
USE nozzle_data, ONLY: MAP,pmapID,rcy_pt_g,flippy
IMPLICIT NONE

INTEGER, INTENT(IN) :: nvar, yproc
LOGICAL, INTENT(IN) :: pdonor,precv
DOUBLE PRECISION, DIMENSION(nvar,ay,az), INTENT(IN)  :: Lplane
DOUBLE PRECISION, DIMENSION(nvar,ay*yproc,az), INTENT(OUT) :: Gplane
DOUBLE PRECISION, DIMENSION(nvar,ay,az) :: Lbuff
DOUBLE PRECISION, DIMENSION(nvar,ay*yproc,az) :: Gbuff
INTEGER :: sproc,rproc
INTEGER :: iy1,iyn,iz1,izn
INTEGER :: ii,cycle,plane,side,zprocID,yprocID,recv,donor,zproc,top,bot
LOGICAL :: onMAP
INTEGER :: req,rstat(MPI_STATUS_SIZE),istart




!! Initialize the plane to zero
Gplane = 0.0D0
Lbuff = Lplane

recv = 1
donor = 2
top = 2
bot = 1
 

!! Are we on the map?
onMAP = .FALSE.
IF(pmapID(1)>=1) onMAP = .TRUE.

IF(onMAP) THEN
      
      !! Get the proc's ID information
      plane = pmapID(1)
      side = pmapID(2)
      zprocID = pmapID(3)
      yprocID = pmapID(4)  !! goes from (1-yproc)

      IF(flippy) side = MOD(side,2)+1      ! flip the sides

      !! Loop over the y-procs in this plane and z-proc
      istart = 1
      DO ii=istart,yproc
         

         IF(pdonor) THEN
         
            cycle = MOD( yprocID+ii-2 , yproc ) + 1 !! goes from (1-yproc)... cycle upwards
            rproc = MAP( recv , side, zprocID, cycle , 1)

            IF (rproc .NE. xyzcom_id) THEN
               !IF(side==top) print*,xyzcom_id,' -> ',rproc
               CALL MPI_iSEND(Lplane,ay*az*nvar,MPI_DOUBLE_PRECISION,rproc,mpitag,xyzcom,req,mpierr)
            ELSE  ! Send nothing if you're your own proc
            END IF
         END IF

         IF(precv) THEN
            
            cycle = MOD( yproc + yprocID - (ii-1) -1 , yproc ) + 1 !! goes from (yproc-1)... cycle downwards
            sproc = MAP( donor, side, zprocID, cycle , 1) !! This is just xyzcom_id IF your a donor
            
            IF (sproc .NE. xyzcom_id) THEN
               !IF(side==top) print*,xyzcom_id,' <- ',sproc
               CALL MPI_RECV(Lbuff,ay*az*nvar,MPI_DOUBLE_PRECISION,sproc,0,xyzcom,rstat,mpierr)
            ELSE
               Lbuff = Lplane   ! Data is local so just reassign array
            END IF

            
            IF(side==bot) THEN  !! Bottom
               iy1 = (cycle-1)*ay + 1
               iyn = iy1 + ay - 1
            END IF
            
            IF(side==top) THEN  !! Top
               iy1 = yproc*ay - cycle*ay + 1
               iyn = iy1 + ay - 1
            END IF

            !! Assign to proper location in the larger data array
            Gplane(:,iy1:iyn,:) = Lbuff

         END IF
      END DO
END IF

!! Now, re-order top BL data to have BL start at j=1, not ay*yproc
!! NOTE: we will have to undo this before the data is moved back to each proc.
IF(precv .AND. onMAP .AND. side==top) THEN
      Gplane = Gplane(:,ay*yproc:1:-1,:)
END IF

END SUBROUTINE



!! This is a helper routine for send_plane to get some proc_ids
!! call once at prob_setup
!! This will return a flag to determine if the proc is a donor or a
!! receiver and also a map for the yz-plane of procs for both planes
!! of processors.  These maps will be used later to communicate 

SUBROUTINE setup_maps(MAP,pmapID,pdonor,idonor,precv,irecv,yproc)
USE inputs, ONLY: px,py,pz
USE globals, ONLY: ax,ay,az,ix,iy,iz
USE interfaces, ONLY: SUMproc
USE mpi
IMPLICIT NONE
INTEGER, INTENT(IN) :: idonor,irecv,yproc
INTEGER, DIMENSION(2,2,pz,yproc,3), INTENT(OUT) :: MAP
INTEGER, DIMENSION(4), INTENT(OUT) :: pmapID
LOGICAL, INTENT(OUT) :: pdonor,precv

INTEGER :: i,j,k,pgrid
DOUBLE PRECISION :: tmp1D,tmp2D,tmp3D
INTEGER :: recv,donor,bottom,top,yy1,zz1

!! Map constants to avoid confusion
recv = 1
donor = 2
bottom = 1
top = 2


!! Check the bounds to see if this proc is donor or receipient plane (based on x index number)
pdonor = .FALSE.
IF( idonor >= ix(1) .AND. idonor <= ix(ax) ) pdonor = .TRUE.

precv = .FALSE.
IF( irecv >= ix(1) .AND. irecv <= ix(ax) ) precv = .TRUE.


!! Initialize the map
MAP = 0
pmapID = 0

DO k=1,pz
      !! Count (1-nz) in az units
      zz1 = (k-1)*az + 1
      DO j=1,yproc

         !! Receiving plane-bottom
         
         !! Count (1-ny) in ay units
         yy1 = (j-1)*ay + 1
         tmp1D = 0.0;tmp2D = 0.0;tmp3D = 0.0;
         IF(precv .AND. yy1 == iy(1) .AND. zz1 == iz(1)) THEN
            tmp1D = xyzcom_id
            tmp2D = yy1
            tmp3D = zz1
            pmapID = (/ recv, bottom, k, j /) 
         END IF
         MAP(recv,bottom,k,j,1) = SUMproc(tmp1D)
         MAP(recv,bottom,k,j,2) = SUMproc(tmp2D)
         MAP(recv,bottom,k,j,3) = SUMproc(tmp3D)
      
         !! Donor plane-bottom
         
         !! Count (1-ny) in ay units
         yy1 = (j-1)*ay + 1
         tmp1D = 0.0;tmp2D = 0.0;tmp3D = 0.0;
         IF(pdonor .AND. yy1 == iy(1) .AND. zz1 == iz(1)) THEN
            tmp1D = xyzcom_id
            tmp2D = yy1
            tmp3D = zz1
            pmapID = (/ donor, bottom, k, j /) 
         END IF
         MAP(donor,bottom,k,j,1) = SUMproc(tmp1D)
         MAP(donor,bottom,k,j,2) = SUMproc(tmp2D)
         MAP(donor,bottom,k,j,3) = SUMproc(tmp3D)         


         !! Receiving plane-top
         
         !! Count (ny-1) in ay units
         yy1 = (py - j) * ay + 1
         tmp1D = 0.0;tmp2D = 0.0;tmp3D = 0.0;
         IF(precv .AND. yy1 == iy(1) .AND. zz1 == iz(1)) THEN
            tmp1D = xyzcom_id
            tmp2D = yy1
            tmp3D = zz1
            pmapID = (/ recv, top, k, j /) 
         END IF
         MAP(recv,top,k,j,1) = SUMproc(tmp1D)
         MAP(recv,top,k,j,2) = SUMproc(tmp2D)
         MAP(recv,top,k,j,3) = SUMproc(tmp3D)
      
         !! Donor plane-top
         
         !! Count (ny-1) in ay units
         yy1 = (py - j) * ay + 1
         tmp1D = 0.0;tmp2D = 0.0;tmp3D = 0.0;
         IF(pdonor .AND. yy1 == iy(1) .AND. zz1 == iz(1)) THEN
            tmp1D = xyzcom_id
            tmp2D = yy1
            tmp3D = zz1
            pmapID = (/ donor, top, k, j /) 
         END IF
         MAP(donor,top,k,j,1) = SUMproc(tmp1D)
         MAP(donor,top,k,j,2) = SUMproc(tmp2D)
         MAP(donor,top,k,j,3) = SUMproc(tmp3D)


      END DO
END DO
         
         


END SUBROUTINE


!! Figure out the number of procs that contain the Boundary Layer
SUBROUTINE yprocBL()
USE mpi
USE globals
USE inputs
USE nozzle_data
USE metrics
USE interfaces, ONLY: SUMproc
IMPLICIT NONE
DOUBLE PRECISION :: tmpD,sumY,ywall


!! Get bottom left-corner y-value
tmpD = 0.0
IF (x1proc .AND. iz(1) == 1 .AND. y1proc) THEN
      tmpD = y_c(1,1,1)
END IF
sumY = SUMproc(tmpD)
ywall = sumY

!! Check to see if you are on top or bottom of grid
tmpD = 0.0
IF (x1proc .AND. iz(1) == 1) THEN
      IF( (y_c(1,1,1)-ywall) < del_BL*3.0d0) THEN
         tmpD = 1.0
      END IF
END IF

!! Sum on all the procs that contain BL and send to all procs
sumY = SUMproc(tmpD)

!! Cast as an INTEGER
yproc = INT(sumY)

END SUBROUTINE
