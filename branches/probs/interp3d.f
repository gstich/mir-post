! you must re-define this module name to be unique across all problems
MODULE interp3d_data
  USE mpi
  SAVE
  ! ---------------------------------------------------------------------------
  ! Place Problem dependent state data here, similar to inputs.f and globals.f
  ! ---------------------------------------------------------------------------  
  CHARACTER(LEN=90) 	:: FROM_dir  = '/path/to/directory'
  CHARACTER(LEN=90) 	:: TO_dir  = '/path/to/directory'

  CHARACTER(LEN=90) 	:: TO_grid     = 'nozzle_TO.grid'
  CHARACTER(LEN=90) 	:: FROM_grid   = 'nozzle_FROM.grid'
  DOUBLE PRECISION      :: FROM_dz  = 7.0D-4
  DOUBLE PRECISION      :: TO_dz    = 7.0D-4

   
  INTEGER               :: FROM_num  = 200  
  INTEGER               :: TO_num    = 200  

  INTEGER               :: varDIM = 8                  ! Variable dimension of the incoming data array (12 is for compressible, single fluid)
  CHARACTER(LEN=3)      :: ftype = 'res'                ! Select file type: 'res'-restart, 'vis'-viz file
  INTEGER               :: buff = 4                ! Number of points on the [F] grid to buffer for local interp

  CHARACTER(LEN=90) 	:: MIRname = 'post'
  CHARACTER(LEN=90) 	:: Gname = 'grid'

  ! Not set in Namelist, these variables are read in and used to set the inflow BC
  CHARACTER(LEN=30) 	:: prob_jobname = 'interp3d'


  ! Processor Maps for the the [F]rom data, [T]o data and the [L]ocal data.
  INTEGER, DIMENSION(:,:),   ALLOCATABLE :: procmapT,procmapF,procmapL
  INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: invprocmapT,invprocmapF,invprocmapL
  DOUBLE PRECISION, DIMENSION(:,:,:,:), ALLOCATABLE :: Tiodata,Fiodata,TMPiodata,GRDdata

  ! Grid info
  DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE :: Txc,Tyc,Tzc,Fxc,Fyc,Fzc  

  INTEGER               :: Tax, Tay, Taz                ! Grid point per-proc for dataset
  INTEGER               :: Fax, Fay, Faz                ! Grid point per-proc for dataset
  INTEGER               :: Tnx, Tny, Tnz                ! Total grid point
  INTEGER               :: Fnx, Fny, Fnz                ! Total grid points
  INTEGER               :: Tpx = 1                      ! Processors in x-direction for dataset
  INTEGER               :: Tpy = 1                      ! Processors in y-direction for dataset
  INTEGER               :: Tpz = 1                      ! Processors in z-direction for dataset

  INTEGER               :: Fpx = 1                      ! Processors in x-direction for dataset
  INTEGER               :: Fpy = 1                      ! Processors in y-direction for dataset
  INTEGER               :: Fpz = 1                      ! Processors in z-direction for dataset


END MODULE interp3d_data

! -----------------------------------------------------------------------------
! prob_inputs
! Read problem-dependent inputs and save in module variables
! must also define 'jobname'
! Called by the Miranda main program.
! -----------------------------------------------------------------------------
SUBROUTINE prob_inputs(fileName)
  USE mpi
  USE prob_interface, ONLY : jobname
  USE interp3d_data
  USE globals, ONLY : inputUnit,molwts
  USE constants, ONLY: Runiv,one
  USE inputs, ONLY : nx,ny,nz,gamma
  USE interfaces, ONLY : SUMproc
  IMPLICIT NONE
  CHARACTER(LEN=*), INTENT(IN) :: fileName
  INTEGER :: dsI,bsx,bsy,bsz
  INTEGER :: in_unit,i 
  CHARACTER(LEN=90) :: plotmir,comments


  ! Uncomment to define a namelist
  NAMELIST /interp3d_vars/ FROM_dir, TO_dir, FROM_grid, TO_grid, ftype, FROM_num, TO_num, & 
      TO_dz, FROM_dz



  ! Uncomment to open and read a namelist file
  OPEN(UNIT=inputUnit,FILE=TRIM(fileName),FORM='FORMATTED',STATUS='OLD')
  READ(UNIT=inputUnit,NML=interp3d_vars)
  CLOSE(inputUnit)
  
  ! give the right name
  jobname = TRIM(prob_jobname)
  
  
  SELECT CASE(ftype)
      CASE('res')
         varDIM = 8     ! Restart file.  8 doubles/grid point
      CASE('vis')
         varDIM = 12    ! Vis file. 12 real*4/grid point
  END SELECT


  ! Read in the procmap file and set Tpx,Tpy,Tpz
  dsI = 5
  nx = 0
  ny = 0
  nz = 0


  ! Read the TO data file  (interp3d onto this grid!)
  WRITE(plotmir,'(2A)') TRIM(TO_dir),'/plot.mir'
  OPEN(UNIT=inputUnit,FILE=TRIM(plotmir),FORM='FORMATTED',STATUS='OLD')
  DO i=1,dsI
      READ(inputUnit,*) comments
  END DO
  READ(inputUnit,*) comments, Tnx, Tny, Tnz
  READ(inputUnit,*) comments, Tax, Tay, Taz
  CLOSE(inputUnit)  

  Tpx = Tnx / Tax
  Tpy = Tny / Tay
  Tpz = Tnz / Taz  

  ! Read the FROM data file (interp3d using this data/grid)
  WRITE(plotmir,'(2A)') TRIM(FROM_dir),'/plot.mir'
  OPEN(UNIT=inputUnit,FILE=TRIM(plotmir),FORM='FORMATTED',STATUS='OLD')
  DO i=1,dsI
      READ(inputUnit,*) comments
  END DO
  READ(inputUnit,*) comments, Fnx, Fny, Fnz
  READ(inputUnit,*) comments, Fax, Fay, Faz
  CLOSE(inputUnit)  

  Fpx = Fnx / Fax
  Fpy = Fny / Fay
  Fpz = Fnz / Faz 


  ! Set the local nx,ny,nz to the [F]rom size  (used for ghosting and such
  nx = Fnx
  ny = Fny
  nz = Fnz



END SUBROUTINE prob_inputs



! -------------------------------------------------------------------------------------
! Set up problem geometry. (This is called on startup and restart before timestepping.)
! -------------------------------------------------------------------------------------
 SUBROUTINE prob_setup()
  USE mpi
  USE constants
  USE inputs, ONLY: px,py,pz
  !USE globals, ONLY:
  !USE extend
  USE interp3d_data
  USE interfaces, ONLY: newdir, restart, viz
  IMPLICIT NONE
  INTEGER :: i,j,funit,isync
  CHARACTER(LEN=90) :: flipfile,fform,plotdir,vizdir,grdir
  LOGICAL :: fexist
  DOUBLE PRECISION :: Rgas,U_in,e_in,T_in
  INTEGER :: TFstart,TFcount,TFend,Nstride,Nsten,ii,TS_old
  DOUBLE PRECISION :: sig,Gnorm
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: Gauss
  DOUBLE PRECISION, DIMENSION(:,:,:,:), ALLOCATABLE :: fil_data
  DOUBLE PRECISION :: Tin
  INTEGER :: nx,ny,nz,ax,ay,az

  ! Set-up the grid and metric terms for the calculation
  !CALL prob_geom()


  ! Get the procmaps and invprocmaps
  CALL setup_dataread()  

  

  !!!!!!!!!! Read in the [F]rom data   !!!!!!!!!!!!
  nx = Fax*Fpx
  ny = Fay*Fpy
  nz = Faz*Fpz
  ax = nx/px
  ay = ny/py
  az = nz/pz
  ALLOCATE(Fxc(ax,ay,az))
  ALLOCATE(Fyc(ax,ay,az))
  ALLOCATE(Fzc(ax,ay,az))
  ALLOCATE(Fiodata(ax,ay,az,varDIM))
  
  WRITE(vizdir,'(3A,I4.4)') TRIM(FROM_dir),'/',TRIM(ftype),FROM_num  
  WRITE(grdir,'(2A)') TRIM(FROM_dir),'/grid'
  ALLOCATE(TMPiodata(ax,ay,az,varDIM))
  ALLOCATE(GRDdata(ax,ay,az,3))
  CALL read_data(vizdir,Fax,Fay,Faz,Fpx,Fpy,Fpz,invprocmapF,ftype)
  CALL read_data(grdir,Fax,Fay,Faz,Fpx,Fpy,Fpz,invprocmapF,'grd')
  Fiodata = TMPiodata
  Fxc = GRDdata(:,:,:,1)
  Fyc = GRDdata(:,:,:,2)
  Fzc = GRDdata(:,:,:,3)
  DEALLOCATE(TMPiodata,GRDdata)




  !!!!!!!!! Read in the [T]o data  !!!!!!!!!!!
  nx = Tax*Tpx
  ny = Tay*Tpy
  nz = Taz*Tpz
  ax = nx/px
  ay = ny/py
  az = nz/pz
  ALLOCATE(Txc(ax,ay,az))
  ALLOCATE(Tyc(ax,ay,az))
  ALLOCATE(Tzc(ax,ay,az))
  ALLOCATE(Tiodata(ax,ay,az,varDIM))
  WRITE(vizdir,'(3A,I4.4)') TRIM(TO_dir),'/',TRIM(ftype),TO_num 
  WRITE(grdir,'(2A)') TRIM(TO_dir),'/grid' 

  !! Read the Viz data on Tdata
  ALLOCATE(TMPiodata(ax,ay,az,varDIM))
  !CALL read_data(vizdir,Tax,Tay,Taz,Tpx,Tpy,Tpz,invprocmapT,ftype)
  !Tiodata = TMPiodata
  DEALLOCATE(TMPiodata)

  !! Read the grid data on Tdata
  CALL prob_grid(Txc,Tyc,Tzc,nx,ny,nz,ax,ay,az,TO_DZ,TO_grid)
  ALLOCATE(GRDdata(ax,ay,az,3))  
  CALL read_data(grdir,Tax,Tay,Taz,Tpx,Tpy,Tpz,invprocmapT,'grd')
  Txc = GRDdata(:,:,:,1)
  Tyc = GRDdata(:,:,:,2)
  Tzc = GRDdata(:,:,:,3)
  DEALLOCATE(GRDdata)
  



  ! Transfer the solution and grid
  !Fiodata = Tiodata
  !Fxc = Txc
  !Fyc = Tyc
  !Fzc = Tzc  
  
  CALL trilinear(Tiodata,Txc,Tyc,Tzc,ax,ay,az,Fiodata,Fxc,Fyc,Fzc,Fax*Fpx/px,Fay*Fpy/py,Faz*Fpz/pz)


  ! Write data to [T]o dir
  ALLOCATE(GRDdata(ax,ay,az,3))  
  ALLOCATE(TMPiodata(ax,ay,az,varDIM))

  TMPiodata = Tiodata  
  GRDdata(:,:,:,1) = Txc
  GRDdata(:,:,:,2) = Tyc
  GRDdata(:,:,:,3) = Tzc

  IF (xyzcom_id == 0) PRINT*,'Here'
  CALL write_data(vizdir,Tax,Tay,Taz,Tpx,Tpy,Tpz,invprocmapT,ftype)
  CALL write_data(grdir,Tax,Tay,Taz,Tpx,Tpy,Tpz,invprocmapT,'grd')
  
  DEALLOCATE(TMPiodata,GRDdata)




  ! Make the program stop here.
  IF (xyzcom_id == 0) PRINT*,'Post-proc has finished... exiting Miranda'
  CALL EXIT(0)


END SUBROUTINE prob_setup



! -----------------------------------------------------------------------------
! prob_init
! Initialize state variables at start of run
! Called by the Miranda main program.
! -----------------------------------------------------------------------------
SUBROUTINE prob_init(rho,u,v,w,e,Y,p,T)



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
  USE interp3d_data
  USE metrics, ONLY: x_c,y_c,z_c
  IMPLICIT NONE
  CHARACTER(LEN=*), INTENT(IN) :: plotdir
  CHARACTER(LEN=flen) :: plotFile,var
  INTEGER, PARAMETER :: plotUnit=21
  INTEGER :: i,k
  DOUBLE PRECISION :: x_a,x_b,y_a,y_b,z_a,z_b,core,x0
  DOUBLE PRECISION, DIMENSION(ny) :: ubar,vbar,wbar,rhobar,pbar,ybar
  DOUBLE PRECISION, DIMENSION(nx) :: p_dwn,p_up,p_cen,xbar
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: corr
  CHARACTER(LEN=flen) :: myformat
  DOUBLE PRECISION, DIMENSION(nz) :: corrT
  INTEGER :: nxs,nys
  DOUBLE PRECISION, DIMENSION(ax,ay,az) :: rhop,up,vp,wp,pp
  DOUBLE PRECISION, DIMENSION(ax,ay,az) :: tmp0,tmp1,box
  



END SUBROUTINE prob_plots

SUBROUTINE prob_source()
END SUBROUTINE




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
  USE interp3d_data
  USE metrics
  USE interfaces, ONLY : filter !, ran1, SUMproc
  IMPLICIT NONE
  DOUBLE PRECISION, DIMENSION(:,:,:), INTENT(INOUT) :: rho,u,v,w
  DOUBLE PRECISION, DIMENSION(:,:,:), INTENT(INOUT), OPTIONAL :: e
  DOUBLE PRECISION, DIMENSION(:,:,:,:), INTENT(INOUT), OPTIONAL :: Y
  DOUBLE PRECISION, DIMENSION(ax,ay,az) :: dumF,dumT
  DOUBLE PRECISION, DIMENSION(size(rho,1),size(rho,3))  :: tmpw,blowx,blowz,blow
  DOUBLE PRECISION, DIMENSION(SIZE(u,1),SIZE(u,3)) :: dxdA,dydA,dzdA,t_mag,u_mag	!  A-column of Jacobian Tensor
  !DOUBLE PRECISION, DIMENSION(1) :: randT
  DOUBLE PRECISION :: del
  DOUBLE PRECISION :: filpt,thick
  INTEGER :: i,k


  
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
  IMPLICIT NONE


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
 USE interp3d_data
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
  OPEN(UNIT=re_unit,FILE=TO_grid,STATUS='OLD')
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

  ! Make sure Jacobian is 2d... extruded mesh
  nz_tmp = nz
  nz = 1
  CALL get_jacobian()
  nz = nz_tmp

END SUBROUTINE

!-------------------------------------------------------------------------------
! Assign problem-specific geometry. This routine is only called if
! upon initialization only
! ------------------------------------------------------------------------------ 
SUBROUTINE prob_grid(XL,YL,ZL,Tnx,Tny,Tnz,Lax,Lay,Laz,LDZ,gridfile)
 USE constants
 USE metrics, ONLY: x_c,y_c,z_c
 USE mpi
 USE constants
 USE interfaces, ONLY : filter
 IMPLICIT NONE
 INTEGER, INTENT(IN) :: Tnx,Tny,Tnz,Lax,Lay,Laz,LDZ
 CHARACTER(LEN=90), INTENT(IN) :: gridfile
 DOUBLE PRECISION, DIMENSION(Lax,Lay,Laz), INTENT(OUT) :: XL,YL,ZL

 DOUBLE PRECISION, DIMENSION(Tnx,Tny) :: tmpx,tmpy
 INTEGER, DIMENSION(Lax) :: ix
 INTEGER, DIMENSION(Lay) :: iy
 INTEGER, DIMENSION(Laz) :: iz

 INTEGER :: i,j,k,nz_tmp
 INTEGER :: re_unit


  re_unit = 19
  
  
  !  Put file read routine here 
  re_unit = 32
  OPEN(UNIT=re_unit,FILE=gridfile,STATUS='OLD')
  DO J=1,Tny
      DO I = 1,Tnx
         READ(re_unit,*) tmpx(I,J),tmpy(I,J)
      END DO
  END DO
  CLOSE(re_unit)

  
  FORALL (i=1:Lax) ix(i) = xcom_id*Lax + i
  FORALL (j=1:Lay) iy(j) = ycom_id*Lay + j
  FORALL (k=1:Laz) iz(k) = zcom_id*Laz + k
   
  XL(:,:,1) = tmpx(ix,iy)
  YL(:,:,1) = tmpy(ix,iy)


  !  Z-direction extrusion
  DO k=1,Laz
      ZL(:,:,k) = zcom_id*Laz*Ldz + Ldz*dble(k-1)
      XL(:,:,k) = XL(:,:,1)
      YL(:,:,k) = YL(:,:,1)
  END DO

  !  For Nozzle-grid convert mm -> cm
  !  Z is left out because it is specified through dz
  XL = XL*1.0d-1
  YL = YL*1.0d-1


  !! Make sure everyone waits till script if finished 
  CALL MPI_BARRIER(MPI_COMM_WORLD,mpierr)

END SUBROUTINE




!! This sets up the variables needed for each processor to read and set the appropriate data.
!! This is based on nx,ny,nz and px,py,pz and Tpx,Tpy,Tpz and this proc's ID
SUBROUTINE setup_dataread()
  USE globals, ONLY: ix,iy,iz,ax,ay,az
  USE inputs, ONLY: nx,ny,nz
  USE interp3d_data
  IMPLICIT NONE

  INTEGER :: i,j,k,p,proc
  DOUBLE PRECISION :: tmp

  !! Set up the proc-map for the data
  proc = Tpx*Tpy*Tpz
  ALLOCATE(procmapT(proc,4))
  ALLOCATE(invprocmapT(Tpx,Tpy,Tpz))
  
  !  Get the proc IJK layout
  p = 1
  DO k=1,Tpz
     DO j=1,Tpy
        DO i=1,Tpx
           procmapT(p,1)=p-1
           procmapT(p,2)=i-1
           procmapT(p,3)=j-1
           procmapT(p,4)=k-1
           invprocmapT(i,j,k) = p-1
           p = p + 1
        END DO
     END DO
  END DO

  !! Set up the proc-map for the data
  proc = Fpx*Fpy*Fpz
  ALLOCATE(procmapF(proc,4))
  ALLOCATE(invprocmapF(Fpx,Fpy,Fpz))
  
  !  Get the proc IJK layout
  p = 1
  DO k=1,Fpz
     DO j=1,Fpy
        DO i=1,Fpx
           procmapF(p,1)=p-1
           procmapF(p,2)=i-1
           procmapF(p,3)=j-1
           procmapF(p,4)=k-1
           invprocmapF(i,j,k) = p-1
           p = p + 1
        END DO
     END DO
  END DO


END SUBROUTINE setup_dataread


!! Each processor gets the data it needs for this particular time step/restart file (should handle both)
!! At the end of this routine iodata, should be fully set.
SUBROUTINE read_data(vizdir,Tax,Tay,Taz,Tpx,Tpy,Tpz,invprocmap,ftype)
  USE mpi
  USE globals, ONLY: iodata,ix,iy,iz,ax,ay,az
  USE constants, ONLY: zero, one
  USE inputs, ONLY: px,py,pz
  USE interp3d_data, ONLY: varDIM,TMPiodata,GRDdata
  IMPLICIT NONE
  CHARACTER(LEN=90), INTENT(IN) :: vizdir
  CHARACTER(LEN=3), INTENT(IN) :: ftype
  INTEGER, INTENT(IN) :: Tax,Tay,Taz,Tpx,Tpy,Tpz
  INTEGER, DIMENSION(Tpx,Tpy,Tpz), INTENT(IN) :: invprocmap
  INTEGER :: i1g,ifg,j1g,jfg,k1g,kfg

  DOUBLE PRECISION :: tmp
  INTEGER :: total,count,i,j,k,p,v
  INTEGER :: ig,jg,kg,iL,jL,kL
  INTEGER :: i1p,ifp,j1p,jfp,k1p,kfp
  INTEGER :: x1off,xnoff,y1off,ynoff,z1off,znoff
  INTEGER :: Lax,Lay,Laz
  CHARACTER(LEN=90)  :: procdir

  INTEGER :: punit=42
  DOUBLE PRECISION, DIMENSION(:,:,:,:), ALLOCATABLE :: Gdata,RESdata
  REAL(KIND=4), DIMENSION(:,:,:,:), ALLOCATABLE :: VISdata,MSHdata


  ! Get the global indice bounds for current proc allocation
  Lax = Tax*Tpx / px   !! Local proc
  Lay = Tay*Tpy / py   
  Laz = Taz*Tpz / pz
  
  ! Global indices for Local proc. layout
  i1g = xcom_id*Lax + 1
  ifg = xcom_id*Lax + Lax

  j1g = ycom_id*Lay + 1
  jfg = ycom_id*Lay + Lay

  k1g = zcom_id*Laz + 1
  kfg = zcom_id*Laz + Laz

  ! Get the needed proc bounds
  tmp = dble(i1g)/dble(Tax)
  i1p = CEILING(tmp) - 1
  tmp = dble(ifg)/dble(Tax)
  ifp = CEILING(tmp) - 1

  tmp = dble(j1g)/dble(Tay)
  j1p = CEILING(tmp) - 1
  tmp = dble(jfg)/dble(Tay)
  jfp = CEILING(tmp) - 1

  tmp = dble(k1g)/dble(Taz)
  k1p = CEILING(tmp) - 1
  tmp = dble(kfg)/dble(Taz)
  kfp = CEILING(tmp) - 1


  ! Get the offsets for this block of data
  x1off =  i1g - (i1p*Tax + 1) + 1
  xnoff =  x1off + Lax -1  !  -(ifp*Tax + 1) + ifg + 1
  y1off =  j1g - (j1p*Tay + 1) + 1
  ynoff =  y1off + Lay - 1 !-(jfp*Tay + 1) + jfg + 1
  z1off =  k1g - (k1p*Taz + 1) + 1
  znoff =  z1off + Laz -1  ! -(kfp*Taz + 1) + kfg + 1

  !IF (xyzcom_id == 0) PRINT*, Tax, Tay, Taz
  !            PRINT*, i1g,i1p,ifp
  !            PRINT*, j1g,j1p,jfp
  !            PRINT*, k1g,k1p,kfp
  !            PRINT*,'X-(',x1off,' ',xnoff,')  '
  !            PRINT*,'Y-(',Y1off,' ',Ynoff,')  '
  !            PRINT*,'Z-(',z1off,' ',znoff,')  '
  !END IF
  

  SELECT CASE(ftype)
      CASE('res')
         ALLOCATE(RESdata(Tax,Tay,Taz,varDIM))
      CASE('vis')
         ALLOCATE(VISdata(Tax,Tay,Taz,varDIM))
      CASE('grd')
         ALLOCATE(MSHdata(Tax,Tay,Taz,3))
  END SELECT

  ALLOCATE(Gdata( (ifp-i1p+1)*Tax, (jfp-j1p+1)*Tay, (kfp-k1p+1)*Taz, varDIM))
  Gdata = zero

  ! Loop only over PROCS in pertainate range
  DO i=i1p,ifp
     DO j=j1p,jfp
        DO k=k1p,kfp
           
           count = count + 1
           

           ! Global indices on corner of proc grid
           ig = i*Tax + 1
           jg = j*Tay + 1
           kg = k*Taz + 1

           ! Local indices on corner of proc grid
           iL = (i-i1p)*Tax + 1
           jL = (j-j1p)*Tay + 1
           kL = (k-k1p)*Taz + 1
           
           !! i,j,k are (0,px0-1).... p is (0,nproc-1)
           p = invprocmap(i+1,j+1,k+1)

           WRITE(procdir,'(2A,I6.6)') TRIM(vizdir),'/p',p
           !PRINT*,procdir

           SELECT CASE(ftype)

           CASE('vis')
           ! Read from visualization data
           OPEN(UNIT=punit,FILE=TRIM(procdir),FORM='UNFORMATTED',STATUS='OLD')
           DO v=1,varDIM
              READ(punit) VISdata(:,:,:,v)
           END DO
           CLOSE(punit)
           ! Slight re-ordering for vis data...
           Gdata(iL:iL+Tax-1,jL:jL+Tay-1,kL:kL+Taz-1,1:5)  = VISdata(:,:,:,1:5)     ! u,v,w,rho,e
           Gdata(iL:iL+Tax-1,jL:jL+Tay-1,kL:kL+Taz-1,6)    = VISdata(:,:,:,12)      ! Y
           Gdata(iL:iL+Tax-1,jL:jL+Tay-1,kL:kL+Taz-1,7:12) = VISdata(:,:,:,6:11)    ! p,T,c,mu,bulk,ktc

           CASE('res')
           ! Read from restart data
           OPEN(UNIT=punit,FILE=TRIM(procdir),FORM='UNFORMATTED',STATUS='OLD',ACTION='READ')
           READ(punit) RESdata
           CLOSE(punit)
           Gdata(iL:iL+Tax-1,jL:jL+Tay-1,kL:kL+Taz-1,:) = RESdata

           CASE('grd')
           ! Read from grid data
           OPEN(UNIT=punit,FILE=TRIM(procdir),FORM='UNFORMATTED',STATUS='OLD',ACTION='READ')
           READ(punit) MSHdata(:,:,:,1)
           READ(punit) MSHdata(:,:,:,2)
           READ(punit) MSHdata(:,:,:,3)
           CLOSE(punit)
           Gdata(iL:iL+Tax-1,jL:jL+Tay-1,kL:kL+Taz-1,1:3) = MSHdata
           END SELECT

        END DO
     END DO
  END DO

  SELECT CASE(ftype)
      CASE('grd')
           GRDdata = Gdata(x1off:xnoff,y1off:ynoff,z1off:znoff,1:3)
      ! Set to the global iodata variable pointer array
      CASE DEFAULT ! For res and vis
         TMPiodata = Gdata(x1off:xnoff,y1off:ynoff,z1off:znoff,:)
  END SELECT

  ! Free the memory
  SELECT CASE(ftype)
      CASE('res')
         DEALLOCATE(Gdata,RESdata)
      CASE('vis')
         DEALLOCATE(Gdata,VISdata)
      CASE('grd')
         DEALLOCATE(Gdata,MSHdata)
  END SELECT


  IF (xyzcom_id == 0) PRINT*, 'Done reading ',TRIM(vizdir)

END SUBROUTINE read_data




!! Each processor gets the data it needs for this particular time step/restart file (should handle both)
!! At the end of this routine iodata, should be fully set.
SUBROUTINE write_data(vizdir,Tax,Tay,Taz,Tpx,Tpy,Tpz,invprocmap,ftype)
  USE mpi
  USE globals, ONLY: iodata,ix,iy,iz,ax,ay,az
  USE inputs, ONLY: px,py,pz
  USE interp3d_data, ONLY: varDIM,TMPiodata,GRDdata
  USE constants, ONLY: zero,one
  USE interfaces, ONLY: newdir
  IMPLICIT NONE
  CHARACTER(LEN=90), INTENT(IN) :: vizdir
  CHARACTER(LEN=3), INTENT(IN) :: ftype
  INTEGER, INTENT(IN) :: Tax,Tay,Taz,Tpx,Tpy,Tpz
  INTEGER, DIMENSION(Tpx,Tpy,Tpz), INTENT(IN) :: invprocmap
  INTEGER :: i1g,ifg,j1g,jfg,k1g,kfg

  DOUBLE PRECISION :: tmp
  INTEGER :: total,count,i,j,k,p,v
  INTEGER :: ig,jg,kg,iL,jL,kL
  INTEGER :: i1p,ifp,j1p,jfp,k1p,kfp
  INTEGER :: x1off,xnoff,y1off,ynoff,z1off,znoff
  INTeGER :: Lax,Lay,Laz
  CHARACTER(LEN=90)  :: procdir

  INTEGER :: punit=42
  DOUBLE PRECISION, DIMENSION(:,:,:,:), ALLOCATABLE :: Gdata,RESdata,Ddata
  REAL(KIND=4), DIMENSION(:,:,:,:), ALLOCATABLE :: VISdata,Rdata,MSHdata
  INTEGER :: Ndata,isync


  ! Get the global indice bounds for current proc allocation
  Lax = Tax*Tpx / px   !! Local proc
  Lay = Tay*Tpy / py   
  Laz = Taz*Tpz / pz

  i1g = xcom_id*Lax + 1
  ifg = xcom_id*Lax + Lax

  j1g = ycom_id*Lay + 1
  jfg = ycom_id*Lay + Lay

  k1g = zcom_id*Laz + 1
  kfg = zcom_id*Laz + Laz

  ! Get the needed proc bounds
  tmp = dble(i1g)/dble(Tax)
  i1p = CEILING(tmp) - 1
  tmp = dble(ifg)/dble(Tax)
  ifp = CEILING(tmp) - 1

  tmp = dble(j1g)/dble(Tay)
  j1p = CEILING(tmp) - 1
  tmp = dble(jfg)/dble(Tay)
  jfp = CEILING(tmp) - 1

  tmp = dble(k1g)/dble(Taz)
  k1p = CEILING(tmp) - 1
  tmp = dble(kfg)/dble(Taz)
  kfp = CEILING(tmp) - 1


  ! Get the offsets for this block of data
  x1off =  i1g - (i1p*Tax + 1) + 1
  xnoff = x1off + Lax - 1 
  y1off =  j1g - (j1p*Tay + 1) + 1
  ynoff = y1off + Lay - 1 
  z1off =  k1g - (k1p*Taz + 1) + 1
  znoff = z1off + Laz - 1 

  SELECT CASE(ftype)
      CASE('res')
         ALLOCATE(RESdata(Tax,Tay,Taz,varDIM))
         ALLOCATE(Ddata(Tax,Tay,Taz,varDIM))
         RESdata = zero
      CASE('vis')
         ALLOCATE(VISdata(Tax,Tay,Taz,varDIM))
         ALLOCATE(Rdata(Tax,Tay,Taz,varDIM))
         VISdata = zero
      CASE('grd')
         ALLOCATE(MSHdata(Tax,Tay,Taz,3))
         ALLOCATE(Rdata(Tax,Tay,Taz,3))
         MSHdata = zero
  END SELECT


  ALLOCATE(Gdata( (ifp-i1p+1)*Tax, (jfp-j1p+1)*Tay, (kfp-k1p+1)*Taz, varDIM))

  CALL newdir(LEN_TRIM(vizdir),TRIM(vizdir),isync)

  ! Loop over ALL procs in the write data set
  DO i=0,Tpx-1
     DO j=0,Tpy-1
        DO k=0,Tpz-1

           ! Reset the Gdata every proc loop
           Gdata = zero
           
           !! i,j,k are (0,px0-1).... p is (0,nproc-1)
           p = invprocmap(i+1,j+1,k+1)
           WRITE(procdir,'(2A,I6.6)') TRIM(vizdir),'/p',p
           IF (xyzcom_id == 0) PRINT*, TRIM(procdir)


           !! Does this new proc, need my local data?  If so, put it in place
           IF ( (i1p<=i .AND. i<=ifp) .AND. (j1p<=j .AND. j<=jfp) .AND. (k1p<=k .AND. k<=kfp)) THEN

              count = count + 1
           
              ! Global indices on corner of proc grid
              ig = i*Tax + 1
              jg = j*Tay + 1
              kg = k*Taz + 1
              
              ! Local indices on corner of proc grid
              iL = (i-i1p)*Tax + 1
              jL = (j-j1p)*Tay + 1
              kL = (k-k1p)*Taz + 1
              
              !  Write to Gdata buffer
              SELECT CASE(ftype)
              CASE ('grd')
                 Gdata(x1off:xnoff,y1off:ynoff,z1off:znoff,1:3) = GRDdata
              CASE DEFAULT  !! ('res' .OR. 'vis')
                 Gdata(x1off:xnoff,y1off:ynoff,z1off:znoff,:) = TMPiodata 
              END SELECT
           END IF


           ! MPI-BARRIER and REDUCE ALL TO MASTER.. THEN WRITE
           SELECT CASE(ftype)
           
              CASE('vis')
                 ! Slight re-ordering for vis data...
                 VISdata(:,:,:,1:5)  = Gdata(iL:iL+Tax-1,jL:jL+Tay-1,kL:kL+Taz-1,1:5)  ! u,v,w,rho,e
                 VISdata(:,:,:,12)   = Gdata(iL:iL+Tax-1,jL:jL+Tay-1,kL:kL+Taz-1,6)    ! Y
                 VISdata(:,:,:,6:11) = Gdata(iL:iL+Tax-1,jL:jL+Tay-1,kL:kL+Taz-1,7:12) ! p,T,c,mu,bulk,ktc

                 Ndata = Tax*Tay*Taz*varDIM
                 CALL MPI_BARRIER(xyzcom,mpierr)
                 CALL MPI_REDUCE(VISdata,Rdata,Ndata,MPI_REAL4,MPI_SUM,0,xyzcom,mpierr)
                 CALL MPI_BARRIER(xyzcom,mpierr)
                 
                 IF (xyzcom_id == 0) THEN
                    ! WRITE visualization data
                    !PRINT*, MAXVAL( Gdata(:,:,:,8) )
                    OPEN(UNIT=punit,FILE=TRIM(procdir),FORM='UNFORMATTED',STATUS='REPLACE')
                    DO v=1,varDIM
                       WRITE(punit) Rdata(:,:,:,v)
                    END DO
                    CLOSE(punit)
                 END IF
                 
              CASE('res')
                 RESdata = Gdata(iL:iL+Tax-1,jL:jL+Tay-1,kL:kL+Taz-1,:)

                 Ndata = Tax*Tay*Taz*varDIM
                 CALL MPI_BARRIER(xyzcom,mpierr)
                 CALL MPI_REDUCE(RESdata,Ddata,Ndata,MPI_DOUBLE_PRECISION,MPI_SUM,0,xyzcom,mpierr)
                 CALL MPI_BARRIER(xyzcom,mpierr)

                 IF (xyzcom_id == 0) THEN
                    ! WRITE restart data
                    OPEN(UNIT=punit,FILE=TRIM(procdir),FORM='UNFORMATTED',STATUS='UNKNOWN',ACTION='WRITE')
                    WRITE(punit) Ddata
                    CLOSE(punit)
                 END IF
              CASE('grd')
                 MSHdata = Gdata(iL:iL+Tax-1,jL:jL+Tay-1,kL:kL+Taz-1,1:3)

                 Ndata = Tax*Tay*Taz*3
                 CALL MPI_BARRIER(xyzcom,mpierr)
                 CALL MPI_REDUCE(MSHdata,Rdata,Ndata,MPI_REAL4,MPI_SUM,0,xyzcom,mpierr)
                 CALL MPI_BARRIER(xyzcom,mpierr)

                 IF (xyzcom_id == 0) THEN
                    ! WRITE restart data
                    OPEN(UNIT=punit,FILE=TRIM(procdir),FORM='UNFORMATTED',STATUS='UNKNOWN',ACTION='WRITE')
                    WRITE(punit) Rdata(:,:,:,1)  ! X-location
                    WRITE(punit) Rdata(:,:,:,2)  ! Y-location
                    WRITE(punit) Rdata(:,:,:,3)  ! Z-location
                    CLOSE(punit)
                 END IF

           END SELECT



        END DO
     END DO
  END DO

  ! Set to the global iodata variable pointer array



  ! Free the memory
  SELECT CASE(ftype)
      CASE('res')
         DEALLOCATE(Gdata,RESdata,Ddata)
      CASE('vis')
         DEALLOCATE(Gdata,VISdata,Rdata)
  END SELECT

  IF (xyzcom_id == 0) PRINT*, 'Done writing to ',TRIM(vizdir)

END SUBROUTINE write_data


SUBROUTINE trilinear(Tvar,Txc,Tyc,Tzc,Tax,Tay,Taz,Fvar,Fxc,Fyc,Fzc,Fax,Fay,Faz)
USE globals, ONLY: x1proc,xnproc,y1proc,ynproc,z1proc,znproc
USE interfaces, ONLY: ghost
USE interp3d_data, ONLY: varDIM,buff
USE mpi, ONLY: xyzcom_id,xyzcom,mpierr
IMPLICIT NONE
INTEGER, INTENT(IN) :: Tax,Tay,Taz,Fax,Fay,Faz
DOUBLE PRECISION, DIMENSION(Tax,Tay,Taz), INTENT(IN) :: Txc,Tyc,Tzc
DOUBLE PRECISION, DIMENSION(Tax,Tay,Taz,varDIM), INTENT(OUT) :: Tvar   ! This is the only output
DOUBLE PRECISION, DIMENSION(Fax,Fay,Faz), INTENT(IN) :: Fxc,Fyc,Fzc
DOUBLE PRECISION, DIMENSION(Fax,Fay,Faz,varDIM), INTENT(IN) :: Fvar

INTEGER :: i,j,k,ii,jj,kk,cmax,count,errorFlag,pmax,pount
INTEGER :: IBEGIN,IEND,ISTEP,JBEGIN,JEND,JSTEP
DOUBLE PRECISION :: anew,bnew,gnew
DOUBLE PRECISION :: x,y,z,da,db,dg
DOUBLE PRECISION :: aax,bbx,ccx,ffx
DOUBLE PRECISION :: aay,bby,ccy,ffy
DOUBLE PRECISION :: aaz,bbz,ccz,ffz
DOUBLE PRECISION, DIMENSION (3,3) :: AAA,iAAA
DOUBLE PRECISION, DIMENSION (3) :: BBB,XXX
DOUBLE PRECISION, DIMENSION (varDIM) :: VARS
DOUBLE PRECISION, DIMENSION (2,2,2,varDIM) :: fVARS
DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE :: pFxc,pFyc,pFzc,GsTmp
DOUBLE PRECISION, DIMENSION(:,:,:,:), ALLOCATABLE :: pFvar
INTEGER :: x1buff,xnbuff,y1buff,ynbuff,z1buff,znbuff

! Do some padding incase patch exceeds bounds
x1buff = buff; xnbuff = buff; y1buff = buff; ynbuff = buff; z1buff = buff; znbuff = buff;

IF(x1proc) x1buff = 0
IF(xnproc) xnbuff = 0
IF(y1proc) y1buff = 0
IF(ynproc) ynbuff = 0
IF(z1proc) z1buff = 0
IF(znproc) znbuff = 0

ALLOCATE( pFxc(1-x1buff:Fax+xnbuff,1-y1buff:Fay+ynbuff,1-z1buff:Faz+znbuff))
ALLOCATE( pFyc(1-x1buff:Fax+xnbuff,1-y1buff:Fay+ynbuff,1-z1buff:Faz+znbuff))
ALLOCATE( pFzc(1-x1buff:Fax+xnbuff,1-y1buff:Fay+ynbuff,1-z1buff:Faz+znbuff))
ALLOCATE( pFvar(1-x1buff:Fax+xnbuff,1-y1buff:Fay+ynbuff,1-z1buff:Faz+znbuff,varDIM))

ALLOCATE( GsTmp(1-buff:Fax+buff,1-buff:Fay+buff,1-buff:Faz+buff))

GsTmp(1:Fax,1:Fay,1:Faz) = Fxc
CALL ghost(GsTmp)
pFxc = GsTmp(1-x1buff:Fax+xnbuff,1-y1buff:Fay+ynbuff,1-z1buff:Faz+znbuff)

GsTmp(1:Fax,1:Fay,1:Faz) = Fyc
CALL ghost(GsTmp)
pFyc = GsTmp(1-x1buff:Fax+xnbuff,1-y1buff:Fay+ynbuff,1-z1buff:Faz+znbuff)

GsTmp(1:Fax,1:Fay,1:Faz) = Fzc
CALL ghost(GsTmp)
pFzc = GsTmp(1-x1buff:Fax+xnbuff,1-y1buff:Fay+ynbuff,1-z1buff:Faz+znbuff)

DO i=1,varDIM
      GsTmp(1:Fax,1:Fay,1:Faz) = Fvar(:,:,:,i)
      CALL ghost(GsTmp)
      pFvar(:,:,:,i) = GsTmp(1-x1buff:Fax+xnbuff,1-y1buff:Fay+ynbuff,1-z1buff:Faz+znbuff)
END DO


ii=1;jj=1;kk=1;

anew = 0; bnew = 0; gnew = 0;

cmax = 20;
pmax = 20;

IBEGIN = Tax;
IEND   = 1;
ISTEP  = -1;
JBEGIN = Tay;
JEND   = 1;
JSTEP  = -1;

DO k=1,Taz
    IF (xyzcom_id == 0) PRINT*,'Interp is ',100.*REAL(k)/REAL(Taz),'% done'
    JBEGIN = Tay + 1 - JBEGIN;
    JEND   = Tay + 1 - JEND;
    JSTEP  = -JSTEP;
    DO j = JBEGIN, JEND, JSTEP
        IBEGIN = Tax + 1 - IBEGIN;
        IEND   = Tax + 1 - IEND;
        ISTEP  = -ISTEP;
        DO i = IBEGIN, IEND, ISTEP

            x = Txc(i,j,k);
            y = Tyc(i,j,k);
            z = Tzc(i,j,k);
            
            DO count=1,cmax
            
                da=1;db=1;dg=1;
                anew = 0; bnew = 0; gnew = 0;
                pount = 0
                DO WHILE ( ABS(da + db + dg) > .0001 .AND. pount < pmax)
                    pount = pount + 1
                    CALL get_Istep(aax,bbx,ccx,ffx,pFxc(ii:ii+1,jj:jj+1,kk:kk+1),x,anew,bnew,gnew)
                    CALL get_Istep(aay,bby,ccy,ffy,pFyc(ii:ii+1,jj:jj+1,kk:kk+1),y,anew,bnew,gnew)
                    CALL get_Istep(aaz,bbz,ccz,ffz,pFzc(ii:ii+1,jj:jj+1,kk:kk+1),z,anew,bnew,gnew)

                    AAA(1,:) = (/ AAx, BBx, CCx /)
                    AAA(2,:) = (/ AAy, BBy, CCy /)
                    AAA(3,:) = (/ AAz, BBz, CCz /)


                    if (xyzcom_id == 0) then
                    end if

                    !CALL FINDInv(AAA, iAAA, 3, errorflag)
                    CALL getinv(AAA, iAAA, errorflag)
                    IF (xyzcom_id == 0) THEN
                    IF (errorflag == -1) print*,'Singular'
                    END IF
                    IF (errorflag == -1) EXIT

                    BBB = -(/ FFx, FFy, FFz /)
                    !XXX =  MATMUL(AAA,BBB)

                    !da = XXX(1);
                    !db = XXX(2);
                    !dg = XXX(3);
                    
                    da = -(iAAA(1,1)*FFx + iAAA(1,2)*FFy + iAAA(1,3)*FFz) 
                    db = -(iAAA(2,1)*FFx + iAAA(2,2)*FFy + iAAA(2,3)*FFz) 
                    dg = -(iAAA(3,1)*FFx + iAAA(3,2)*FFy + iAAA(3,3)*FFz) 

                    !if (xyzcom_id == 0) then
                       !print*,'============== pcount ( ',pount,' ) ===================='
                       !print*,'===== A ====='
                       !print*,real(AAA(1,:))
                       !print*,real(AAA(2,:))
                       !print*,real(AAA(3,:))
                       !print*,'===== inv(A) ====='
                       !print*,real(iAAA(1,:))
                       !print*,real(iAAA(2,:))
                       !print*,real(iAAA(3,:))
                       !print*,'===== B ====='
                       !print*,real(BBB)
                       !print*,'===== X ====='
                       !print*,real(da),real(db),real(dg)
                    !end if

                    anew = anew + da;
                    bnew = bnew + db;
                    gnew = gnew + dg;
                
                END DO
                
                !if (xyzcom_id == 0) then
                 !  print*,'=========== Counter ( ',count,' ) ========='
                 !  print*,anew,bnew,gnew
                 !  print*,ii,jj,kk
                 !  print*,real(x),real(y),real(z)
                 !  print*,real(pFxc(ii,jj,kk)),real(pFyc(ii,jj,kk)),real(pFzc(ii,jj,kk))
                !end if

                !! Decide what direction to move
                IF (anew < -1.0)  ii = max( 1-x1buff, ii-1);
                IF (anew >  1.0)  ii = min(Fax-1+xnbuff, ii+1);
                IF (bnew < -1.0)  jj = max( 1-y1buff, jj-1);
                IF (bnew >  1.0)  jj = min(Fay-1+ynbuff, jj+1);
                IF (gnew < -1.0)  kk = max( 1-z1buff, kk-1);
                IF (gnew >  1.0)  kk = min(Faz-1+znbuff, kk+1);
                
                IF ( (abs(anew) <= 1 .AND. abs(bnew) <= 1 .AND. abs(gnew) <= 1) ) THEN
                    ! Get the interpolanted value.
                    !fVARS =  pFvar(ii:ii+1,jj:jj+1,kk:kk+1,:)
                    !CALL get_Ival( VARS , fVARS , anew , bnew , gnew )
                    Tvar(i,j,k,:) = VARS
                    EXIT
                END IF
                
                IF (count >= cmax) THEN
                    !fVARS =  pFvar(ii:ii+1,jj:jj+1,kk:kk+1,:)
                    !CALL get_Ival( VARS , fVARS , anew , bnew , gnew )
                    Tvar(i,j,k,:) = VARS
                    EXIT
                END IF                    
                
            END DO
        END DO
    END DO
END DO


CALL MPI_BARRIER(xyzcom,mpierr)
if (xyzcom_id == 0) print*,'Done with interpolation routine'

END SUBROUTINE


SUBROUTINE get_Istep(aa,bb,cc,ff,Xg,xp,a,b,g)
USE constants, ONLY: eight
IMPLICIT NONE
DOUBLE PRECISION, INTENT(OUT) :: aa,bb,cc,ff
DOUBLE PRECISION, INTENT(IN) :: xp,a,b,g
DOUBLE PRECISION, DIMENSION(2,2,2), INTENT(IN) :: Xg
INTEGER :: i,j,k
DOUBLE PRECISION :: x1,x2,x3,x4,x5,x6,x7,x8
DOUBLE PRECISION :: f0,f1,f2,f3,f4,f5,f6,f7


i=1;j=1;k=1;

x1 = Xg(i  ,j  ,k  );
x2 = Xg(i+1,j  ,k  );
x3 = Xg(i  ,j+1,k  );
x4 = Xg(i+1,j+1,k  );
x5 = Xg(i  ,j  ,k+1);
x6 = Xg(i+1,j  ,k+1);
x7 = Xg(i  ,j+1,k+1);
x8 = Xg(i+1,j+1,k+1);


f0 = (x8 + x7 + x6 + x5 + x4 + x3 + x2 + x1) / eight - xp;
f1 = (x8 - x7 + x6 - x5 + x4 - x3 + x2 - x1) / eight;
f2 = (x8 + x7 - x6 - x5 + x4 + x3 - x2 - x1) / eight;
f3 = (x8 + x7 + x6 + x5 - x4 - x3 - x2 - x1) / eight;
f4 = (x8 - x7 - x6 + x5 + x4 - x3 - x2 + x1) / eight;
f5 = (x8 - x7 + x6 - x5 - x4 + x3 - x2 + x1) / eight;
f6 = (x8 + x7 - x6 - x5 - x4 - x3 + x2 + x1) / eight;
f7 = (x8 - x7 - x6 + x5 - x4 + x3 + x2 - x1) / eight;

ff = f0     + f1*a   + f2*b   + f3*g + &
       f4*a*b + f5*a*g + f6*b*g + f7*a*b*g;


aa = f1 + f4*b + f5*g + f7*b*g;
bb = f2 + f4*a + f6*g + f7*a*g;
cc = f3 + f5*a + f6*b + f7*a*b;

END SUBROUTINE

SUBROUTINE get_Ival(val,Vg,a,b,g)
USE constants, ONLY: eight
USE interp3d_data, ONLY: varDIM
IMPLICIT NONE
DOUBLE PRECISION, DIMENSION(varDIM), INTENT(OUT) :: val
DOUBLE PRECISION, DIMENSION(2,2,2,varDIM), INTENT(IN) :: Vg
DOUBLE PRECISION, INTENT(IN) :: a,b,g
DOUBLE PRECISION, DIMENSION(varDIM) :: v1,v2,v3,v4,v5,v6,v7,v8
INTEGER :: i,j,k

i=1;j=1;k=1;

v1 = Vg(i  ,j  ,k  ,:)
v2 = Vg(i+1,j  ,k  ,:)
v3 = Vg(i  ,j+1,k  ,:)
v4 = Vg(i+1,j+1,k  ,:)
v5 = Vg(i  ,j  ,k+1,:)
v6 = Vg(i+1,j  ,k+1,:)
v7 = Vg(i  ,j+1,k+1,:)
v8 = Vg(i+1,j+1,k+1,:)


val  =( (1-g)*(1-b)*(1-a)*v1 + (1-g)*(1-b)*(1+a)*v2 + &
      (1-g)*(1+b)*(1-a)*v3 + (1-g)*(1+b)*(1+a)*v4 + &
      (1+g)*(1-b)*(1-a)*v5 + (1+g)*(1-b)*(1+a)*v6 + &
      (1+g)*(1+b)*(1-a)*v7 + (1+g)*(1+b)*(1+a)*v8) / eight;

END SUBROUTINE 


!Subroutine to find the inverse of a square matrix
SUBROUTINE FINDInv(matrix, inverse, n, errorflag)
	IMPLICIT NONE
	!Declarations
	INTEGER, INTENT(IN) :: n
	INTEGER, INTENT(OUT) :: errorflag  !Return error status. -1 for error, 0 for normal
	REAL*8, INTENT(IN), DIMENSION(n,n) :: matrix  !Input matrix
	REAL*8, INTENT(OUT), DIMENSION(n,n) :: inverse !Inverted matrix
	
	LOGICAL :: FLAG = .TRUE.
	INTEGER :: i, j, k, l
	REAL :: m
	REAL, DIMENSION(n,2*n) :: augmatrix !augmented matrix
	
	!Augment input matrix with an identity matrix
	DO i = 1, n
		DO j = 1, 2*n
			IF (j <= n ) THEN
				augmatrix(i,j) = matrix(i,j)
			ELSE IF ((i+n) == j) THEN
				augmatrix(i,j) = 1
			Else
				augmatrix(i,j) = 0
			ENDIF
		END DO
	END DO
	
	!Reduce augmented matrix to upper traingular form
	DO k =1, n-1
		IF (augmatrix(k,k) == 0) THEN
			FLAG = .FALSE.
			DO i = k+1, n
				IF (augmatrix(i,k) /= 0) THEN
					DO j = 1,2*n
						augmatrix(k,j) = augmatrix(k,j)+augmatrix(i,j)
					END DO
					FLAG = .TRUE.
					EXIT
				ENDIF
				IF (FLAG .EQV. .FALSE.) THEN
					!PRINT*, "Matrix is non - invertible"
					inverse = 0
					errorflag = -1
					return
				ENDIF
			END DO
		ENDIF
		DO j = k+1, n			
			m = augmatrix(j,k)/augmatrix(k,k)
			DO i = k, 2*n
				augmatrix(j,i) = augmatrix(j,i) - m*augmatrix(k,i)
			END DO
		END DO
	END DO
	
	!Test for invertibility
	DO i = 1, n
		IF (augmatrix(i,i) == 0) THEN
			!PRINT*, "Matrix is non - invertible"
			inverse = 0
			errorflag = -1
			return
		ENDIF
	END DO
	
	!Make diagonal elements as 1
	DO i = 1 , n
		m = augmatrix(i,i)
		DO j = i , (2 * n)				
			   augmatrix(i,j) = (augmatrix(i,j) / m)
		END DO
	END DO
	
	!Reduced right side half of augmented matrix to identity matrix
	DO k = n-1, 1, -1
		DO i =1, k
		m = augmatrix(i,k+1)
			DO j = k, (2*n)
				augmatrix(i,j) = augmatrix(i,j) -augmatrix(k+1,j) * m
			END DO
		END DO
	END DO				
	
	!store answer
	DO i =1, n
		DO j = 1, n
			inverse(i,j) = augmatrix(i,j+n)
		END DO
	END DO
	errorflag = 0
END SUBROUTINE FINDinv


subroutine getinv(M,M1,err)
implicit none
double precision, dimension(3,3), intent(in) :: M
double precision, dimension(3,3), intent(out) :: M1
double precision, dimension(3,3) :: MTEMP
integer, intent(out) :: err
double precision :: det

err = 0
det = M(1,1)*( M(2,2)*M(3,3)-M(2,3)*M(3,2) ) - M(1,2)*( M(2,1)*M(3,3) - M(3,1)*M(2,3) ) & 
     + M(1,3)*( M(2,1)*M(3,2) - M(3,1)*M(2,2) )

if ( DABS(det) < 1.0D-16 ) then
   M1 = 0.0D0
   err = -1
   return
end if

MTEMP(1,1) =  (M(2,2)*M(3,3) - M(3,2)*M(2,3))
MTEMP(1,2) = -(M(1,2)*M(3,3) - M(3,2)*M(1,3))
MTEMP(1,3) =  (M(1,2)*M(2,3) - M(2,2)*M(1,3))
MTEMP(2,1) = -(M(2,1)*M(3,3) - M(3,1)*M(2,3))
MTEMP(2,2) =  (M(1,1)*M(3,3) - M(3,1)*M(1,3))
MTEMP(2,3) = -(M(1,1)*M(2,3) - M(2,1)*M(1,3))
MTEMP(3,1) =  (M(2,1)*M(3,2) - M(3,1)*M(2,2))
MTEMP(3,2) = -(M(1,1)*M(3,2) - M(3,1)*M(1,2))
MTEMP(3,3) =  (M(1,1)*M(2,2) - M(2,1)*M(1,2))

M1 = MTEMP/det

end subroutine getinv
