! you must re-define this module name to be unique across all problems
MODULE post_proc_data
  USE mpi
  SAVE
  ! ---------------------------------------------------------------------------
  ! Place Problem dependent state data here, similar to inputs.f and globals.f
  ! ---------------------------------------------------------------------------  
  DOUBLE PRECISION      :: Tinitial = 300.0D0
  DOUBLE PRECISION      :: Tout = 500.0D0    ! This is the correction factor for the outflow
  DOUBLE PRECISION      :: Pinitial = 1.0D6
  CHARACTER(LEN=90) 	:: restart_dir  = '/path/to/directory'
  CHARACTER(LEN=90) 	:: gridfile    = 'nozzle.grid'
  CHARACTER(LEN=90) 	:: dump_file  = '/path/to/directory'
  INTEGER               :: Tcount = 1                   ! Start index for this loop <= Tstart
  INTEGER               :: Tstart = 1                   ! Start index for database
  INTEGER               :: Tend = 10                    ! End index for file loop
  INTEGER               :: varDIM = 12                  ! Variable dimension of the incoming data array (12 is for compressible, single fluid)
  DOUBLE PRECISION      :: Re_BL = 10.0D3
  DOUBLE PRECISION      :: del_BL = 2.0D0
  CHARACTER(LEN=90) 	:: flowfile    = 'nozzle.inflow'
  INTEGER               :: Tpx = 1                      ! Processors in x-direction for dataset
  INTEGER               :: Tpy = 1                      ! Processors in y-direction for dataset
  INTEGER               :: Tpz = 1                      ! Processors in z-direction for dataset
  CHARACTER(LEN=3)      :: ftype = 'res'                ! Select file type: 'res'-restart, 'vis'-viz file
  CHARACTER(LEN=90) 	:: post_type = 'normal'         ! Select the post-proc routines to run

  CHARACTER(LEN=90) 	:: MIRname = 'post'
  CHARACTER(LEN=90) 	:: Gname = 'grid'
  INTEGER               :: ixwI = 300                   ! For the walls case, the subdomain in x-global coords (min)
  INTEGER               :: ixwF = 500                   ! For the walls case, the subdomain in x-global coords (max)

  ! Not set in Namelist, these variables are read in and used to set the inflow BC
  CHARACTER(LEN=30) 	:: prob_jobname = 'post_proc'
  INTEGER               :: Tax, Tay, Taz                ! Grid point per-proc for dataset
  INTEGER               :: i1g,ifg,j1g,jfg,k1g,kfg      ! Global index numbers of this proc.
  INTEGER               :: i1p,ifp,j1p,jfp,k1p,kfp      ! Global proc index range for this chunk
  INTEGER               :: T_iter

  DOUBLE PRECISION      :: Mach,NPR,P_in,rho_in
  DOUBLE PRECISION      :: mu_0

  INTEGER, DIMENSION(:,:), ALLOCATABLE :: procmap
  INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: invprocmap


  !! Correlation stuff
  CHARACTER(LEN=90) 	:: corrfile    = 'corr.dat'
  INTEGER               :: ncorr = 5
  INTEGER, DIMENSION(:), ALLOCATABLE :: corrX,corrY
  CHARACTER(LEN=90) 	:: corrV = 'u'
  CHARACTER(LEN=90) 	:: cortag = 'spec-u-'
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: corr_data


  !! Flags
  LOGICAL :: corr_flag = .FALSE.
  LOGICAL :: plane_flag = .TRUE.
  LOGICAL :: mir_flag = .FALSE.

END MODULE post_proc_data

! -----------------------------------------------------------------------------
! prob_inputs
! Read problem-dependent inputs and save in module variables
! must also define 'jobname'
! Called by the Miranda main program.
! -----------------------------------------------------------------------------
SUBROUTINE prob_inputs(fileName)
  USE mpi
  USE prob_interface, ONLY : jobname
  USE post_proc_data
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
  NAMELIST /post_proc_vars/ restart_dir,gridfile,flowfile,Tcount,Tstart,Tend,varDIM,Re_BL,del_BL, & 
      ftype,Tinitial,Pinitial,corrfile,corr_flag,cortag,corrV,plane_flag,mir_flag,post_type,ixwI,ixwF,Tout, &
      dump_file

  ! Uncomment to open and read a namelist file
  OPEN(UNIT=inputUnit,FILE=TRIM(fileName),FORM='FORMATTED',STATUS='OLD')
  READ(UNIT=inputUnit,NML=post_proc_vars)
  CLOSE(inputUnit)
  
  ! give the right name
  jobname = TRIM(prob_jobname)
  
      

  ! Read in the procmap file and set Tpx,Tpy,Tpz
  dsI = 5
  nx = 0
  ny = 0
  nz = 0
  WRITE(plotmir,'(2A)') TRIM(restart_dir),'/plot.mir'
  OPEN(UNIT=inputUnit,FILE=TRIM(plotmir),FORM='FORMATTED',STATUS='OLD')
  DO i=1,dsI
      READ(inputUnit,*) comments
  END DO
  READ(inputUnit,*) comments, nx, ny, nz
  READ(inputUnit,*) comments, bsx,bsy,bsz
  CLOSE(inputUnit)  

  Tpx = nx / bsx
  Tpy = ny / bsy
  Tpz = nz / bsz  


  ! Read the inflow boundary conditions                                     
  in_unit = 26
  OPEN(UNIT=in_unit, file=TRIM(flowfile))
  READ(in_unit,*) NPR,P_in,rho_in,Mach
  CLOSE(in_unit)


  
  ! Re-arrange the correlation data  
  IF(corr_flag) THEN
  in_unit = 26
  OPEN(UNIT=in_unit, file=TRIM(corrfile))
  READ(in_unit,*) ncorr
  ALLOCATE(corrX(ncorr))
  ALLOCATE(corrY(ncorr))
  DO i=1,ncorr
      READ(in_unit,*) corrX(i),corrY(i)
  END DO
  CLOSE(in_unit)
  END IF


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
  USE post_proc_data
  USE metrics, ONLY: x_c
  USE interfaces, ONLY: newdir, restart, viz
  IMPLICIT NONE
  INTEGER :: i,j,funit,isync
  CHARACTER(LEN=90) :: flipfile,fform,plotdir
  LOGICAL :: fexist
  DOUBLE PRECISION :: Rgas,U_in,e_in,T_in
  INTEGER :: TFstart,TFcount,TFend,Nstride,Nsten,ii,TS_old
  DOUBLE PRECISION :: sig,Gnorm
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: Gauss
  DOUBLE PRECISION, DIMENSION(:,:,:,:), ALLOCATABLE :: fil_data
  DOUBLE PRECISION, DIMENSION(ax,ay,az) :: dum,dudy
  DOUBLE PRECISION :: Tin

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


  ! Read in arbitrary restart file -or- vis file
  ! nx,ny,nz must be the same for both cases... check that here
  CALL setup_dataread()  

  

  SELECT CASE(post_type)
      
      !! Modify the data for this viz-file and then dump a restart file.
      CASE('moddata')

         CALL read_data(Tcount)

         !! Fix the temperature so the outflow is always 300K
         dum = T
         T = ( dum - T_in ) * ( one / (Tout-T_in) )
         dum = T
         T = dum * ( Tinitial - T_in ) + T_in

         rho = p / Rgas / T
         e = (p/(gamma-one))/rho


         ! Dump restart data here...
         WRITE(iodir,'(2A,I4.4)') TRIM(jobdir),'/res',0
         IF (world_id==master) WRITE(6,'(2A)') ' Writing data to ',TRIM(iodir)
         CALL newdir(LEN_TRIM(iodir),TRIM(iodir),isync)
         CALL restart('w',TRIM(iodir),iodata(:,:,:,1:nres))


         ! Dump the viz, to test
         WRITE(iodir,'(2A)') TRIM(jobdir),'/grid'
         CALL newdir(LEN_TRIM(iodir),TRIM(iodir),isync)
         CALL grid()


         WRITE(iodir,'(2A,I4.4)') TRIM(jobdir),'/vis',0
         IF (world_id==master) WRITE(6,'(2A)') ' Writing data to ',TRIM(iodir)
         CALL newdir(LEN_TRIM(iodir),TRIM(iodir),isync)
         CALL MPI_BARRIER(MPI_COMM_WORLD,mpierr)
         CALL viz(TRIM(iodir),iodata)




      CASE('walls')

         CALL get_walls


      CASE('filter')
         MIRname = 'filter'

         ! Set up filter array
         ALLOCATE(fil_data(size(iodata,1),size(iodata,2),size(iodata,3),size(iodata,4)))

         ! Loop over the new filtered data set
         Nstride = 5
         Nsten = 3   ! Point to left and right... total stencil is 2*N+1
         TFstart = 1
         TFcount = Tcount
         TFend = FLOOR( dble(Tend - Tstart - (2*Nsten+1 ))/ dble(Nstride) )

         !! Set up the Gaussian over Nsten points
         ALLOCATE(Gauss(2*Nsten+1))
         sig = - DBLE(Nsten**2)/ LOG( .01D0 )
         DO i=1,2*Nsten+1
            Gauss(i) = exp( - DBLE(i-(Nsten+1))**2 / sig )
         END DO
         Gnorm = SUM( Gauss )
         Gauss = Gauss / Gnorm

         Tend = TFend
         TS_old = Tstart
         Tstart = TFstart
         tviz = tviz*Nstride


         DO i=TFcount,TFend

            T_iter = i

            fil_data = zero
            ! Loop over points in Gaussian
            DO j=1,2*Nsten+1

               ii = TS_old + Nstride*i + (j - (Nsten+1) )
               IF(xyzcom_id == 0) THEN
                  PRINT*,ii,(j - (Nsten+1) )
               END IF

               CALL read_data(ii)
               fil_data = fil_data + iodata * Gauss(j)

            END DO

            iodata = fil_data

            WRITE(plotdir,'(3A,I4.4)') TRIM(jobdir),'/',TRIM(MIRname),i
            CALL newdir(LEN_TRIM(plotdir),TRIM(plotdir),isync)
            CALL prob_plots(plotdir)

         END DO

         IF(xyzcom_id == 0) PRINT*,Gauss,SUM(Gauss)


      CASE DEFAULT
         ! Loop over the available files and perform operations on them as needed.
         MIRname = 'post'
         DO i=Tcount,Tend
            T_iter = i

            IF (xyzcom_id == 0) PRINT*,'Reading file.....# ',i
            CALL read_data(i)
            
            IF (xyzcom_id == 0) PRINT*,'Processing file..# ',i
            WRITE(plotdir,'(3A,I4.4)') TRIM(jobdir),'/',TRIM(MIRname),i
            CALL newdir(LEN_TRIM(plotdir),TRIM(plotdir),isync)
            
            CALL prob_plots(plotdir)
            
         END DO

      END SELECT


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


SUBROUTINE get_walls
  USE mpi
  USE constants
  USE inputs
  USE globals
  USE extend
  USE post_proc_data
  USE metrics, ONLY: x_c,y_c,z_c
  USE interfaces, ONLY: newdir,grad
  IMPLICIT NONE
  INTEGER :: i,j,funit,isync,nx_wall,tt,ii
  LOGICAL :: proc_side
  INTEGER :: indx_side,iL,iG
  CHARACTER(LEN=90) :: plotdir,griddir
  DOUBLE PRECISION, DIMENSION(ax,ay,az) :: dum,dudy,dTdy,vtmp
  DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE :: wall,wall2,wall_local
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: walltmp
  INTEGER :: nvars,T_tmp,cdim
  CHARACTER(LEN=90), DIMENSION(:), ALLOCATABLE :: vars
  INTEGER, DIMENSION(:), ALLOCATABLE :: minvar,maxvar


      nx_wall = ixwF - ixwI + 1
      nvars = 4

      ALLOCATE(wall(nx_wall,nz,nvars+3))
      ALLOCATE(wall2(nx_wall,nz,nvars+3))
      ALLOCATE(wall_local(ax,az,nvars+3))
      ALLOCATE(walltmp(nx_wall,nz))

      ALLOCATE(vars(nvars))
      ALLOCATE(minvar(50))
      ALLOCATE(maxvar(50))

      DO tt=Tcount,Tend
         T_iter = tt

         IF (xyzcom_id == 0) PRINT*,'Reading file.....# ',tt
         CALL read_data(tt)
            
         IF (xyzcom_id == 0) PRINT*,'Processing file..# ',tt


         CALL grad(u,dum,dudy,dum)
         CALL grad(T,dum,dTdy,dum)

         DO j=1,2  !! Get top and bottom walls

           SELECT CASE(j)
           CASE(1)
              proc_side = ynproc
              indx_side = ay
              MIRname = 'walls_top'
              IF(xyzcom_id==0) PRINT*,'Getting top wall'
           CASE(2)
              proc_side = y1proc
              indx_side = 1
              MIRname = 'walls_bot'
              IF(xyzcom_id==0) PRINT*,'Getting bottom wall'
           CASE DEFAULT
              IF(xyzcom_id==0) PRINT*,"No side picked for wall output"
              EXIT
           END SELECT


           WRITE(plotdir,'(3A,I4.4)') TRIM(jobdir),'/',TRIM(MIRname),tt
           CALL newdir(LEN_TRIM(plotdir),TRIM(plotdir),isync)

           wall = zero
           wall2 = zero
           
           IF (proc_side) THEN
              wall_local(:,:,1) = p(:,indx_side,:)
              wall_local(:,:,2) = T(:,indx_side,:)
              wall_local(:,:,3) = dudy(:,indx_side,:)
              wall_local(:,:,4) = dTdy(:,indx_side,:)
              wall_local(:,:,5) = x_c(:,indx_side,:)
              wall_local(:,:,6) = y_c(:,indx_side,:)
              wall_local(:,:,7) = z_c(:,indx_side,:)

              ! Loop over the global ix of sub-data set
              DO ii=ixwI,ixwF
                 IF( ii >= ix(1) .AND. ii <= ix(ax)) THEN ! In range on this proc
                    iL = ii - ix(1) + 1
                    iG = ii - ixwI + 1
                    wall( iG,iz,:) = wall_local( iL, :, :)
                    !PRINT*,iL,iG,ii
                 END IF
              END DO

           END IF

           !! Sum this array up... put it on wall2 array
           cdim = nx_wall*nz*(nvars+3)
           !IF(xycom_id==0) PRINT*,cdim,size(wall),size(wall2)
           CALL MPI_REDUCE(wall,wall2,cdim,MPI_DOUBLE_PRECISION,MPI_SUM,0,xyzcom,mpierr)
           !CALL MPI_BARRIER(xyzcom)


           i = 1
           vars(i) = 'Pwall'
           walltmp = wall2(:,:,i)
           minvar(i) = minval(walltmp); maxvar(i) = maxval(walltmp)
           CALL write_VIZ2d(walltmp,nx_wall,nz,.TRUE.,plotdir)
         
           i = i+1
           vars(i) = 'Twall'
           walltmp = wall2(:,:,i)
           minvar(i) = minval(walltmp); maxvar(i) = maxval(walltmp)
           CALL write_VIZ2d(walltmp,nx_wall,nz,.FALSE.,plotdir)
           
           i = i+1
           vars(i) = 'dudy'
           walltmp = wall2(:,:,i)
           minvar(i) = minval(walltmp); maxvar(i) = maxval(walltmp)
           CALL write_VIZ2d(walltmp,nx_wall,nz,.FALSE.,plotdir)
         
           i = i+1
           vars(i) = 'dTdy'
           walltmp = wall2(:,:,i)
           minvar(i) = minval(walltmp); maxvar(i) = maxval(walltmp)
           CALL write_VIZ2d(walltmp,nx_wall,nz,.FALSE.,plotdir)

         
           !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           !! Count up the total variables and    !!
           !!  write the *.mir file               !!
           !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           CALL write_MIR(vars,nvars,minvar,maxvar,nx_wall,1,nz)
           
           !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           !! Write the grid the first time !!!!!!!!
           !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           Gname = 'grid2dw'
           IF (T_iter==Tcount .AND. j==1) THEN
              WRITE(griddir,'(3A)') TRIM(plotdir),'/../',TRIM(Gname)
              CALL newdir(LEN_TRIM(griddir),TRIM(griddir),isync)
            
              i = nvars+1
              walltmp = wall2(:,:,i)
              CALL write_VIZ2d(walltmp,nx_wall,nz,.TRUE.,griddir)
              i = i+1
              walltmp = wall2(:,:,i)
              CALL write_VIZ2d(walltmp,nx_wall,nz,.FALSE.,griddir)
              i = i+1
              walltmp = wall2(:,:,i)
              CALL write_VIZ2d(walltmp,nx_wall,nz,.FALSE.,griddir)
           END IF
         
           !! IF the mir-flag is on, after the first data set, write out the entire file time in plot.mir
           IF(mir_flag) THEN
              T_tmp = T_iter
              T_iter = Tend
              CALL write_MIR(vars,nvars,minvar,maxvar,nx_wall,1,nz)
              T_iter = T_tmp
           END IF


        END DO  ! End loop over sides (top-bottom)

      END DO  ! End loop over Viz files



END SUBROUTINE


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
  USE post_proc_data
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
  

  ! Pressure down length of nozzle
  z_a = 0.0d0 !! dble(nz-1)*dz / 2.0d0 
  z_b = z_a + dble(nz-1)*dz

  where (z_c >= z_a .and. z_c <= z_a)
      box = 1.0d0
  elsewhere
      box = 0.0d0
  end where

  tmp0 = 0.0d0
  if (iy(1)==1) tmp0(:,1,:) = 1.0d0
  tmp1 = box*tmp0

  p_dwn = SUBSUM3YZ(p*tmp1)/SUBSUM3YZ(tmp1)
  xbar = SUBSUM3YZ(x_c*tmp1)/SUBSUM3YZ(tmp1)

  tmp0 = 0.0d0
  if (iy(ay)==ny) tmp0(:,ay,:) = 1.0d0
  tmp1 = box*tmp0
  p_up = SUBSUM3YZ(p*tmp1)/SUBSUM3YZ(tmp1)

  core = del_BL*2.0D0 
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


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Get planar averaged values and stores them as 2d/single proc/ *.mir files
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  IF (plane_flag) THEN
    CALL get_planes(plotdir)
  END IF
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Compute a span averaged line of data for the given variable
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  IF (corr_flag) THEN
    var = corrV
    CALL get_corr(var,plotdir)
  END IF
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE prob_plots

SUBROUTINE prob_source()
END SUBROUTINE


SUBROUTINE get_corr(var,plotdir)
  USE mpi
  USE globals, ONLY: flen,jobdir
  USE inputs, ONLY: nz
  USE post_proc_data
  IMPLICIT NONE
  CHARACTER(LEN=flen), INTENT(IN) :: var
  CHARACTER(LEN=*), INTENT(IN) :: plotdir
  CHARACTER(LEN=flen) :: plotFile
  INTEGER, PARAMETER :: plotUnit=21
  INTEGER :: i,j,k
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: corr
  CHARACTER(LEN=flen) :: myformat
  DOUBLE PRECISION, DIMENSION(nz) :: corrT
  DOUBLE PRECISION, DIMENSION(2*nz,ncorr) :: dtmp
  DOUBLE PRECISION, DIMENSION(nz,ncorr) :: tmp,Scorr
  DOUBLE PRECISION, DIMENSION(ncorr) :: ave
  
INTEGER :: nxs,nys

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Compute a span averaged line of data for the given variable
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  IF (.NOT. ALLOCATED(corr_data)) THEN
      ALLOCATE(corr_data(nz,ncorr))
      corr_data = 0.0D0
  END IF


  ALLOCATE(corr(nz,ncorr))
  DO i=1,ncorr
      nxs = corrX(i)
      nys = corrY(i)
      CALL get_span(corrT,nxs,nys,var)
      corr(:,i) = corrT
      ave(i) = SUM(corr(:,i))/DBLE(nz)
      tmp(:,i) = corr(:,i) - ave(i)
  END DO

  ! Periodic copy
  dtmp(1:nz,:) = tmp(1:nz,:)
  dtmp(nz+1:2*nz,:) = tmp(1:nz,:)

  Scorr = 0.0D0
  DO i=1,nz ! Loop over lengths  
    DO j=1,nz
        Scorr(i,:) = Scorr(i,:) + dtmp(j,:)*dtmp(j+i-1,:);        
    END DO    
  END DO

  ! Normalize the plot
  DO i=1,ncorr
      Scorr(:,i) = Scorr(:,i) / Scorr(1,i)
  END DO

  corr_data = corr_data + Scorr

  WRITE(myformat,'(A,I2,A)') '(', ncorr , 'ES12.4)'  
  ! master cpu writes plot file--------------------------------------------------------------------
  SELECT CASE(xyzcom_id)
  CASE(master)
    WRITE(plotFile,'(5A)') TRIM(plotdir),'/',TRIM(var),TRIM(cortag),'.dat'
    OPEN(UNIT=plotUnit,FILE=TRIM(plotFile),FORM='FORMATTED',STATUS='REPLACE')
    WRITE(plotUnit,*) "%# Var=",TRIM(var),', pts=',ncorr
    DO i=1,ncorr
    DO k=1,nz
      WRITE(plotUnit,TRIM(myformat)) corr(k,i)
    END DO
    END DO
    CLOSE(plotUnit)
  END SELECT
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  
  !! Also write the average of the corr_data to a file in jobdir
  IF (T_iter == Tend) THEN
      corr_data = corr_data / ( DBLE(Tend-Tcount+1) )

      SELECT CASE(xyzcom_id)
      CASE(master)
         WRITE(plotFile,'(5A)') TRIM(jobdir),'/',TRIM(var),TRIM(cortag),'total.dat'
         OPEN(UNIT=plotUnit,FILE=TRIM(plotFile),FORM='FORMATTED',STATUS='REPLACE')
         WRITE(plotUnit,*) "%# Var=",TRIM(var),', pts=',ncorr
         !DO i=1,ncorr
            DO k=1,nz
               WRITE(plotUnit,TRIM(myformat)) corr_data(k,:)
            END DO
         !END DO
         CLOSE(plotUnit)
      END SELECT

  END IF
  



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
  USE post_proc_data
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
  USE constants, ONLY : zero,one,half,three,Runiv
  USE inputs, ONLY : ns,viscous,conductive,diffusive,resistive,diffusive,gamma
  USE globals, ONLY : materials,T
  USE globals, ONLY : molwts
    USE globals, ONLY : mu,bulk,ktc,Diff,eta ! *** OUTPUTS ***
  USE leosvars, ONLY : part_r
  USE post_proc_data
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
  !ktc = cp * mu / Pr

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
 USE post_proc_data
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

  ! Make sure Jacobian is 2d... extruded mesh
  nz_tmp = nz
  nz = 1
  CALL get_jacobian()
  nz = nz_tmp

END SUBROUTINE




SUBROUTINE get_2d_z(var,plane)
      USE globals, ONLY: ax,ay,az,ix,iy,iz
      USE inputs, ONLY: nx,ny,nz
      USE interfaces, ONLY: SUBSUM3YZ
      IMPLICIT NONE
      DOUBLE PRECISION, DIMENSION(ax,ay,az), INTENT(IN) :: var
      DOUBLE PRECISION, DIMENSION(nx,ny) ,INTENT(OUT)   :: plane
      
      DOUBLE PRECISION, DIMENSION(ax,1,az) :: xlineL
      DOUBLE PRECISION, DIMENSION(nx) :: xlineG
      INTEGER :: i,ii

      DO i=1,ny
         ii = MOD( (i-1), ay ) + 1
         
         xlineL = 0.0D0
         IF( iy(ii) == i ) xlineL(:,1,:) = var(:,ii,:)
           
         xlineG = SUBSUM3YZ(xlineL)
         plane(:,i) = xlineG
      END DO

      plane = plane / dble(nz)

END SUBROUTINE


SUBROUTINE get_planes(plotdir)
  USE mpi
  USE globals
  USE inputs
  USE constants, ONLY: zero,one,two,half
  USE interfaces, ONLY: newdir,grad,Dhydro,filter
  USE metrics, ONLY: x_c,y_c,z_c,muA,muB,muC,del_A,del_B,del_C
  USE RKscheme, ONLY : nRKsteps
  USE post_proc_data
  IMPLICIT NONE
  CHARACTER(LEN=90), INTENT(IN) :: plotdir
  DOUBLE PRECISION, DIMENSION(50)     :: minvar,maxvar
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: mean,flc
  DOUBLE PRECISION, DIMENSION(SIZE(u,1),SIZE(u,2),SIZE(u,3)) :: vtmp,Lmean,Lflc,L1flc,L2flc,L3flc
  DOUBLE PRECISION, DIMENSION(SIZE(u,1),SIZE(u,2),SIZE(u,3)) :: rho_grad,Ptot
  CHARACTER(LEN=90), DIMENSION(:), ALLOCATABLE :: vars
  CHARACTER(LEN=90), DIMENSION(50) :: varstmp
  CHARACTER(LEN=90) :: griddir
  DOUBLE PRECISION :: ddt,st
  INTEGER :: i,nvars,isync,k,T_tmp

  !! Get the z-averaged planes and save them as *.mir format... these can be used for plotting in VisIt
  !! or used in matlab to get profiles.  The 2d files are written out on ONE proc.

  Ptot = p + half*rho*(u*u + v*v + w*w)   ! Stagnation/Total pressure
  CALL grad(rho,L1flc,L2flc,L3flc)
  rho_grad = sqrt( L1flc**2 + L2flc**2 + L3flc**2 )
  vtmp = rho_grad
  CALL filter(lfilter,vtmp,rho_grad)


  ALLOCATE(mean(nx,ny))
  ALLOCATE(flc(nx,ny))


  !! Mean variables
  i = 1
  vtmp = u
  varstmp(i) = 'U-vel'
  CALL get_2d_z(vtmp,mean)
  minvar(i) = minval(mean); maxvar(i) = maxval(mean)
  CALL write_VIZ2d(mean,nx,ny,.TRUE.,plotdir)

  i = i + 1
  vtmp = v
  varstmp(i) = 'V-vel'
  CALL get_2d_z(vtmp,mean)
  minvar(i) = minval(mean); maxvar(i) = maxval(mean)
  CALL write_VIZ2d(mean,nx,ny,.FALSE.,plotdir)

  i = i + 1
  vtmp = w
  varstmp(i) = 'W-vel'
  CALL get_2d_z(vtmp,mean)
  minvar(i) = minval(mean); maxvar(i) = maxval(mean)
  CALL write_VIZ2d(mean,nx,ny,.FALSE.,plotdir)

  i = i + 1
  vtmp = p
  varstmp(i) = 'pressure'
  CALL get_2d_z(vtmp,mean)
  minvar(i) = minval(mean); maxvar(i) = maxval(mean)
  CALL write_VIZ2d(mean,nx,ny,.FALSE.,plotdir)

  i = i + 1
  vtmp = rho
  varstmp(i) = 'density'
  CALL get_2d_z(vtmp,mean)
  minvar(i) = minval(mean); maxvar(i) = maxval(mean)
  CALL write_VIZ2d(mean,nx,ny,.FALSE.,plotdir)

  i = i + 1
  vtmp = T
  varstmp(i) = 'temperature'
  CALL get_2d_z(vtmp,mean)
  minvar(i) = minval(mean); maxvar(i) = maxval(mean)
  CALL write_VIZ2d(mean,nx,ny,.FALSE.,plotdir)
  
  !! Reynolds Stress terms
  i = i + 1
  vtmp = u
  varstmp(i) = 'uu'
  CALL get_2d_z(vtmp,mean)
  DO k=1,az; Lmean(:,:,k) = mean(ix,iy); END DO
  L1flc = vtmp - Lmean
  Lflc = L1flc**2
  CALL get_2d_z(Lflc,flc)
  minvar(i) = minval(flc); maxvar(i) = maxval(flc)
  CALL write_VIZ2d(flc,nx,ny,.FALSE.,plotdir)

  i = i + 1
  vtmp = v
  varstmp(i) = 'vv'
  CALL get_2d_z(vtmp,mean)
  DO k=1,az; Lmean(:,:,k) = mean(ix,iy); END DO
  L2flc = vtmp - Lmean
  Lflc = L2flc**2
  CALL get_2d_z(Lflc,flc)
  minvar(i) = minval(flc); maxvar(i) = maxval(flc)
  CALL write_VIZ2d(flc,nx,ny,.FALSE.,plotdir)

  i = i + 1
  vtmp = w
  varstmp(i) = 'ww'
  CALL get_2d_z(vtmp,mean)
  DO k=1,az; Lmean(:,:,k) = mean(ix,iy); END DO
  L3flc = vtmp - Lmean
  Lflc = L3flc**2
  CALL get_2d_z(Lflc,flc)
  minvar(i) = minval(flc); maxvar(i) = maxval(flc)
  CALL write_VIZ2d(flc,nx,ny,.FALSE.,plotdir)

  i = i + 1
  varstmp(i) = 'uv'
  Lflc = L1flc*L2flc
  CALL get_2d_z(Lflc,flc)
  minvar(i) = minval(flc); maxvar(i) = maxval(flc)
  CALL write_VIZ2d(flc,nx,ny,.FALSE.,plotdir)

  i = i + 1
  varstmp(i) = 'uw'
  Lflc = L1flc*L3flc
  CALL get_2d_z(Lflc,flc)
  minvar(i) = minval(flc); maxvar(i) = maxval(flc)
  CALL write_VIZ2d(flc,nx,ny,.FALSE.,plotdir)

  i = i + 1
  varstmp(i) = 'vw'
  Lflc = L2flc*L3flc
  CALL get_2d_z(Lflc,flc)
  minvar(i) = minval(flc); maxvar(i) = maxval(flc)
  CALL write_VIZ2d(flc,nx,ny,.FALSE.,plotdir)

  !! Total Pressure
  i = i + 1
  vtmp = Ptot
  varstmp(i) = 'Ptot'
  CALL get_2d_z(vtmp,mean)
  minvar(i) = minval(mean); maxvar(i) = maxval(mean)
  CALL write_VIZ2d(mean,nx,ny,.FALSE.,plotdir)  

  i = i + 1
  vtmp = Ptot
  varstmp(i) = 'PtPt'
  CALL get_2d_z(vtmp,mean)
  DO k=1,az; Lmean(:,:,k) = mean(ix,iy); END DO
  Lflc = vtmp - Lmean
  Lflc = Lflc**2
  CALL get_2d_z(Lflc,flc)
  minvar(i) = minval(flc); maxvar(i) = maxval(flc)
  CALL write_VIZ2d(flc,nx,ny,.FALSE.,plotdir)

  !! Density Gradient
  i = i + 1
  vtmp = rho_grad
  varstmp(i) = 'gradRHO'
  CALL get_2d_z(vtmp,mean)
  minvar(i) = minval(mean); maxvar(i) = maxval(mean)
  CALL write_VIZ2d(mean,nx,ny,.FALSE.,plotdir)  

  !! Get the artificial terms
  CALL compute_AFLES()

  !! Shear viscosity
  i = i + 1
  vtmp = mu
  varstmp(i) = 'mu'
  CALL get_2d_z(vtmp,mean)
  minvar(i) = minval(mean); maxvar(i) = maxval(mean)
  CALL write_VIZ2d(mean,nx,ny,.FALSE.,plotdir)  

  i = i + 1
  vtmp = muA
  varstmp(i) = 'muA'
  CALL get_2d_z(vtmp,mean)
  minvar(i) = minval(mean); maxvar(i) = maxval(mean)
  CALL write_VIZ2d(mean,nx,ny,.FALSE.,plotdir) 

  i = i + 1
  vtmp = muB
  varstmp(i) = 'muB'
  CALL get_2d_z(vtmp,mean)
  minvar(i) = minval(mean); maxvar(i) = maxval(mean)
  CALL write_VIZ2d(mean,nx,ny,.FALSE.,plotdir) 

  i = i + 1
  vtmp = muC
  varstmp(i) = 'muC'
  CALL get_2d_z(vtmp,mean)
  minvar(i) = minval(mean); maxvar(i) = maxval(mean)
  CALL write_VIZ2d(mean,nx,ny,.FALSE.,plotdir) 

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Count up the total variables and    !!
  !!  write the *.mir file               !!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  nvars = i
  ALLOCATE(vars(nvars))
  vars = varstmp(1:nvars)
  CALL write_MIR(vars,nvars,minvar,maxvar,nx,ny,1)



  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Write the grid the first time !!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  Gname = 'grid2d'
  IF (T_iter==Tcount) THEN
      WRITE(griddir,'(3A)') TRIM(plotdir),'/../',TRIM(Gname)
      CALL newdir(LEN_TRIM(griddir),TRIM(griddir),isync)

      vtmp = x_c
      CALL get_2d_z(vtmp,mean)
      CALL write_VIZ2d(mean,nx,ny,.TRUE.,griddir)
      
      vtmp = y_c
      CALL get_2d_z(vtmp,mean)
      CALL write_VIZ2d(mean,nx,ny,.FALSE.,griddir)
      
      vtmp = z_c
      CALL get_2d_z(vtmp,mean)
      CALL write_VIZ2d(mean,nx,ny,.FALSE.,griddir)
  END IF

  !! IF the mir-flag is on, after the first data set, write out the entire file time in plot.mir
  IF(mir_flag) THEN
      T_tmp = T_iter
      T_iter = Tend
      CALL write_MIR(vars,nvars,minvar,maxvar,nx,ny,1)
      T_iter = T_tmp
  END IF


END SUBROUTINE

SUBROUTINE write_VIZ2d(var,d1,d2,init,vizdir)
  USE mpi
  USE globals
  USE inputs
  USE post_proc_data
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: d1,d2
  DOUBLE PRECISION, DIMENSION(d1,d2), INTENT(IN) :: var
  LOGICAL, INTENT(IN) :: init
  CHARACTER(LEN=90), INTENT(IN) :: vizdir
  REAL(KIND=4), DIMENSION(d1,d2) :: tmp


  !! Master CPU does the writing
  IF(xyzcom_id == 0) THEN

      WRITE(ioFile,'(2A,I6.6)') TRIM(vizdir),'/p',xyzcom_id
      IF(init) THEN
         OPEN(UNIT=ioUnit,FILE=TRIM(ioFile),FORM='UNFORMATTED',STATUS='REPLACE')
      ELSE
         OPEN(UNIT=ioUnit,FILE=TRIM(ioFile),FORM='UNFORMATTED',STATUS='OLD',POSITION='APPEND')
      END IF
      
      tmp = var
      WRITE(ioUnit) tmp

      CLOSE(ioUnit)
  END IF


END SUBROUTINE



SUBROUTINE write_MIR(vars,nvars,minvar,maxvar,nxx,nyy,nzz)
  USE mpi
  USE globals, ONLY: jobdir
  USE inputs, ONLY: nx,ny,nz,dx,dy,dz,x1,y1,z1,tviz,curvlin,fileorder
  USE post_proc_data
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: nvars,nxx,nyy,nzz
  CHARACTER(LEN=90), DIMENSION(nvars), INTENT(IN) :: vars
  DOUBLE PRECISION, DIMENSION(50), INTENT(IN)     :: minvar,maxvar
  CHARACTER(LEN=90) :: plotinfo
  CHARACTER(LEN=90) :: datafiles
  CHARACTER(LEN=90) :: gridfiles
  INTEGER :: i



  IF( xyzcom_id == 0) THEN
  WRITE(plotinfo,'(4A)') TRIM(jobdir),'/', TRIM(MIRname), '.mir'
  WRITE(datafiles,'(3A)') 'datafiles: ',TRIM(MIRname),'%04d/p%06d # data files'
  WRITE(gridfiles,'(3A)') 'gridfiles: ',TRIM(Gname),'/p%06d # data files'
  OPEN(UNIT=10,FILE=TRIM(plotinfo),FORM='FORMATTED',STATUS='REPLACE')
  WRITE(10,'(A)') 'VERSION 1.2'
  IF(curvlin) WRITE(10,'(A)') 'curvilinear: yes  # grid is curvilinear'		
  !!WRITE(10,'(A)') 'gridfiles: grid2d/p%06d     # grid files'
  WRITE(10,'(A)') TRIM(gridfiles)
  WRITE(10,'(A)') TRIM(datafiles)
  WRITE(10,'(3A)') 'fileorder: ',fileorder,' # processor order'
  WRITE(10,'(A,3I6,A)') 'domainsize: ',nxx,nyy,nzz,'  # overall grid dimensions'
  WRITE(10,'(A,3I6,A)') 'blocksize: ',nxx,nyy,nzz,'  # grid dim per proc'
  WRITE(10,'(A,3(X,D22.15),A)') 'origin:',x1,y1,z1,'  # xmin, ymin, zmin'
  WRITE(10,'(A,3(X,D22.15),A)') 'spacing:',dx,dy,dz,'  # dx, dy, dz'
  ! scalar and vector variables--------------------------------------------------------------
  WRITE(10,'(A,I4,A)')  'variables: ',nvars,'  # number of variables'
  i=1
  DO i=1,nvars
      WRITE(10,'(2X,2A,2(X,ES14.6))') TRIM(vars(i)), ' 1',minvar(i),maxvar(i)
  END DO 

  ! times--------------------------------------------------------------------------------------
  WRITE(10,'(A,I4,A)') 'timesteps: ',T_iter-Tstart+1,'  # number of times to plot'
  DO i=Tstart,T_iter
      WRITE(10,'(2X,I4.4,1(ES17.6))') i,REAL(i)*tviz ! dump numbers and times
  END DO
  CLOSE(10)
  END IF

END SUBROUTINE write_MIR




SUBROUTINE get_span(corr,nxs,nys,var)
  USE inputs, ONLY: nx,ny,nz
  USE globals, ONLY: flen,u,v,w,rho,p,ax,ay,az,ix,iy,iz
  USE constants, ONLY: zero,one,two,half
  USE metrics, ONLY: x_c,y_c,z_c
  USE interfaces, ONLY : SUMproc,subsum3xy,subsum3yz,subsum3xz,filter
  IMPLICIT NONE
  DOUBLE PRECISION, DIMENSION(nz), INTENT(OUT) :: corr
  INTEGER, INTENT(IN) :: nxs,nys
  CHARACTER(LEN=flen), INTENT(IN) :: var

  DOUBLE PRECISION, DIMENSION(ax,ay,az) :: cvar,tmp
  DOUBLE PRECISION, DIMENSION(ax) :: xline
  DOUBLE PRECISION, DIMENSION(ay) :: yline
  DOUBLE PRECISION, DIMENSION(nz) :: zdata
  DOUBLE PRECISION :: zmean
  INTEGER :: i,j,k

  !! Get the span-wise correlation at this specific xs,ys station.
  SELECT CASE(TRIM(var))
      CASE('u')
         cvar = u
      CASE('v')
         cvar = v
      CASE('w')
         cvar = w
      CASE('rho')
         cvar = rho
      CASE('p')
         cvar = p
      CASE('ptot')
         cvar = p + half*rho*sqrt( u*u + v*v + w*w )  !! Incompressible only
         
  END SELECT   

  xline = zero
  WHERE(ix == nxs) xline = one

  yline = zero
  WHERE(iy == nys) yline = one

  FORALL(i=1:ax,j=1:ay,k=1:az)
      tmp(i,j,k) = cvar(i,j,k)*xline(i)*yline(j)
  END FORALL

  zdata = SUBSUM3XY(tmp)  
  zmean = SUM(zdata)/DBLE(nz)
  corr = zdata !- zmean  ! Subtract out the mean



END SUBROUTINE




SUBROUTINE compute_AFLES()
  USE globals, ONLY: u,v,w,rho,p,T,e,c
  USE inputs, ONLY : radiative,magnetic,accpss,condsolver,dirB,fluidprops
  USE globals, ONLY : ax,ay,az,nRHS
  USE globals, ONLY : momentum_index,energy_index,mass_index,hydro_index,radiation_index
  USE RKscheme, ONLY : nRKsteps
  USE interfaces, ONLY : acceleration,stress,conduction,diffusion,momentum,energy,mass
  USE interfaces, ONLY : radenergy,induction,DAFLES,viscGrad,Dmomentum,Denergy
  USE interfaces, ONLY : timestep,eos,nscbc,boundary
  USE metrics, ONLY : bulkA,bulkB,bulkC
  IMPLICIT NONE
  DOUBLE PRECISION, DIMENSION(SIZE(rho,1),SIZE(rho,2),SIZE(rho,3)) :: Sxx,Syy,Szz,Sxy,Syz,Sxz             ! viscous stress tensor
  DOUBLE PRECISION, DIMENSION(SIZE(rho,1),SIZE(rho,2),SIZE(rho,3)) :: qx,qy,qz                            ! conductive heat flux



     CALL properties()         ! compute viscosity, thermal conductivity, magnetic diffusivity and species diffusivities
     CALL viscGrad(T,u,v,w,Sxx,Syy,Szz,Sxy,Sxz,Syz,qx,qy,qz)
     CALL DAFLES(rho,T,c,u,v,w,Sxx,Syy,Szz,Sxy,Sxz,Syz,qx,qy,qz)
     CALL properties()         ! compute viscosity, thermal conductivity, magnetic diffusivity and species diffusivities

END SUBROUTINE compute_AFLES



!! This sets up the variables needed for each processor to read and set the appropriate data.
!! This is based on nx,ny,nz and px,py,pz and Tpx,Tpy,Tpz and this proc's ID
SUBROUTINE setup_dataread()
  USE globals, ONLY: ix,iy,iz,ax,ay,az
  USE inputs, ONLY: nx,ny,nz
  USE post_proc_data
  IMPLICIT NONE

  INTEGER :: i,j,k,p,proc
  DOUBLE PRECISION :: tmp

  !! Set up the proc-map for the data
  Tax = nx / Tpx
  Tay = ny / Tpy
  Taz = nz / Tpz
  proc = Tpx*Tpy*Tpz
  ALLOCATE(procmap(proc,4))
  ALLOCATE(invprocmap(Tpx,Tpy,Tpz))
  
  !  Get the proc IJK layout
  p = 1
  DO k=1,Tpz
     DO j=1,Tpy
        DO i=1,Tpx
           procmap(p,1)=p-1
           procmap(p,2)=i-1
           procmap(p,3)=j-1
           procmap(p,4)=k-1
           invprocmap(i,j,k) = p-1
           p = p + 1
        END DO
     END DO
  END DO


  ! Get the x-direction indices in this chunk
  !! Set the processor bounds for the data
  i1g = ix(1)
  ifg = ix(ax)
  tmp = dble(i1g)/dble(Tax)
  i1p = CEILING(tmp) - 1
  tmp = dble(ifg)/dble(Tax)
  ifp = CEILING(tmp) - 1

  j1g = iy(1)
  jfg = iy(ay)
  tmp = dble(j1g)/dble(Tay)
  j1p = CEILING(tmp) - 1
  tmp = dble(jfg)/dble(Tay)
  jfp = CEILING(tmp) - 1

  k1g = iz(1)
  kfg = iz(az)
  tmp = dble(k1g)/dble(Taz)
  k1p = CEILING(tmp) - 1
  tmp = dble(kfg)/dble(Taz)
  kfp = CEILING(tmp) - 1

END SUBROUTINE setup_dataread


!! Each processor gets the data it needs for this particular time step/restart file (should handle both)
!! At the end of this routine iodata, should be fully set.
SUBROUTINE read_data(viz_num)
  USE mpi
  USE globals, ONLY: iodata,ix,iy,iz,ax,ay,az
  USE post_proc_data
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: viz_num

  INTEGER :: total,count,i,j,k,p,v
  INTEGER :: ig,jg,kg,iL,jL,kL
  INTEGER :: x1off,xnoff,y1off,ynoff,z1off,znoff
  CHARACTER(LEN=90) :: procdir,vizdir
  INTEGER :: punit=37
  DOUBLE PRECISION, DIMENSION(:,:,:,:), ALLOCATABLE :: Gdata,RESdata
  REAL, DIMENSION(:,:,:,:), ALLOCATABLE :: VISdata

  ! Get the offsets for this block of data
  x1off =  ix(1) - (i1p*Tax + 1) + 1
  xnoff = (ifp*Tax + 1) - ix(ax) + 1
  y1off =  iy(1) - (j1p*Tay + 1) + 1
  ynoff = (jfp*Tay + 1) - iy(ay) + 1
  z1off =  iz(1) - (k1p*Taz + 1) + 1
  znoff = (kfp*Taz + 1) - iz(az) + 1

  SELECT CASE(ftype)
      CASE('res')
         ALLOCATE(RESdata(Tax,Tay,Taz,varDIM))
         WRITE(vizdir,'(2A,I4.4)') TRIM(restart_dir),'/res',viz_num
      CASE('vis')
         ALLOCATE(VISdata(Tax,Tay,Taz,varDIM))
         WRITE(vizdir,'(2A,I4.4)') TRIM(restart_dir),'/vis',viz_num
  END SELECT

  ALLOCATE(Gdata( (ifp-i1p+1)*Tax, (jfp-j1p+1)*Tay, (kfp-k1p+1)*Taz, varDIM))


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

           END SELECT

        END DO
     END DO
  END DO

  ! Set to the global iodata variable pointer array
  iodata = Gdata(x1off:xnoff,y1off:ynoff,z1off:znoff,:)


  ! Free the memory
  SELECT CASE(ftype)
      CASE('res')
         DEALLOCATE(Gdata,RESdata)
      CASE('vis')
         DEALLOCATE(Gdata,VISdata)
  END SELECT


END SUBROUTINE read_data



