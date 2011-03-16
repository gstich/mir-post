!  Britton J. Olson
!  Stanford University
!  March. 4th, 2011

!  Subroutines and framework for getting at a Miranda data set outside of visit
!  As of now, this will be a serial operation but can/should be extended to parallel later.


!  Routines to use that are identical in Miranda source
!  SUBSUM3IJ - Sum Kplane and return scalar, F(k), (options: stride(1:4) = [i1,in,j1,jn], sets width in i and j of IJplane)
!  SUBSUM3JK - Sum Iplane and return scalar, F(i), (options: stride(1:4) = [j1,jn,k1,kn], sets width in j and k of JKplane)
!  SUBSUM3IK - Sum Jplane and return scalar, F(j), (options: stride(1:4) = [i1,in,k1,kn], sets width in i and k of IKplane)

!  SUBSUMI - Makes a JK plane of data at given global I
!  SUBSUMJ - Makes a IK plane of data at given global J
!  SUBSUMK - Makes a IJ plane of data at given global K

!  SUMtime - Sum quantity over all time

!  PROC_level - Loop over processor files for a given vis step 
!  VIS_level  - Loop over all time files




!  Program will get 2-d and 1-d planes and lines based on global i,j,k coordinates
!  given by the user.  This is helpful for boundary layer grids

!  Also, 
MODULE post_routines


INTERFACE SUBSUM3IK
   SUBROUTINE SUBSUM3IK(stride,vars,prof1D)
     IMPLICIT NONE
     DOUBLE PRECISION, DIMENSION(:,:) :: prof1D
     INTEGER, DIMENSION(:) :: vars
     INTEGER, DIMENSION(2,2) :: stride
   END SUBROUTINE SUBSUM3IK
END INTERFACE


INTERFACE SUBSUM3IJ
   SUBROUTINE SUBSUM3IJ(stride,vars,prof1D)
     IMPLICIT NONE
     DOUBLE PRECISION, DIMENSION(:,:) :: prof1D
     INTEGER, DIMENSION(:) :: vars
     INTEGER, DIMENSION(2,2) :: stride
   END SUBROUTINE SUBSUM3IJ
END INTERFACE


INTERFACE SUBSUM3JK
   SUBROUTINE SUBSUM3JK(stride,vars,prof1D)
     IMPLICIT NONE
     DOUBLE PRECISION, DIMENSION(:,:) :: prof1D
     INTEGER, DIMENSION(:) :: vars
     INTEGER, DIMENSION(2,2) :: stride
   END SUBROUTINE SUBSUM3JK
END INTERFACE


INTERFACE PLAN_AVE
   SUBROUTINE PLAN_AVE(direction,vars,plane2D)
     IMPLICIT NONE
     CHARACTER(LEN=10) :: direction
     INTEGER, DIMENSION(:) :: vars
     DOUBLE PRECISION, DIMENSION(:,:,:) :: plane2D
   END SUBROUTINE PLAN_AVE
END INTERFACE



END MODULE post_routines

MODULE globals
  IMPLICIT NONE
  INTEGER, PARAMETER :: flen=90          ! file length

  INTEGER :: nx,ny,nz
  INTEGER :: t1,tf,tloop
  CHARACTER(LEN=flen) :: jobdir

  INTEGER :: ax,ay,az
  INTEGER :: px,py,pz
  INTEGER :: varDIM
  INTEGER :: proc
  INTEGER, DIMENSION(:,:), ALLOCATABLE :: procmap
  INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: invprocmap
  
  DOUBLE PRECISION, PARAMETER :: pi=3.14159265358979323846d0,zero=0.0d0,one=1.0d0,two=2.0d0

END MODULE globals

MODULE sumdata
  IMPLICIT NONE

  INTEGER :: i1p,ifp,j1p,jfp,k1p,kfp
  INTEGER :: i1g,ifg,j1g,jfg,k1g,kfg
  LOGICAL :: fSS3IK = .FALSE.
  LOGICAL :: fSS3JK = .FALSE.
  LOGICAL :: fSS3IJ = .FALSE.
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: out
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: gout
  DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE :: pout,pgrid

  DOUBLE PRECISION :: norm


END MODULE sumdata



PROGRAM post
  USE post_routines 
  USE globals, ONLY: varDIM,ax,ay,az,nx,ny,nz,t1
  !USE sumdata, ONLY: out,gout
  IMPLICIT NONE
  INTEGER :: i
  INTEGER, DIMENSION(2,2) :: stride
  INTEGER, DIMENSION(:), ALLOCATABLE :: vars

  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: output
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: Uvd
  CHARACTER(LEN=50) :: ofile

  DOUBLE PRECISION :: dudy,tauw,mu,rho_w,rho_0,del,utau,off,dUp


  CALL read_input
  varDIM = 12

 
  stride(1,1) = 1
  stride(1,2) = 1
  stride(2,1) = 1
  stride(2,2) = 1
  ALLOCATE(vars(9))
  vars = (/ 1,2,3,4,5,6,7,8,9 /)
  ALLOCATE(output( ny , SIZE(vars)+3) )
  CALL SUBSUM3IK(stride,vars,output)

  !! Dump the unscaled values
  WRITE(ofile,'(A,I4.4)') 'Uwall.',t1
  OPEN(UNIT=33,FILE=ofile,FORM='FORMATTED',STATUS='UNKNOWN')
  WRITE(33,*) '%% Y (cm), U(cm/s), rho(g/cm^3),mu( g/ s*cm )'
  WRITE(33,*) '%% Ny = 256'
  DO i=1,ny
     WRITE(33,'(4ES12.4)') output(i,2),output(i,4),output(i,7),output(i,12)
  END DO



  !! Get some y_plus units
  rho_w = output(1,7)
  rho_0 = output(ny/2,7)
  mu = output(ny/2,12)
  !mu = output(1,12)
  !mu = mu * 2.0d0
  print*,mu
  !mu = rho_0 * out(ny/2,1) * 0.1D0 / 50.0D3
  

  dudy = ( output(2,4) - output(1,4) )/ ( output(2,2) - output(1,2) )
  tauw = mu * dudy
  utau = sqrt( tauw / rho_w )
  print*,'U_tau',utau

  !! Scale y
  del = mu / ( rho_w * utau )
  output(:,2) = output(:,2) / del
  off = output(1,2)
  output(:,2) = output(:,2) - off
  print*,'Del',del

  


  !! Scale velocity
  output(:,4) = output(:,4) / utau


  !! Van Driest transformation
  ALLOCATE(Uvd(ny))
  Uvd(1) = 0.0D0
  DO i=2,ny/2
     dUp = output(i,4)-output(i-1,4)
     Uvd(i) = Uvd(i-1) + sqrt(output(i,7)/rho_w) * dUp
  END DO

  Uvd(ny) = 0.0D0
  DO i=ny-1,ny/2+1,-1
     dUp = output(i,4)-output(i+1,4)
     Uvd(i) = Uvd(i+1) + sqrt(output(i,7)/rho_w) * dUp
  END DO



  OPEN(UNIT=33,FILE='Uwall.dat',FORM='FORMATTED',STATUS='UNKNOWN')
  DO i=1,ny/2
     WRITE(33,'(4ES12.4)') output(i,2),output(i,4),output(i,7),output(i,12)
     !WRITE(33,*)  output(i,2),Uvd(i)     ,output(i,4)
  END DO
  DO i=ny/2+1,ny
     WRITE(33,'(4ES12.4)') (output(ny,2)-output(i,2)),output(i,4),output(i,7),output(i,12)
     !WRITE(33,*) (output(ny,2)-output(i,2)),Uvd(i)      ,output(i,4)
  END DO


  !! Dump the unscaled values


  CLOSE(33)


!!$  DEALLOCATE(vars)
!!$  DEALLOCATE(output)
!!$  stride(1,1) = 120
!!$  stride(1,2) = 140
!!$  stride(2,1) = 1
!!$  stride(2,2) = 32
!!$  ALLOCATE(vars(3))
!!$  vars = (/ 1 , 2 , 3 /)
!!$  ALLOCATE(output( nx , SIZE(vars)+3) )
!!$  CALL SUBSUM3JK(stride,vars,output)
!!$
!!$  OPEN(UNIT=33,FILE='Unoz.dat',FORM='FORMATTED',STATUS='UNKNOWN')
!!$  DO i=1,nx
!!$     WRITE(33,*) output(i,1),output(i,4)
!!$  END DO
!!$  CLOSE(33)
!!$
!!$
!!$  DEALLOCATE(vars)
!!$  DEALLOCATE(output)
!!$  stride(1,1) = 200
!!$  stride(1,2) = 225
!!$  stride(2,1) = 3
!!$  stride(2,2) = 20
!!$  ALLOCATE(vars(3))
!!$  vars = (/ 1 , 2 , 3 /)
!!$  ALLOCATE(output( nz , SIZE(vars)+3) )
!!$  CALL SUBSUM3IJ(stride,vars,output)
!!$
!!$  OPEN(UNIT=33,FILE='U_z.dat',FORM='FORMATTED',STATUS='UNKNOWN')
!!$  DO i=1,nz
!!$     WRITE(33,*) output(i,3),output(i,4)
!!$  END DO
!!$  CLOSE(33)

  !! Main loop over viz_dumpfiles
  !CALL viz_loop()





  !print*,'y+=',gout(2,2),' x+=',(gout(2,1)-gout(1,1))/del,' z+=',(gout(2,3)-gout(1,3))/del


END PROGRAM post


SUBROUTINE viz_ops(vizdir)
  USE globals, ONLY: flen
  IMPLICIT NONE
  CHARACTER(LEN=flen),INTENT(IN) :: vizdir

  PRINT*,'Reading from',TRIM(vizdir)


END SUBROUTINE viz_ops



SUBROUTINE viz_loop()
  USE globals
  USE sumdata, ONLY: out,gout,norm
  IMPLICIT NONE
  INTEGER :: i,tviz
  CHARACTER(LEN=flen) :: vizdir

  out = zero
  gout = zero

  DO i=1,tloop
     tviz = i + t1 - 1
     WRITE(vizdir,'(2A,I4.4)') TRIM(jobdir),'/vis',tviz
     CALL viz_ops(vizdir)
     CALL proc_loop(vizdir)
  END DO


  !! Normalize the sums here over space
  out = out / norm
  gout = gout / norm

  !! Normalize the sums here over space
  out = out / dble(tf-t1+1)
  gout = gout / dble(tf-t1+1)
  


END SUBROUTINE viz_loop


SUBROUTINE proc_loop(vizdir)
  USE globals
  USE sumdata, ONLY: i1p,ifp,j1p,jfp,k1p,kfp,out,gout,norm
  IMPLICIT NONE
  CHARACTER(LEN=flen), INTENT(IN) :: vizdir
  REAL(KIND=4), DIMENSION(ax,ay,az,varDIM) :: iodata
  REAL(KIND=4), DIMENSION(ax,ay,az,3) :: grid
  CHARACTER(LEN=flen) :: procdir
  INTEGER :: i,j,k,p,v
  INTEGER :: punit=23
  INTEGER :: ig,jg,kg
  

  !ALLOCATE(iodata(ax,ay,ax,varDIM))

  ! Loop only over PROCS in pertainate range
  DO i=i1p,ifp
     DO j=j1p,jfp
        DO k=k1p,kfp
           
           ! Global indices on corner of proc grid
           ig = i*ax + 1
           jg = j*ay + 1
           kg = k*az + 1
           
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
           WRITE(procdir,'(2A,I6.6)') TRIM(jobdir),'/grid/p',p
           OPEN(UNIT=punit,FILE=TRIM(procdir),FORM='UNFORMATTED',STATUS='OLD')
           DO v=1,3
              READ(punit) grid(:,:,:,v)
           END DO
           CLOSE(punit)
           

           CALL SUMSUM(iodata,grid,ig,jg,kg)
           

        END DO
     END DO
  END DO



END SUBROUTINE proc_loop


SUBROUTINE SUMSUM(iodata,grid,ig,jg,kg)
  USE globals
  USE sumdata
  IMPLICIT NONE
  REAL(KIND=4), DIMENSION(ax,ay,az,varDIM), INTENT(IN) :: iodata
  REAL(KIND=4), DIMENSION(ax,ay,az,3), INTENT(IN) :: grid
  INTEGER, INTENT(IN) :: ig,jg,kg
  REAL(KIND=4), DIMENSION(varDIM) :: temp
  REAL(KIND=4), DIMENSION(3) :: gtemp
  INTEGER :: i,j,k,v
  INTEGER :: irngL,irngU,jrngL,jrngU,krngL,krngU


  temp = zero
  gtemp = zero

  IF(fSS3IK) THEN
     irngL = max(i1g-ig+1,1)
     irngU = min(ifg-ig+1,ax)

     krngL = max(k1g-kg+1,1)
     krngU = min(kfg-kg+1,az)
     
     print*,irngL,irngU,krngL,krngU

     DO j = 1,ay
        DO v = 1,varDIM
           temp(v) = SUM(iodata( irngL:irngU ,j , krngL:krngU , v) )
        END DO
        out(jg+j-1,:) = out(jg+j-1,:) + temp
        DO v = 1,3
           gtemp(v) = SUM(grid( irngL:irngU ,j , krngL:krngU , v) )
        END DO
        gout(jg+j-1,:) = gout(jg+j-1,:) + gtemp
     END DO
  END IF


  IF(fSS3JK) THEN
     jrngL = max(j1g-jg+1,1)
     jrngU = min(jfg-jg+1,ay)

     krngL = max(k1g-kg+1,1)
     krngU = min(kfg-kg+1,az)
     
     print*,jrngL,jrngU,krngL,krngU

     DO i = 1,ax
        DO v = 1,varDIM
           temp(v) = SUM(iodata( i, jrngL:jrngU , krngL:krngU , v) )
        END DO
        out(ig+i-1,:) = out(ig+i-1,:) + temp
        DO v = 1,3
           gtemp(v) = SUM(grid( i, jrngL:jrngU , krngL:krngU , v) )
        END DO
        gout(ig+i-1,:) = gout(ig+i-1,:) + gtemp
     END DO
  END IF


  IF(fSS3IJ) THEN
     irngL = max(i1g-ig+1,1)
     irngU = min(ifg-ig+1,ax)

     jrngL = max(j1g-jg+1,1)
     jrngU = min(jfg-jg+1,ay)
     
     print*,irngL,irngU,jrngL,jrngU

     DO k = 1,az
        DO v = 1,varDIM
           temp(v) = SUM(iodata( irngL:irngU , jrngL:jrngU , k , v) )
        END DO
        out(kg+k-1,:) = out(kg+k-1,:) + temp
        DO v = 1,3
           gtemp(v) = SUM(grid( irngL:irngU , jrngL:jrngU ,k , v) )
        END DO
        gout(kg+k-1,:) = gout(kg+k-1,:) + gtemp
     END DO
  END IF

  IF(PLNAVE) THEN


     !! Reduce 3d data to a plane and add to global 2d array
     IF(z) THEN
        DO i=1,ax
           DO j=1,ay
              iig = ig+i-1
              jjg = jg+j-1
              DO v=1,varDIM
                 pout(iig,jjg,:) = pout(iig,jjg,:) + SUM(iodata(i,j,:,v))
              END DO
              DO v=1,3
                 pgrid(iig,jjg,:) = pgrid(iig,jjg,:) + SUM(grid(i,j,:,v))
              END DO
           END DO
        END DO
     END IF




  END IF

END SUBROUTINE SUMSUM

SUBROUTINE SUBSUM3IK(stride,vars,prof1D)
  USE globals
  USE sumdata
  IMPLICIT NONE
  DOUBLE PRECISION, DIMENSION(:,:) :: prof1D
  INTEGER, DIMENSION(:) :: vars
  INTEGER, DIMENSION(2,2) :: stride

  INTEGER, DIMENSION(SIZE(vars)) :: voff
  DOUBLE PRECISION :: tmp

  i1g = stride(1,1)
  ifg = stride(1,2)
  
  k1g = stride(2,1)
  kfg = stride(2,2)

  fSS3IK = .TRUE.
  ALLOCATE(out(ny,varDIM))
  ALLOCATE(gout(ny,3))
  
  tmp = dble(i1g)/dble(ax)
  i1p = CEILING(tmp) - 1
  tmp = dble(ifg)/dble(ax)
  ifp = CEILING(tmp) - 1

  j1p = 0       !j1g/ay
  jfp = py-1    !jfg/ay

  tmp = dble(k1g)/dble(az)
  k1p = CEILING(tmp) - 1
  tmp = dble(kfg)/dble(az)
  kfp = CEILING(tmp) - 1

  CALL print_map
   
  norm = dble( ifg - i1g + 1) * dble( kfg - k1g + 1)

  CALL viz_loop()


  prof1D(:,1:3) = gout
  voff = vars + 3
  prof1D(:,voff) = out(:,vars)
  

  CALL clean_sum

END SUBROUTINE SUBSUM3IK


SUBROUTINE SUBSUM3JK(stride,vars,prof1D)
  USE globals
  USE sumdata
  IMPLICIT NONE
  DOUBLE PRECISION, DIMENSION(:,:) :: prof1D
  INTEGER, DIMENSION(:) :: vars
  INTEGER, DIMENSION(2,2), INTENT(IN) :: stride

  INTEGER, DIMENSION(SIZE(vars)) :: voff
  DOUBLE PRECISION :: tmp

  j1g = stride(1,1)
  jfg = stride(1,2)
  
  k1g = stride(2,1)
  kfg = stride(2,2)

  fSS3JK = .TRUE.
  ALLOCATE(out(nx,varDIM))
  ALLOCATE(gout(nx,3))
  
  i1p = 0       
  ifp = px-1    

  tmp = dble(j1g)/dble(ay)
  j1p = CEILING(tmp) - 1
  tmp = dble(jfg)/dble(ay)
  jfp = CEILING(tmp) - 1

  tmp = dble(k1g)/dble(az)
  k1p = CEILING(tmp) - 1
  tmp = dble(kfg)/dble(az)
  kfp = CEILING(tmp) - 1

  CALL print_map
    
  
  norm = dble( jfg - j1g + 1) * dble( kfg - k1g + 1)

  CALL viz_loop()


  prof1D(:,1:3) = gout
  voff = vars + 3
  prof1D(:,voff) = out(:,vars)
  

  CALL clean_sum

END SUBROUTINE SUBSUM3JK


SUBROUTINE SUBSUM3IJ(stride,vars,prof1D)
  USE globals
  USE sumdata
  IMPLICIT NONE
  DOUBLE PRECISION, DIMENSION(:,:) :: prof1D
  INTEGER, DIMENSION(2,2), INTENT(IN) :: stride
  INTEGER, DIMENSION(:) :: vars

  INTEGER, DIMENSION(SIZE(vars)) :: voff
  DOUBLE PRECISION :: tmp

  i1g = stride(1,1)
  ifg = stride(1,2)
  
  j1g = stride(2,1)
  jfg = stride(2,2)

  fSS3IJ = .TRUE.
  ALLOCATE(out(nz,varDIM))
  ALLOCATE(gout(nz,3))
  
  tmp = dble(i1g)/dble(ax)
  i1p = CEILING(tmp) - 1
  tmp = dble(ifg)/dble(ax)
  ifp = CEILING(tmp) - 1

  tmp = dble(j1g)/dble(az)
  j1p = CEILING(tmp) - 1
  tmp = dble(jfg)/dble(az)
  jfp = CEILING(tmp) - 1

  k1p = 0       !j1g/ay
  kfp = pz-1    !jfg/ay

  CALL print_map

  norm = dble( ifg - i1g + 1) * dble( jfg - j1g + 1)

  CALL viz_loop()


  prof1D(:,1:3) = gout
  voff = vars + 3
  prof1D(:,voff) = out(:,vars)
  

  CALL clean_sum

END SUBROUTINE SUBSUM3IJ


SUBROUTINE PLANE_AVE(direction,vars,plane2D)
  USE globals
  USE sumdata
  IMPLICIT NONE
  CHARACTER(LEN=10) :: direction
  INTEGER, DIMENSION(:) :: vars
  DOUBLE PRECISION, DIMENSION(:,:,:) :: plane2D


  INTEGER, DIMENSION(SIZE(vars)) :: voff
  DOUBLE PRECISION :: tmp


  PLNAVE = .TRUE.
  ALLOCATE(pout(nx,ny,varDIM))
  ALLOCATE(pgrid(nx,ny,3))

  ! Loop over ALL procs
  i1p = 0
  ifp = px-1
  j1p = 0       !j1g/ay
  jfp = py-1    !jfg/ay
  k1p = 0
  kfp = pz-1
  
  norm = dble(nz) 

  CALL viz_loop()


  plane2D(:,:,1:3) = pgrid
  prof1D(:,4:size(vars)+3) = pout(:,vars)

  CALL clean_sum

END SUBROUTINE PLANE_AVE

SUBROUTINE clean_sum
  USE sumdata
  IMPLICIT NONE
  IF(ALLOCATED(out)) DEALLOCATE(out)
  IF(ALLOCATED(gout)) DEALLOCATE(gout)
  fSS3IK = .FALSE.
  fSS3JK = .FALSE.
  fSS3IJ = .FALSE.
END SUBROUTINE clean_sum




SUBROUTINE read_input
  USE globals
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
  READ(ioUnit,*) comments
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

  

  tloop = tf-t1+1


END SUBROUTINE read_input



SUBROUTINE print_map
  USE sumdata
  IMPLICIT NONE

  PRINT*,"Global index ranges"
  PRINT*,'Iproc: (',i1g,',',ifg,')'
  PRINT*,'Jproc: (',j1g,',',jfg,')'
  PRINT*,'Kproc: (',k1g,',',kfg,')'

  PRINT*,"Proc index ranges"
  PRINT*,'Iproc: (',i1p,',',ifp,')'
  PRINT*,'Jproc: (',j1p,',',jfp,')'
  PRINT*,'Kproc: (',k1p,',',kfp,')'

END SUBROUTINE print_map

