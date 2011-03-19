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


PROGRAM post
  USE operators
  USE globals, ONLY: flen,varDIM,ax,ay,az,nx,ny,nz,t1,jobdir
  !USE sumdata, ONLY: out,gout
  IMPLICIT NONE
  INTEGER :: i
  INTEGER, DIMENSION(2,2) :: stride
  INTEGER, DIMENSION(:), ALLOCATABLE :: vars

  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: output
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: Uvd
  CHARACTER(LEN=flen) :: ofile
  CHARACTER(LEN=flen) :: vfile

  DOUBLE PRECISION :: dudy,tauw,mu,rho_w,rho_0,del,utau,off,dUp


  CALL read_input

 
  stride(1,1) = 100
  stride(1,2) = 150
  stride(2,1) = 1
  stride(2,2) = 64
  ALLOCATE(vars(9))
  vars = (/ 1,2,3,4,5,6,7,8,9 /)
  ALLOCATE(output( ny , SIZE(vars)+3) )

  CALL viz_name(jobdir,941,vfile)
  CALL SUBSUM3IK(stride,output,vfile)

  !! Dump the unscaled values
  WRITE(ofile,'(A,I4.4)') 'Uwall.',t1
  OPEN(UNIT=33,FILE=ofile,FORM='FORMATTED',STATUS='UNKNOWN')
  WRITE(33,*) '%% Y (cm), U(cm/s), rho(g/cm^3),mu( g/ s*cm )'
  WRITE(33,*) '%% Ny = 256'
  DO i=1,ny
     WRITE(33,'(4ES12.4)') output(i,2),output(i,4),output(i,7),output(i,12)
  END DO



!!$  !! Get some y_plus units
!!$  rho_w = output(1,7)
!!$  rho_0 = output(ny/2,7)
!!$  mu = output(ny/2,12)
!!$  !mu = output(1,12)
!!$  !mu = mu * 2.0d0
!!$  print*,mu
!!$  !mu = rho_0 * out(ny/2,1) * 0.1D0 / 50.0D3
!!$  
!!$
!!$  dudy = ( output(2,4) - output(1,4) )/ ( output(2,2) - output(1,2) )
!!$  tauw = mu * dudy
!!$  utau = sqrt( tauw / rho_w )
!!$  print*,'U_tau',utau
!!$
!!$  !! Scale y
!!$  del = mu / ( rho_w * utau )
!!$  output(:,2) = output(:,2) / del
!!$  off = output(1,2)
!!$  output(:,2) = output(:,2) - off
!!$  print*,'Del',del
!!$
!!$  
!!$
!!$
!!$  !! Scale velocity
!!$  output(:,4) = output(:,4) / utau
!!$
!!$
!!$  !! Van Driest transformation
!!$  ALLOCATE(Uvd(ny))
!!$  Uvd(1) = 0.0D0
!!$  DO i=2,ny/2
!!$     dUp = output(i,4)-output(i-1,4)
!!$     Uvd(i) = Uvd(i-1) + sqrt(output(i,7)/rho_w) * dUp
!!$  END DO
!!$
!!$  Uvd(ny) = 0.0D0
!!$  DO i=ny-1,ny/2+1,-1
!!$     dUp = output(i,4)-output(i+1,4)
!!$     Uvd(i) = Uvd(i+1) + sqrt(output(i,7)/rho_w) * dUp
!!$  END DO
!!$
!!$
!!$
!!$  OPEN(UNIT=33,FILE='Uwall.dat',FORM='FORMATTED',STATUS='UNKNOWN')
!!$  DO i=1,ny/2
!!$     WRITE(33,'(4ES12.4)') output(i,2),output(i,4),output(i,7),output(i,12)
!!$     !WRITE(33,*)  output(i,2),Uvd(i)     ,output(i,4)
!!$  END DO
!!$  DO i=ny/2+1,ny
!!$     WRITE(33,'(4ES12.4)') (output(ny,2)-output(i,2)),output(i,4),output(i,7),output(i,12)
!!$     !WRITE(33,*) (output(ny,2)-output(i,2)),Uvd(i)      ,output(i,4)
!!$  END DO
!!$
!!$
!!$  !! Dump the unscaled values
!!$
!!$
!!$  CLOSE(33)


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



!!$SUBROUTINE viz_loop()
!!$  USE globals
!!$  USE sumdata, ONLY: out,gout,norm
!!$  IMPLICIT NONE
!!$  INTEGER :: i,tviz
!!$  CHARACTER(LEN=flen) :: vizdir
!!$
!!$  out = zero
!!$  gout = zero
!!$
!!$  DO i=1,tloop
!!$     tviz = i + t1 - 1
!!$     WRITE(vizdir,'(2A,I4.4)') TRIM(jobdir),'/vis',tviz
!!$     CALL viz_ops(vizdir)
!!$     CALL proc_loop(vizdir)
!!$  END DO
!!$
!!$
!!$  !! Normalize the sums here over space
!!$  out = out / norm
!!$  gout = gout / norm
!!$
!!$  !! Normalize the sums here over space
!!$  out = out / dble(tf-t1+1)
!!$  gout = gout / dble(tf-t1+1)
!!$  
!!$
!!$
!!$END SUBROUTINE viz_loop





!!$SUBROUTINE PLANE_AVE(direction,vars,plane2D)
!!$  USE globals
!!$  USE sumdata
!!$  IMPLICIT NONE
!!$  CHARACTER(LEN=10) :: direction
!!$  INTEGER, DIMENSION(:) :: vars
!!$  DOUBLE PRECISION, DIMENSION(:,:,:) :: plane2D
!!$
!!$
!!$  INTEGER, DIMENSION(SIZE(vars)) :: voff
!!$  DOUBLE PRECISION :: tmp
!!$
!!$
!!$  PLNAVE = .TRUE.
!!$  ALLOCATE(pout(nx,ny,varDIM))
!!$  ALLOCATE(pgrid(nx,ny,3))
!!$
!!$  ! Loop over ALL procs
!!$  i1p = 0
!!$  ifp = px-1
!!$  j1p = 0       !j1g/ay
!!$  jfp = py-1    !jfg/ay
!!$  k1p = 0
!!$  kfp = pz-1
!!$  
!!$  norm = dble(nz) 
!!$
!!$  CALL viz_loop()
!!$
!!$
!!$  plane2D(:,:,1:3) = pgrid
!!$  prof1D(:,4:size(vars)+3) = pout(:,vars)
!!$
!!$  CALL clean_sum
!!$
!!$END SUBROUTINE PLANE_AVE

