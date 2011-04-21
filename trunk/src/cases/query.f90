!  Britton J. Olson
!  Stanford University
!  March. 4th, 2011

!  Description: Compute the mean velocity profiles for a turbulent boundary layer
!       and make the appropriate scaling in y+ and u_tau units
MODULE query_data
  SAVE
  
  INTEGER :: ix1 = 1        ! X-index 
  INTEGER :: iy1 = 1        ! Y-index 
  INTEGER :: iz1 = 1        ! Z-index 

END MODULE query_data


PROGRAM query
  USE operators
  USE globals
  USE query_data
  IMPLICIT NONE
  INTEGER :: i

  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: output

  REAL, DIMENSION(DIM) :: Rout
  CHARACTER(LEN=flen) :: ofile
  CHARACTER(LEN=flen) :: vfile
  CHARACTER(LEN=flen) :: inputFile

  INTEGER :: nviz,iviz,funit=34
  NAMELIST /query_vars/ ix1,iy1,iz1


  !! General input file for t1,tf,nx,ny,nz,filename
  CALL read_input

  ! Read in namelist for specific case stuff  
  CALL GETARG(1,inputFile)
  OPEN(UNIT=funit,FILE=TRIM(inputFile),FORM='FORMATTED',STATUS='OLD')
  READ(UNIT=funit,NML=query_vars)
  CLOSE(funit)


  nviz=tf-t1+1

  !! Initialize the data arrays
  ALLOCATE(output( DIM) )



  !! Main loop over the viz-files
  DO i=1,nviz
     iviz = i + t1 - 1
     output = 0.0D0
     
     ! Get the next viz-directory and call kernal
     CALL viz_name(jobdir,iviz,vfile)
     CALL POINT(output,vfile,ix1,iy1,iz1)
     Rout = output
     PRINT*,'Query time:',iviz

     IF (i==1) THEN
        OPEN(UNIT=33,FILE='query.dat',FORM='FORMATTED',STATUS='UNKNOWN')
        WRITE(33,*) '# X, Y, Z:',Rout(x_c),Rout(y_c),Rout(z_c) 
        WRITE(33,*) '# 1-4    : step,U,V,W'
        WRITE(33,*) '# 5-7    : P,T,rho'
        WRITE(33,*) '# 8-10   : mu,beta,ktc'
     END IF

     WRITE(33,'(10ES12.4)') REAL(iviz),Rout(U),Rout(V),Rout(W), &
          & Rout(P),Rout(T),Rout(rho),Rout(mu),Rout(bulk),Rout(ktc)

  END DO


  PRINT*,'Point query: (I,J,K) = (',ix1,',',iy1,',',iz1,')'
  PRINT*,'X, Y, Z    :',Rout(x_c),Rout(y_c),Rout(z_c)
  PRINT*,'U, V, W    :',Rout(U),Rout(V),Rout(W)
  PRINT*,'P, T, rho  :',Rout(P),Rout(T),Rout(rho)
  PRINT*,'Mu,Beta,Ktc:',Rout(mu),Rout(bulk),Rout(ktc)
  
  DEALLOCATE(output)


END PROGRAM query

