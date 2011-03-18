


SUBROUTINE block_ops(iodata,corner,stride,Gdata,FUNK)
  USE globals, ONLY: varDIM,coordDIM,DIM,flen,zero,one,ax,ay,az
  IMPLICIT NONE
  REAL(KIND=4), DIMENSION(:,:,:,:), INTENT(IN) :: iodata
  INTEGER, DIMENSION(3), INTENT(IN) :: corner
  INTEGER, DIMENSION(2,2), INTENT(IN) :: stride
  DOUBLE PRECISION, DIMENSION(:,:,:,:), INTENT(INOUT) :: Gdata
  CHARACTER(LEN=flen), INTENT(IN) :: FUNK

  REAL(KIND=4), DIMENSION(varDIM) :: temp
  REAL(KIND=4), DIMENSION(coordDIM) :: gtemp
  DOUBLE PRECISION :: ddum
  INTEGER :: i,j,k,v
  INTEGER :: irngL,irngU,jrngL,jrngU,krngL,krngU
  INTEGER :: ig,jg,kg
  INTEGER :: i1g,ifg,j1g,jfg,k1g,kfg

  temp = zero
  gtemp = zero

  ig=corner(1);jg=corner(2);kg=corner(3);


!  IF(fSS3IK) THEN
  SELECT CASE(FUNK)
  

  CASE('SUBSUM3IK')

     !! Set the index box size for the averaging kernal
     i1g = stride(1,1)
     ifg = stride(1,2)
     k1g = stride(2,1)
     kfg = stride(2,2)
     
     irngL = max(i1g-ig+1,1)
     irngU = min(ifg-ig+1,ax)
     krngL = max(k1g-kg+1,1)
     krngU = min(kfg-kg+1,az)

     DO j = 1,ay
        DO v = 1,DIM
           ddum = SUM(iodata( irngL:irngU ,j , krngL:krngU , v) )
           Gdata(1,jg+j-1,1,v) = Gdata(1,jg+j-1,1,v) + ddum
        END DO
     END DO


  CASE('SUBSUM3JK')

     !! Set the index box size for the averaging kernal
     j1g = stride(1,1)
     jfg = stride(1,2)
     k1g = stride(2,1)
     kfg = stride(2,2)

     jrngL = max(j1g-jg+1,1)
     jrngU = min(jfg-jg+1,ay)
     krngL = max(k1g-kg+1,1)
     krngU = min(kfg-kg+1,az)

     DO i = 1,ax
        DO v = 1,DIM
           ddum = SUM(iodata( i, jrngL:jrngU , krngL:krngU , v) )
           Gdata(ig+i-1,1,1,v) = Gdata(ig+i-1,1,1,v) + ddum
        END DO
     END DO


  CASE('SUBSUM3IJ')
     
     !! Set the index box size for the averaging kernal
     i1g = stride(1,1)
     ifg = stride(1,2)
     j1g = stride(2,1)
     jfg = stride(2,2)
     
     irngL = max(i1g-ig+1,1)
     irngU = min(ifg-ig+1,ax)

     jrngL = max(j1g-jg+1,1)
     jrngU = min(jfg-jg+1,ay)
     
     DO k = 1,az
        DO v = 1,DIM
           ddum = SUM(iodata( irngL:irngU , jrngL:jrngU , k , v) )
           Gdata(1,1,kg+k-1,v) = Gdata(1,1,kg+k-1,v) + ddum
        END DO
     END DO

  END SELECT




!!$
!!$  IF(fSS3IJ) THEN
!!$     irngL = max(i1g-ig+1,1)
!!$     irngU = min(ifg-ig+1,ax)
!!$
!!$     jrngL = max(j1g-jg+1,1)
!!$     jrngU = min(jfg-jg+1,ay)
!!$     
!!$     print*,irngL,irngU,jrngL,jrngU
!!$
!!$     DO k = 1,az
!!$        DO v = 1,varDIM
!!$           temp(v) = SUM(iodata( irngL:irngU , jrngL:jrngU , k , v) )
!!$        END DO
!!$        out(kg+k-1,:) = out(kg+k-1,:) + temp
!!$        DO v = 1,3
!!$           gtemp(v) = SUM(grid( irngL:irngU , jrngL:jrngU ,k , v) )
!!$        END DO
!!$        gout(kg+k-1,:) = gout(kg+k-1,:) + gtemp
!!$     END DO
!!$  END IF

!  IF(PLNAVE) THEN


     !! Reduce 3d data to a plane and add to global 2d array
!     IF(z) THEN
!        DO i=1,ax
!           DO j=1,ay
!              iig = ig+i-1
!              jjg = jg+j-1
!              DO v=1,varDIM
!                 pout(iig,jjg,:) = pout(iig,jjg,:) + SUM(iodata(i,j,:,v))
!              END DO
!              DO v=1,3
!                 pgrid(iig,jjg,:) = pgrid(iig,jjg,:) + SUM(grid(i,j,:,v))
!              END DO
!!           END DO
!        END DO
!     END IF




 ! END IF

END SUBROUTINE block_ops
