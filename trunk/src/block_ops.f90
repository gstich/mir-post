


SUBROUTINE block_ops(iodata,corner,stride,Gdata,FUNK)
  USE globals, ONLY: varDIM,coordDIM,DIM,flen,zero,one,ax,ay,az,nx,ny,nz
  IMPLICIT NONE
  REAL(KIND=4), DIMENSION(:,:,:,:), INTENT(IN) :: iodata
  INTEGER, DIMENSION(3), INTENT(IN) :: corner
  INTEGER, DIMENSION(3,2), INTENT(IN) :: stride
  DOUBLE PRECISION, DIMENSION(:,:,:,:), INTENT(INOUT) :: Gdata
  CHARACTER(LEN=flen), INTENT(IN) :: FUNK

  REAL(KIND=4), DIMENSION(varDIM) :: temp
  REAL(KIND=4), DIMENSION(coordDIM) :: gtemp
  DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE :: dpdum
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
     k1g = stride(3,1)
     kfg = stride(3,2)
     
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
     j1g = stride(2,1)
     jfg = stride(2,2)
     k1g = stride(3,1)
     kfg = stride(3,2)

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

  CASE('PLANE_IK')
     
     ALLOCATE(dpdum(ax,az,DIM))

     !! Set the index box size for the averaging kernal
     j1g = stride(2,1)
     jfg = stride(2,2)
     
     jrngL = max(j1g-jg+1,1)
     jrngU = min(jfg-jg+1,ay)

     DO j = jrngL,jrngU
        dpdum = DBLE(iodata( : , j , : , :))
        Gdata(ig:ig+ax-1,1,kg:kg+az-1,:) = Gdata(ig:ig+ax-1,1,kg:kg+az-1,:) + dpdum
     END DO


  CASE('PLANE_JK')
     
     ALLOCATE(dpdum(ay,az,DIM))

     !! Set the index box size for the averaging kernal
     i1g = stride(1,1)
     ifg = stride(1,2)
     
     irngL = max(i1g-ig+1,1)
     irngU = min(ifg-ig+1,ax)

     DO i = irngL,irngU
        dpdum = DBLE(iodata( i , : , : , :))
        Gdata(1,jg:jg+ay-1,kg:kg+az-1,:) = Gdata(1,jg:jg+ay-1,kg:kg+az-1,:) + dpdum
     END DO

  CASE('PLANE_IJ')
     
     ALLOCATE(dpdum(ax,ay,DIM))

     !! Set the index box size for the averaging kernal
     k1g = stride(3,1)
     kfg = stride(3,2)
     
     krngL = max(k1g-kg+1,1)
     krngU = min(kfg-kg+1,az)

     DO k = krngL,krngU
        dpdum = DBLE(iodata( : , : , k , :))
        Gdata(ig:ig+ax-1,jg:jg+ay-1,1,:) = Gdata(ig:ig+ax-1,jg:jg+ay-1,1,:) + dpdum
     END DO


  CASE('BLOCK')

     Gdata(ig:ig+ax-1,jg:jg+ay-1,kg:kg+az-1,:) = iodata


  CASE('POINT')
     i1g = stride(1,1)
     ifg = stride(1,2)
     j1g = stride(2,1)
     jfg = stride(2,2)
     k1g = stride(3,1)
     kfg = stride(3,2)
     
     irngL = max(i1g-ig+1,1)
     irngU = min(ifg-ig+1,ax)
     jrngL = max(j1g-jg+1,1)
     jrngU = min(jfg-jg+1,ay)
     krngL = max(k1g-kg+1,1)
     krngU = min(kfg-kg+1,az)

     Gdata(1,1,1,:) = iodata(irngL,jrngL,krngL,:)

  CASE('POINT_W')
     i1g = stride(1,1)
     ifg = stride(1,2)
     j1g = stride(2,1)
     jfg = stride(2,2)
     k1g = stride(3,1)
     kfg = stride(3,2)
     
     irngL = max(i1g-ig+1,1)
     irngU = min(ifg-ig+1,ax)
     jrngL = max(j1g-jg+1,1)
     jrngU = min(jfg-jg+1,ay)
     krngL = max(k1g-kg+1,1)
     krngU = min(kfg-kg+1,az)

     !! Sum the entire mini-block and average outside this kernal
     DO i=irngL,irngU
        DO j=jrngL,jrngU
           DO k=krngL,krngU
              Gdata(1,1,1,:) = Gdata(1,1,1,:) + iodata(i,j,k,:)
           END DO
        END DO
     END DO

  END SELECT



END SUBROUTINE block_ops
