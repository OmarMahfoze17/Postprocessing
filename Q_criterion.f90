MODULE VARIABLE
!integer, parameter :: NX=6145 , NY=321, nz=128

END MODULE VARIABLE


PROGRAM Q_Criterion

!use VARIABLE
IMPLICIT NONE

integer, parameter :: NX=3073 , NY=321, nz=128
REAL(8),dimension(NX,NY,NZ) :: UX,UY,UZ,Q,VORT
REAL(8) :: vort_x,vort_y,vorT_z
REAL(8),dimension(NX,NY,NZ) :: dUX1,dUX2,dUX3
REAL(8),dimension(NX,NY,NZ) :: dUY1,dUY2,dUY3
REAL(8),dimension(NX,NY,NZ) :: dUZ1,dUZ2,dUZ3
real(8) :: y(Ny),x(NX),Z(NZ)
real(8) :: xlx,yly,zlz,beta
integer :: ncly,istret,i,j,k,ERROR,IFILE,FILE1,FILEN

character(len=200) :: path='/home/om1014/PhD/INCOMPACT3D/TBL/testPressGrad/cannonicalWithPressureGrad/'
character(len=200)  :: fileName


xlx=750
yly=40
zlz=15
istret=3 !2 !3
beta=0.8 !0.259065151 !1.3
ncly=2
FILE1=11
FILEN=11




call stretching(y,yly,ny,beta,istret,ncly,'yp_1.dat')

DO IFILE=FILE1,FILEN

print *, 'Read data file'
!fileName='sauve.dat'
10 format(I3.3)
write(filename, 10) IFILE


CALL readDataFiles(UX,UY,UZ,NX,NY,NZ,path,fileName,2,ERROR) !! fileType=1 << sauve.dat fileType=2 << ux,uy,uz

IF (error /=0) then 
   print *, 'File(s) not exists '
   goto 100
endif
    
print *, 'End reading data'

print *, 'Calculate Q_criterion'
CALL DERIV_X(UX,DUX1,NX,NY,NZ,XLX)
CALL DERIV_X(UY,DUY1,NX,NY,NZ,XLX)
CALL DERIV_X(UZ,DUZ1,NX,NY,NZ,XLX)

CALL DERIV_Y(UX,DUX2,NX,NY,NZ,Y)
CALL DERIV_Y(UY,DUY2,NX,NY,NZ,Y)
CALL DERIV_Y(UZ,DUZ2,NX,NY,NZ,Y)

CALL DERIV_Z(UX,DUX3,NX,NY,NZ,ZLZ)
CALL DERIV_Z(UY,DUY3,NX,NY,NZ,ZLZ)
CALL DERIV_Z(UZ,DUZ3,NX,NY,NZ,ZLZ)


DO K=1,NZ
   DO J=1,NY
      DO I=1,NX
          Q(I,J,K)=-(DUX1(i,j,k)*DUX1(i,j,k)+DUY2(i,j,k)*DUY2(i,j,k)+DUZ3(i,j,k)*DUZ3(i,j,k)&
                    +2*(DUX2(i,j,k)*DUY1(i,j,k)+DUX3(i,j,k)*DUZ1(i,j,k)+DUZ2(i,j,k)*DUY3(i,j,k)))
	  vort_z=DUX2(i,j,k)-DUY1(i,j,k)
	  vort_y=DUX3(i,j,k)-DUZ1(i,j,k)
          vort_x=DUY3(i,j,k)-DUZ2(i,j,k)
          VORT(I,J,K)=(VORT_X**2+VORT_Y**2+VORT_Z**2)**(0.5)
      ENDDO
   ENDDO
END DO

print *, 'Write Q_Criterion to file'
CALL WRITEFILE(Q,NX,NY,NZ,path, 'Q_'//filename)
CALL WRITEFILE(VORT,NX,NY,NZ,path, 'V_'//filename)
CALL WRITEFILE(UX,NX,NY,NZ,path, 'U1'//filename)
CALL WRITEFILE(UY,NX,NY,NZ,path, 'U2'//filename)
CALL WRITEFILE(UZ,NX,NY,NZ,path, 'U3'//filename)
100 ENDDO  !!! DO IFILE=FILE1,FILEN
END PROGRAM Q_Criterion



!*******************************************************************
subroutine DERIV_X(U,DU,NX,NY,NZ,XLX)
 
!*******************************************************************

IMPLICIT NONE
INTEGER :: NX,NY,NZ,I,J,K
REAL(8),dimension(NX,NY,NZ) :: U,DU
real(8) :: XLX,dx,DX2

dx=XLX/(nx-1)
DX2=2*DX
DO K=1,NZ
   DO J=1,NY
      DO I=2,NX-1
          DU(I,J,K)=(U(I+1,J,K)-U(I-1,J,K))/DX2
      ENDDO
   ENDDO
END DO

!! AT THE BOUNDARIES
DO K=1,NZ
   DO J=1,NY
          DU(1,J,K)=(U(2,J,K)-U(1,J,K))/DX 
          DU(NX,J,K)=(U(NX,J,K)-U(NX-1,J,K))/DX     
   ENDDO
END DO

end subroutine DERIV_X

!*******************************************************************
subroutine DERIV_Z(U,DU,NX,NY,NZ,ZLZ)
 
!*******************************************************************

IMPLICIT NONE
INTEGER :: NX,NY,NZ,I,J,K
REAL(8),dimension(NX,NY,NZ) :: U,DU
real(8) :: ZLZ,dZ,DZ2

dZ=ZLZ/(nZ-1)
DZ2=2*DZ
DO K=2,NZ-1
   DO J=1,NY
      DO I=1,NX
          DU(I,J,K)=(U(I,J,K+1)-U(I,J,K-1))/DZ2
      ENDDO
   ENDDO
END DO
!! AT THE BOUNDARIES
DO J=1,NY
   DO I=1,NX
          DU(I,J,1)=(U(I,J,2)-U(I,J,1))/DZ 
          DU(I,J,NZ)=(U(I,J,NZ)-U(I,J,NZ-1))/DZ     
   ENDDO
END DO

end subroutine DERIV_Z


!*******************************************************************
subroutine DERIV_Y(U,DU,NX,NY,NZ,Y)
 
!*******************************************************************

IMPLICIT NONE
INTEGER :: NX,NY,NZ,I,J,K
REAL(8),dimension(NX,NY,NZ) :: U,DU
real(8),dimension(NY) :: Y,dY

DO J=2,NY-1
   DY(J)=(Y(J+1)-Y(J-1))/2.
ENDDO
DY(1)=Y(2)-Y(1)
DY(NY)=Y(NY)-Y(NY-1)

DO K=1,NZ
   DO J=2,NY-1
      DO I=1,NX
          DU(I,J,K)=(U(I,J+1,K)-U(I,J-1,K))/(DY(J)*2.)
      ENDDO
   ENDDO
END DO
!! AT THE BOUNDARIES
DO K=1,NZ
   DO I=1,NX
          DU(I,1,K)=(U(I,2,K)-U(I,1,K))/DY(1) 
          DU(I,NY,K)=(U(I,NY,K)-U(I,NY-1,K))/DY(NY)     
   ENDDO
END DO

end subroutine DERIV_Y


!*******************************************************************
subroutine readDataFiles(UX,UY,UZ,NX,NY,NZ,path,filename,fileType,ERROR)
 
!*******************************************************************
! 


IMPLICIT NONE
INTEGER :: NX,NY,NZ,I,J,K,error,fileType
integer :: io,ii
REAL(8),dimension(nx,ny,nz) :: UX,UY,UZ
character(len=*) :: fileName,path
character(len=10),dimension(3):: fileNameTemp
integer(8) :: COUNT

logical :: exist



if (fileType==1) then
OPEN(10,FILE=trim(path)//trim(fileName),FORM='UNFORMATTED', ACCESS='DIRECT', RECL=8)

print *, nx,ny,nz
! READ UX ------------------------------------
COUNT = 1
  DO K=1,nz
     DO J=1,ny
        DO I=1,nx
           READ(10,REC=COUNT) UX(I,J,K)
           COUNT = COUNT + 1           
        ENDDO
     ENDDO
  ENDDO


! READ UY ------------------------------------
  DO K=1,nz
     DO J=1,ny
        DO I=1,nx
           READ(10,REC=COUNT) UY(I,J,K)
           COUNT = COUNT + 1
        ENDDO
     ENDDO
  ENDDO

! READ UZ ------------------------------------
  DO K=1,nz
     DO J=1,ny
        DO I=1,nx
           READ(10,REC=COUNT) UZ(I,J,K)
           COUNT = COUNT + 1
        ENDDO
     ENDDO
  ENDDO

close(10)

elseif (fileType==2) then
fileNameTemp=(/ 'ux','uy','uz' /)
Do ii=1,3
print *, trim(path)//trim(fileNameTemp(ii))//trim(fileName)
inquire(file=trim(path)//trim(fileNameTemp(ii))//trim(fileName), exist=exist)
if (exist .neqv. .True.) then 
   print *, ' file does not exist  ' , trim(fileNameTemp(ii))//trim(fileName)
   error=1
   goto 100
endif

enddo !! Do ii=1,3


Do ii=1,3 
OPEN(10,FILE=trim(path)//trim(fileNameTemp(ii))//trim(fileName),FORM='UNFORMATTED', ACCESS='DIRECT', RECL=8)
! READ UX ------------------------------------

print *, 'reading file ',trim(fileNameTemp(ii))//trim(fileName)

if (ii==1) then
COUNT = 1
  DO K=1,nz
     DO J=1,ny
        DO I=1,nx
           READ(10,REC=COUNT,IOstat=io) UX(I,J,K)
           COUNT = COUNT + 1 
           if (io .ne. 0) then 
           PRINT *, '======================================================='
           print *, 'Error occurred in reading file' , filename
           PRINT *, 'The error at i, j, k ',i,j,k
           PRINT *, '======================================================='
	   error=2
           goto 100
           endif          
        ENDDO
     ENDDO
  ENDDO
close(10)

! READ UY ------------------------------------
elseif (ii==2) then
COUNT = 1
  DO K=1,nz
     DO J=1,ny
        DO I=1,nx
           READ(10,REC=COUNT,IOstat=io) UY(I,J,K)
           COUNT = COUNT + 1
           if (io .ne. 0) then 
           PRINT *, '======================================================='
           print *, 'Error occurred in reading file' , filename
           PRINT *, 'The error at i, j, k ',i,j,k
           PRINT *, '======================================================='
	   error=2
           goto 100
           endif
        ENDDO
     ENDDO
  ENDDO

! READ UZ ------------------------------------
elseif (ii==3) then
COUNT = 1
  DO K=1,nz
     DO J=1,ny
        DO I=1,nx
           READ(10,REC=COUNT,IOstat=io) UZ(I,J,K)
           COUNT = COUNT + 1
           if (io .ne. 0) then 
           PRINT *, '======================================================='
           print *, 'Error occurred in reading file' , filename
           PRINT *, 'The error at i, j, k ',i,j,k
           PRINT *, '======================================================='
	   error=2
           goto 100
           endif
        ENDDO
     ENDDO
  ENDDO

close(10)
endif !! if (ii==1) then
enddo !! Do ii=1,3

endif !! if (fileType==1) then
100 continue
END subroutine readDataFiles

!*******************************************************************
subroutine WRITEFILE(U,NX,NY,NZ,path,fileName)

!*******************************************************************

INTEGER :: NX,NY,NZ,I,J,K
REAL(8),dimension(nx,ny,nz) :: U
character(len=*) :: fileName,path



open(10, file=trim(path)//trim(filename), FORM='UNFORMATTED', action="write",access='stream')

DO K=1,nz
   DO J=1,ny
      DO I=1,nx
         write (10) U(I,J,K)
      ENDDO
   ENDDO
ENDDO
close(10)


end subroutine WRITEFILE














!*******************************************************************
subroutine stretching(yp,yly,ny,beta,istret,ncly,fileName)
! This subroutine is taken form INCOMAPACT3D code.
! https://www.incompact3d.com/
! It determine the grid points locations when stretching is used.
 
!*******************************************************************
! 
implicit none

real(8) :: yinf,den,xnum,xcx,den1,den2,den3,den4,xnum1,cst
real(8) :: alpha,beta,pi,yly
integer :: j,istret,ncly,ny
real(8),dimension(:), allocatable :: yeta,yetai,ypi
real(8),intent(out) :: yp(*)
character(len=*) :: fileName
allocate(yeta(ny))
allocate(yetai(ny))
allocate(ypi(ny))
pi=acos(-1.)
yinf=-yly/2.
den=2.*beta*yinf
xnum=-yinf-sqrt(pi*pi*beta*beta+yinf*yinf)
alpha=abs(xnum/den)
xcx=1./beta/alpha

if (alpha.ne.0.) then
   if (istret.eq.1) yp(1)=0.
   if (istret.eq.2) yp(1)=0.
   if (istret.eq.1) yeta(1)=0.
   if (istret.eq.2) yeta(1)=-1./2.
   if (istret.eq.3) yp(1)=0.
   if (istret.eq.3) yeta(1)=-1./2.
   do j=2,ny!-1
      if (istret==1) then
         if (ncly.eq.0) yeta(j)=(j-1.)*(1./ny)
         if ((ncly.eq.1).or.(ncly.eq.2)) yeta(j)=(j-1.)*(1./(ny-1.))
      endif
      if (istret==2) then
         if (ncly.eq.0) yeta(j)=(j-1.)*(1./ny)-0.5
         if ((ncly.eq.1).or.(ncly.eq.2)) yeta(j)=(j-1.)*(1./(ny-1.))-0.5
      endif
      if (istret==3) then
         if (ncly.eq.0) yeta(j)=((j-1.)*(1./2./ny)-0.5)
         if ((ncly.eq.1).or.(ncly.eq.2)) yeta(j)=((j-1.)*(1./2./(ny-1.))-0.5)
      endif
      den1=sqrt(alpha*beta+1.)
      xnum=den1/sqrt(alpha/pi)/sqrt(beta)/sqrt(pi)
      den=2.*sqrt(alpha/pi)*sqrt(beta)*pi*sqrt(pi)
      den3=((sin(pi*yeta(j)))*(sin(pi*yeta(j)))/beta/pi)+alpha/pi
      den4=2.*alpha*beta-cos(2.*pi*yeta(j))+1.
      if ((ncly.ne.0).and.(j==ny).and.(istret==1)) then
         xnum1=0.
      else
         xnum1=(atan(xnum*tan(pi*yeta(j))))*den4/den1/den3/den
      endif
      cst=sqrt(beta)*pi/(2.*sqrt(alpha)*sqrt(alpha*beta+1.))
      if (istret==1) then
         if (yeta(j).lt.0.5) yp(j)=xnum1-cst-yinf
         if (yeta(j).eq.0.5) yp(j)=0.-yinf
         if (yeta(j).gt.0.5) yp(j)=xnum1+cst-yinf
      endif
      if (istret==2) then
         if (yeta(j).lt.0.5) yp(j)=xnum1-cst+yly
         if (yeta(j).eq.0.5) yp(j)=0.+yly
         if (yeta(j).gt.0.5) yp(j)=xnum1+cst+yly
      endif
      if (istret==3) then
         if (yeta(j).lt.0.5) yp(j)=(xnum1-cst+yly)*2.
         if (yeta(j).eq.0.5) yp(j)=(0.+yly)*2.
         if (yeta(j).gt.0.5) yp(j)=(xnum1+cst+yly)*2.
      endif
   enddo

endif
if (alpha.eq.0.) then
   yp(1)=-1.e10
   do j=2,ny
      yeta(j)=(j-1.)*(1./ny)
      yp(j)=-beta*cos(pi*yeta(j))/sin(yeta(j)*pi)
   enddo
endif
if (alpha.ne.0.) then
   do j=1,ny
      if (istret==1) then
         if (ncly.eq.0) yetai(j)=(j-0.5)*(1./ny)
         if ((ncly.eq.1).or.(ncly.eq.2)) yetai(j)=(j-0.5)*(1./(ny-1.))
      endif
      if (istret==2) then
         if (ncly.eq.0) yetai(j)=(j-0.5)*(1./ny)-0.5
         if ((ncly.eq.1).or.(ncly.eq.2)) yetai(j)=(j-0.5)*(1./(ny-1.))-0.5
      endif
      if (istret==3) then
         if (ncly.eq.0) yetai(j)=(j-0.5)*(1./2./ny)-0.5
         if ((ncly.eq.1).or.(ncly.eq.2)) yetai(j)=(j-0.5)*(1./2./(ny-1.))-0.5
      endif
      
      den1=sqrt(alpha*beta+1.)
      xnum=den1/sqrt(alpha/pi)/sqrt(beta)/sqrt(pi)
      den=2.*sqrt(alpha/pi)*sqrt(beta)*pi*sqrt(pi)
      den3=((sin(pi*yetai(j)))*(sin(pi*yetai(j)))/beta/pi)+alpha/pi
      den4=2.*alpha*beta-cos(2.*pi*yetai(j))+1.
      xnum1=(atan(xnum*tan(pi*yetai(j))))*den4/den1/den3/den
      cst=sqrt(beta)*pi/(2.*sqrt(alpha)*sqrt(alpha*beta+1.))
      if (istret==1) then
         if (yetai(j).lt.0.5) ypi(j)=xnum1-cst-yinf
         if (yetai(j).eq.0.5) ypi(j)=0.-yinf
         if (yetai(j).gt.0.5) ypi(j)=xnum1+cst-yinf
      endif
      if (istret==2) then
         if (yetai(j).lt.0.5) ypi(j)=xnum1-cst+yly
         if (yetai(j).eq.0.5) ypi(j)=0.+yly
         if (yetai(j).gt.0.5) ypi(j)=xnum1+cst+yly
      endif
      if (istret==3) then
         if (yetai(j).lt.0.5) ypi(j)=(xnum1-cst+yly)*2.
         if (yetai(j).eq.0.5) ypi(j)=(0.+yly)*2.
         if (yetai(j).gt.0.5) ypi(j)=(xnum1+cst+yly)*2.
      endif
   enddo
endif
if (alpha.eq.0.) then
   ypi(1)=-1.e10
   do j=2,ny
      yetai(j)=(j-1.)*(1./ny)
      ypi(j)=-beta*cos(pi*yetai(j))/sin(yetai(j)*pi)
   enddo
endif
open(10,file=fileName, form='formatted')
do j=1,ny
write(10,*)yp(j)
enddo
 close(10)
end subroutine stretching
