module variables
  type string
    character(len=:), allocatable :: str
  end type string
end module variables
!############################################################################
!
program verif
!
!############################################################################
  use variables
  implicit none
  integer, parameter :: mytype = KIND(0.0D0)
  integer :: nx, ny, nz,ivarName
  real(mytype),allocatable,dimension(:, :, :) :: fileData_3D
  real(mytype),allocatable,dimension(:, :) :: fileData_2D
  integer :: i,j,k,count,nfil,io
  character(len=200) :: filename_3D ,filename_2D ,filename_4D

  type(string), allocatable:: varName(:)
  type(string), allocatable:: path(:)
  logical :: exist, skipExisted
  nx=3073
  ny=321
  nz=128
  allocate(path(1))
  path=(/string('/home/om1014/PhD/INCOMPACT3D/TBL/testPressGrad/iteration_13/meanVal_200000_3rd/') /)

  allocate(varName(10))
  varName=(/string('umean.dat'),string('uumean.dat'),string('uvmean.dat'),string('vmean.dat'),&
   string('vvmean.dat'),string('vwmean.dat'),string('wmean.dat'),string('wwmean.dat'),string('pmean.dat') /)
  skipExisted=.true.
!! ***************** Allocating the variables ****************************
allocate(fileData_3D(nx,ny,nz))
allocate(fileData_2D(nx,ny))
!----------------------------------------
!goto 200
!! U Perturbation statistics !!


do ivarName=1,size(varName)


10 format(A,A,A)
 write(filename_2D, 10) path(1)%str,'2D_',varName(ivarName)%str
30 format(A,A)
 write(filename_3D, 30) path(1)%str,varName(ivarName)%str

!! check if the shrinked file exists 
inquire(file=filename_2D, exist=exist)
if (exist .and. skipExisted) then 
   print *, ' file existed  ' , filename_2D, skipExisted
   print *, ' skip the file ' , filename_3D
   goto 100
endif 
!!=======================================

 
print *,
 print *, 'Reading file  <<  ' , filename_3D 

 OPEN(20,FILE=filename_3D,FORM='UNFORMATTED',ACCESS='DIRECT', RECL=8)
 fileData_3D(:,:,:)=0._mytype
 fileData_2D(:,:)=0._mytype
  COUNT = 1
  DO K=1,nz
     DO J=1,ny
        DO I=1,nx
           READ(20,REC=COUNT,IOstat=io) fileData_3D(I,J,K)
           if (io .ne. 0) then 
           PRINT *, '======================================================='
           print *, 'Error occurred in reading file' , filename_3D
           PRINT *, 'The error at i, j, k ',i,j,k
           PRINT *, '======================================================='
           goto 100
           endif
           COUNT = COUNT + 1
        ENDDO
     ENDDO
  ENDDO
 close(20)
print *,
 print *, 'Finish reading file  <<  ' , filename_3D 
  DO K=1,nz
     DO J=1,ny
        DO I=1,nx
  	    fileData_2D(i,j)=fileData_2D(i,j)+fileData_3D(i,j,k)
        ENDDO
     ENDDO
  ENDDO

 fileData_2D(:,:)= fileData_2D(:,:)/nz
 
print *,
print *, 'Writing file  >>  ', filename_2D
 open(40, file=filename_2D, FORM='UNFORMATTED', action="write",access='stream')

  DO J=1,ny

     DO I=1,nx
        write (40) fileData_2D(I,J)
     ENDDO
  ENDDO

  close(40)
100 enddo
 
end program
!!!#############################################################
