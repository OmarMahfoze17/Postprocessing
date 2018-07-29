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
  integer nx, ny, nz,shrinkStep,error
  integer :: file1,filen,nfiles,nclx,ncly,nclz,istret,ivarName
  integer :: nxShr,nyShr,nzShr,jj,jjj
  integer(4) :: icrfile, ifile, dig1, dig2, dig3, dig4
  real(4) :: xlx,yly,zlz,dt,dx,dy,dz,ydumy
  real(mytype),allocatable,dimension(:, :, :) :: fileData
  real(4), allocatable :: yp(:),y1(:),y3(:)
  integer :: i,j,k,count,nfil,ifileNumb,io
  character(3) :: chits !IF THE DATA ARE STORED WITH 3 DIGITS, IE UX001,UX002,ETC.
  character(len=20) :: filename ,filenameNew 
  type(string), allocatable:: varName(:)
  type(string), allocatable:: path(:)
  logical :: exist, skipExisted
  nx=3073
  ny=321
  nz=128
  xlx=750
  yly=40
  zlz=15
  nclx=2
  ncly=2
  nclz=0
  file1=55
  filen=55
  nfiles=filen-file1+1
  istret=1.
  shrinkStep=1
  allocate(path(1))
  path=(/string('/home/om1014/PhD/INCOMPACT3D/TBL/testPressGrad/cannonicalWithPressureGrad/meanVal_204000_2nd/') /)
  allocate(varName(1))
  varName=(/string('pmean.dat')/)
  skipExisted=.true.
!! ***************** Allocating the variables ****************************
allocate(fileData(nx,ny,nz))
!----------------------------------------
!goto 200
!! U Perturbation statistics !!


do ivarName=1,size(varName)
do ifileNumb=file1,filen

10 format(A,A,I3.3)
 write(filenameNew, 10) varName(ivarName)%str,'_Orignal'
30 format(A,I3.3)
 write(filename, 30) varName(ivarName)%str
print *, filename, filenameNew


!! check if the shrinked file exists 
inquire(file=path(1)%str//filenameNew, exist=exist)
if (exist .and. skipExisted) then 
   print *, ' file existed  ' , filenameNew, skipExisted
   print *, ' skip the file ' , filename
   goto 100
endif 
!!=======================================

 print *, path(1)%str
print *,
 print *, 'Reading file  <<  ' , filename 

 OPEN(20,FILE=path(1)%str//filename,FORM='UNFORMATTED',ACCESS='DIRECT', RECL=8)
 fileData(:,:,:)=0._mytype
  COUNT = 1
  DO K=1,nz
     DO J=1,ny
        DO I=1,nx
           READ(20,REC=COUNT,IOstat=io) fileData(I,J,K)
           if (io .ne. 0) then 
           PRINT *, '======================================================='
           print *, 'Error occurred in reading file' , filename
           PRINT *, 'The error at i, j, k ',i,j,k
           PRINT *, '======================================================='
           goto 100
           endif
           COUNT = COUNT + 1
        ENDDO
     ENDDO
  ENDDO
 close(20)
call rename(path(1)%str//filename, path(1)%str//filenameNew,error)

if (error/=0) then 
    print *, 'Error: The orignal file did not renamed'

stop
endif
fileData=fileData/2.  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
print *,
print *, 'Writing file  >>  ', filename
 open(40, file=path(1)%str//filename, FORM='UNFORMATTED', action="write",access='stream')
nzShr=0
  DO K=1,nz,shrinkStep
     nyShr=0
     DO J=1,ny,shrinkStep
        nxShr=0
        DO I=1,nx,shrinkStep
           nxShr=nxShr+1
           write (40) fileData(I,J,K)
        ENDDO
        nyShr=nyShr+1
     ENDDO
     nzShr=nzShr+1
  ENDDO
  close(40)
100 enddo
enddo 
200 nxShr=nxShr 
nyShr=nyShr
nzShr=nzShr
print *,
print *, 'Shrinked domain size ',nxShr,nyShr,nzShr

end program verif
