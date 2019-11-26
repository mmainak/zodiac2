!*|-------------------------------------------------------|
      SUBROUTINE INT_MPI 
!*|-------------------------------------------------------|
use ntypes
use Domain
use Grid
use   mpi_var    
implicit none
   
CHARACTER*35 MPI_IO_NUM


! This subroutine initializes all mpi variables

      CALL MPI_INIT(IERROR)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD,NPROCES,IERROR)
      CALL MPI_COMM_RANK(MPI_COMM_WORLD,RANK,IERROR)

      write(6,*) 'MPI Initialized, NPROCS,RANK: ',NPROCES,RANK,MPI_COMM_WORLD

! Set a string to determine which input/output files to use
! When MPI is used, each process will read/write to files with the
! number of their rank (+1) appended to the end of the file.
! The string MPI_IO_NUM will be used to define the RANK+1 for each process
        IF (NPROCES.le.10) THEN
          MPI_IO_NUM=CHAR(MOD(RANK+1,10)+48)
        ELSE IF (NPROCES.le.100) THEN
          MPI_IO_NUM=CHAR(MOD(RANK+1,100)/10+48)   &
                  //CHAR(MOD(RANK+1,10)+48)
        ELSE IF (NPROCES.le.1000) THEN
          MPI_IO_NUM=CHAR(MOD(RANK+1,1000)/100+48) &
                  //CHAR(MOD(RANK+1,100)/10+48)   &
                  //CHAR(MOD(RANK+1,10)+48)
        ELSE IF (NPROCES.le.10000) THEN
          MPI_IO_NUM=CHAR(MOD(RANK+1,10000)/1000+48) &
                  //CHAR(MOD(RANK+1,1000)/100+48)   &
                  //CHAR(MOD(RANK+1,100)/10+48)     &
                  //CHAR(MOD(RANK+1,10)+48)
        ELSE
           WRITE(6,*) 'ERROR, NPROCES>10,000, Unsupported problem size'
        END IF

      RETURN
      END


      SUBROUTINE MPI_TRANSPOSE_REAL_X_TO_Z(IN_R,OUT_R)

! This subroutine is part of the MPI package for the channel flow
! Diablo package.
! Here, we define a set of ghost cells on each process
! the ghost cells contain information from the neighboring nodes
! and allow us to compute finite differences over the local gridpoints.
! We need to update the contents of the ghost zero cells at the start of
! les process

use ntypes
use Domain
use Grid
use   mpi_var
!INCLUDE 'mpif.h'
implicit none

      integer i,j,k,N
      integer(r8) :: M, SEND_NO 
      real(r8),DIMENSION(0:NXP,0:NZV-1,0:NY+1) :: IN_R
      real(r8),DIMENSION(0:NXV-1,0:NZP,0:NY+1) :: OUT_R  
!      real(r8),DIMENSION(1:NXV*(NY+2)*NZV/NP)  :: M_IN,M_OUT

      DO N=0,NP-1
       DO J=0,NY+1
        DO I=0,NXP
         DO K=N*(NZP+1)+1,(N+1)*(NZP+1) 
          M= (K-N*(NZP+1)) + (NZP+1)*I + (NXP+1)*(NZP+1)*J + (NXP+1)*(NZP+1)*(NY+2)*N
          M_IN(M) = IN_R(I,K-1,J) 
         ENDDO
        ENDDO
       ENDDO
      ENDDO


      SEND_NO = (NXP+1)*(NZP+1)*(NY+2)
!      write(6,*)'The Size of M', SEND_NO,NXV*(NY+2)*NZV/NP,shape(M_IN), shape(M_OUT),MPI_COMM_WORLD 


      CALL MPI_ALLTOALL(M_IN,SEND_NO,MPI_DOUBLE_PRECISION,M_OUT,SEND_NO,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,IERROR) 

      OUT_R(:,:,:)=0.d0

      DO N=0,NP-1
       DO J=0,NY+1
        DO I=N*(NXP+1)+1,(N+1)*(NXP+1)
         Do K=0,NZP
          M= (NZP+1)*(I-N*(NXP+1)-1) + (K+1) + (NXP+1)*(NZP+1)*J + (NXP+1)*(NZP+1)*(NY+2)*N
          OUT_R(I-1,K,J) = M_OUT(M)
         ENDDO
        ENDDO
       ENDDO
      ENDDO


      RETURN
      END

!--|*-----------------------------------------------------
      SUBROUTINE MPI_TRANSPOSE_REAL_Z_TO_X(IN_R,OUT_R)
!--|*-----------------------------------------------------

! This subroutine is part of the MPI package for the duct flow
! Diablo package.
! Here, we define a set of ghost cells on each process

use ntypes
use Domain
use Grid
use   mpi_var
!INCLUDE 'mpif.h'
implicit none

      integer i,j,k,N
      integer(r8) :: M, SEND_NO 
      real(r8),DIMENSION(0:NXP,0:NZV-1,0:NY+1) :: OUT_R
      real(r8),DIMENSION(0:NXV-1,0:NZP,0:NY+1) :: IN_R  
!      real(r8),DIMENSION(1:NXV*(NY+2)*NZV/NP)  :: M_IN,M_OUT

      DO N=0,NP-1
       DO J=0,NY+1
        DO I=N*(NXP+1)+1,(N+1)*(NXP+1)
         Do K=0,NZP
          M= (NZP+1)*(I-N*(NXP+1)-1) + (K+1) + (NXP+1)*(NZP+1)*J + (NXP+1)*(NZP+1)*(NY+2)*N  ! Bishakh
!          M = (I-N*(NXP+1)) + (NXP+1)*K + (NXP+1)*(NZP+1)*J + (NXP+1)*(NZP+1)*(NY+2)*N      ! Eric 
          M_IN(M)=IN_R(I-1,K,J) 
         ENDDO
        ENDDO
       ENDDO
      ENDDO


      SEND_NO = (NXP+1)*(NZP+1)*(NY+2)
!      write(6,*)'The Size of M', SEND_NO,NXV*(NY+2)*NZV/NP,shape(M_IN), shape(M_OUT),MPI_COMM_WORLD 


      CALL MPI_ALLTOALL(M_IN,SEND_NO,MPI_DOUBLE_PRECISION,M_OUT,SEND_NO,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,IERROR) 

      OUT_R(:,:,:)=0.d0

      DO N=0,NP-1
       DO J=0,NY+1
        DO I=0,NXP
         DO K=N*(NZP+1)+1,(N+1)*(NZP+1)
          M= (K-N*(NZP+1)) + (NZP+1)*I + (NXP+1)*(NZP+1)*J + (NXP+1)*(NZP+1)*(NY+2)*N   ! Bishakh
!           M = (I+1) + (NXP+1)*(K-N*(NZP+1)-1) + (NXP+1)*(NZP+1)*J + (NXP+1)*(NZP+1)*(NY+2)*N ! Eric
          OUT_R(I,K-1,J)=M_OUT(M) 
         ENDDO
        ENDDO
       ENDDO
      ENDDO

      RETURN
      END

!--|*-----------------------------------------------------
      SUBROUTINE MPI_TRANSPOSE_COMPLEX_Z_TO_X(IN_CZ,OUT_CX)
!--|*-----------------------------------------------------

! This subroutine is part of the MPI package for the duct flow
! Diablo package.
! Here, we define a set of ghost cells on each process

use ntypes
use Domain
use Grid
use   mpi_var
!INCLUDE 'mpif.h'
implicit none

      integer i,j,k,N
      integer(r8) :: M, SEND_NO 
      complex(r8),DIMENSION(0:NX2P,0:NZV-1,0:NY+1) :: OUT_CX
      complex(r8),DIMENSION(0:NX2V-1,0:NZP,0:NY+1) :: IN_CZ  
!      complex(r8),DIMENSION(1:NX2V*(NY+2)*NZV/NP)  :: M_IN_C,M_OUT_C

      DO N=0,NP-1
       DO J=0,NY+1
        DO I=N*(NX2P+1)+1,(N+1)*(NX2P+1)
         Do K=0,NZP
          M= (NZP+1)*(I-N*(NX2P+1)-1) + (K+1) + (NX2P+1)*(NZP+1)*J + (NX2P+1)*(NZP+1)*(NY+2)*N  ! Bishakh
!          M = (I-N*(NX2P+1)) + (NX2P+1)*K + (NX2P+1)*(NZP+1)*J + (NX2P+1)*(NZP+1)*(NY+2)*N      ! Eric 
          M_IN_C(M)=IN_CZ(I-1,K,J) 
         ENDDO
        ENDDO
       ENDDO
      ENDDO


      SEND_NO = (NX2P+1)*(NZP+1)*(NY+2)
!      write(6,*)'The Size of M', SEND_NO,NXV*(NY+2)*NZV/NP,shape(M_IN), shape(M_OUT),MPI_COMM_WORLD 


      CALL MPI_ALLTOALL(M_IN_C,SEND_NO,MPI_DOUBLE_COMPLEX,M_OUT_C,SEND_NO,MPI_DOUBLE_COMPLEX,MPI_COMM_WORLD,IERROR) 

      OUT_CX(:,:,:)=0.d0

      DO N=0,NP-1
       DO J=0,NY+1
        DO I=0,NX2P
         DO K=N*(NZP+1)+1,(N+1)*(NZP+1)
          M= (K-N*(NZP+1)) + (NZP+1)*I + (NX2P+1)*(NZP+1)*J + (NX2P+1)*(NZP+1)*(NY+2)*N   ! Bishakh
!           M = (I+1) + (NX2P+1)*(K-N*(NZP+1)-1) + (NX2P+1)*(NZP+1)*J + (NX2P+1)*(NZP+1)*(NY+2)*N ! Eric
          OUT_CX(I,K-1,J)=M_OUT_C(M) 
         ENDDO
        ENDDO
       ENDDO
      ENDDO

      RETURN
      END


      SUBROUTINE MPI_TRANSPOSE_COMPLEX_X_TO_Z(IN_C,OUT_C)

! This subroutine is part of the MPI package for the channel flow
! Diablo package.
! Here, we define a set of ghost cells on each process
! the ghost cells contain information from the neighboring nodes
! and allow us to compute finite differences over the local gridpoints.
! We need to update the contents of the ghost zero cells at the start of
! les process

use ntypes
use Domain
use Grid
use   mpi_var
!INCLUDE 'mpif.h'
implicit none

      integer i,j,k,N
      integer(r8) :: M, SEND_NO 
      complex(r8),DIMENSION(0:NX2P,0:NZV-1,0:NY+1) :: IN_C
      complex(r8),DIMENSION(0:NX2V-1,0:NZP,0:NY+1) :: OUT_C  
!      complex(r8),DIMENSION(1:NX2V*(NY+2)*NZV/NP)  :: M_IN_C,M_OUT_C

      DO N=0,NP-1
       DO J=0,NY+1
        DO I=0,NX2P
         DO K=N*(NZP+1)+1,(N+1)*(NZP+1) 
          M= (K-N*(NZP+1)) + (NZP+1)*I + (NX2P+1)*(NZP+1)*J + (NX2P+1)*(NZP+1)*(NY+2)*N
          M_IN_C(M) = IN_C(I,K-1,J) 
         ENDDO
        ENDDO
       ENDDO
      ENDDO


      SEND_NO = (NX2P+1)*(NZP+1)*(NY+2)
!      write(6,*)'The Size of M', SEND_NO,NXV*(NY+2)*NZV/NP,shape(M_IN), shape(M_OUT),MPI_COMM_WORLD 


      CALL MPI_ALLTOALL(M_IN_C,SEND_NO,MPI_DOUBLE_COMPLEX,M_OUT_C,SEND_NO,MPI_DOUBLE_COMPLEX,MPI_COMM_WORLD,IERROR) 

      OUT_C(:,:,:)=(0.d0,0.d0)

      DO N=0,NP-1
       DO J=0,NY+1
        DO I=N*(NX2P+1)+1,(N+1)*(NX2P+1)
         Do K=0,NZP
          M= (NZP+1)*(I-N*(NX2P+1)-1) + (K+1) + (NX2P+1)*(NZP+1)*J + (NX2P+1)*(NZP+1)*(NY+2)*N
          OUT_C(I-1,K,J) = M_OUT_C(M)
         ENDDO
        ENDDO
       ENDDO
      ENDDO


      RETURN
      END
   
      SUBROUTINE REAL_FOURIER_TRANS_U1 (flag) 
use ntypes 
use Domain 
use Grid
use TIME_STEP_VAR
use run_variable 
use   mpi_var
      implicit none
      integer :: i,j,k
      logical :: flag             

     
      if (flag) then
      S1X=U1X
      CALL MPI_TRANSPOSE_REAL_X_TO_Z(S1X,S1Z)
      varp(:,:,:) = 0.d0
      DO I=0,NXM
       varp(I,:,:)=S1Z(I,:,:)
      ENDDO
      CALL FFT_X_TO_FOURIER_OP(varp,cvarp,0,NY+1,0,NZP)
      DO I=0,NKX
      CS1Z(I,:,:)=cvarp(I,:,:)
      ENDDO
      CS1Z(NKX+1:NX2V-1,:,:)=(0.0,0.0)
      CALL MPI_TRANSPOSE_COMPLEX_Z_TO_X(CS1Z,CS1X)
      call allocation_u(.false.)
      CU1X=CS1X

      else 

      
      CS1X=CU1X
      call allocation_u(.true.)
      CALL MPI_TRANSPOSE_COMPLEX_X_TO_Z(CS1X,CS1Z)
      cvarp(:,:,:)=(0.d0,0.d0)
      DO I=0,NKX
      cvarp(I,:,:)=CS1Z(I,:,:)
      ENDDO
      CALL FFT_X_TO_PHYSICAL_OP(cvarp,varp,0,NY+1,0,NZP)
      DO I=0,NXM
      S1Z(I,:,:)=varp(I,:,:)
      ENDDO
      S1Z(NXM+1:NXV-1,:,:)=0.0
      CALL MPI_TRANSPOSE_REAL_Z_TO_X(S1Z,S1X)
      U1X=S1X
      endif

      RETURN 
      END  

      SUBROUTINE REAL_FOURIER_TRANS_U2 (flag) 
use ntypes 
use Domain 
use Grid
use TIME_STEP_VAR
use run_variable 
use   mpi_var
      implicit none
      integer :: i,j,k
      logical :: flag             

     
      if (flag) then
      S1X=U2X
      CALL MPI_TRANSPOSE_REAL_X_TO_Z(S1X,S1Z)
      varp(:,:,:) = 0.d0
      DO I=0,NXM
       varp(I,:,:)=S1Z(I,:,:)
      ENDDO
      CALL FFT_X_TO_FOURIER_OP(varp,cvarp,0,NY+1,0,NZP)
      DO I=0,NKX
      CS1Z(I,:,:)=cvarp(I,:,:)
      ENDDO
      CS1Z(NKX+1:NX2V-1,:,:)=(0.0,0.0)
      CALL MPI_TRANSPOSE_COMPLEX_Z_TO_X(CS1Z,CS1X)
      call allocation_v(.false.)
      CU2X=CS1X

      else 

      
      CS1X=CU2X
      call allocation_v(.true.)
      CALL MPI_TRANSPOSE_COMPLEX_X_TO_Z(CS1X,CS1Z)
      cvarp(:,:,:)=(0.d0,0.d0)
      DO I=0,NKX
      cvarp(I,:,:)=CS1Z(I,:,:)
      ENDDO
      CALL FFT_X_TO_PHYSICAL_OP(cvarp,varp,0,NY+1,0,NZP)
      DO I=0,NXM
      S1Z(I,:,:)=varp(I,:,:)
      ENDDO
      S1Z(NXM+1:NXV-1,:,:)=0.0
      CALL MPI_TRANSPOSE_REAL_Z_TO_X(S1Z,S1X)
      U2X=S1X
      endif

      RETURN 
      END  
 
      SUBROUTINE REAL_FOURIER_TRANS_U3 (flag) 
use ntypes 
use Domain 
use Grid
use TIME_STEP_VAR
use run_variable 
use   mpi_var
      implicit none
      integer :: i,j,k
      logical :: flag             

     
      if (flag) then
      S1X=U3X
      CALL MPI_TRANSPOSE_REAL_X_TO_Z(S1X,S1Z)
      varp(:,:,:) = 0.d0
      DO I=0,NXM
       varp(I,:,:)=S1Z(I,:,:)
      ENDDO
      CALL FFT_X_TO_FOURIER_OP(varp,cvarp,0,NY+1,0,NZP)
      DO I=0,NKX
      CS1Z(I,:,:)=cvarp(I,:,:)
      ENDDO
      CS1Z(NKX+1:NX2V-1,:,:)=(0.0,0.0)
      CALL MPI_TRANSPOSE_COMPLEX_Z_TO_X(CS1Z,CS1X)
      call allocation_w(.false.)
      CU3X=CS1X

      else 

      
      CS1X=CU3X
      call allocation_w(.true.)
      CALL MPI_TRANSPOSE_COMPLEX_X_TO_Z(CS1X,CS1Z)
      cvarp(:,:,:)=(0.d0,0.d0)
      DO I=0,NKX
      cvarp(I,:,:)=CS1Z(I,:,:)
      ENDDO
      CALL FFT_X_TO_PHYSICAL_OP(cvarp,varp,0,NY+1,0,NZP)
      DO I=0,NXM
      S1Z(I,:,:)=varp(I,:,:)
      ENDDO
      S1Z(NXM+1:NXV-1,:,:)=0.0
      CALL MPI_TRANSPOSE_REAL_Z_TO_X(S1Z,S1X)
      U3X=S1X
      endif

      RETURN 
      END  

      SUBROUTINE REAL_FOURIER_TRANS_P (flag) 
use ntypes 
use Domain 
use Grid
use TIME_STEP_VAR
use run_variable 
use   mpi_var
      implicit none
      integer :: i,j,k
      logical :: flag             

     
      if (flag) then
      S1X=PX
      CALL MPI_TRANSPOSE_REAL_X_TO_Z(S1X,S1Z)
      varp(:,:,:) = 0.d0
      DO I=0,NXM
       varp(I,:,:)=S1Z(I,:,:)
      ENDDO
      CALL FFT_X_TO_FOURIER_OP(varp,cvarp,0,NY+1,0,NZP)
      DO I=0,NKX
      CS1Z(I,:,:)=cvarp(I,:,:)
      ENDDO
      CS1Z(NKX+1:NX2V-1,:,:)=(0.0,0.0)
      CALL MPI_TRANSPOSE_COMPLEX_Z_TO_X(CS1Z,CS1X)
      call allocation_p(.false.)
      CPX=CS1X

      else 

      
      CS1X=CPX
      call allocation_p(.true.)
      CALL MPI_TRANSPOSE_COMPLEX_X_TO_Z(CS1X,CS1Z)
      cvarp(:,:,:)=(0.d0,0.d0)
      DO I=0,NKX
      cvarp(I,:,:)=CS1Z(I,:,:)
      ENDDO
      CALL FFT_X_TO_PHYSICAL_OP(cvarp,varp,0,NY+1,0,NZP)
      DO I=0,NXM
      S1Z(I,:,:)=varp(I,:,:)
      ENDDO
      S1Z(NXM+1:NXV-1,:,:)=0.0
      CALL MPI_TRANSPOSE_REAL_Z_TO_X(S1Z,S1X)
      PX=S1X
      endif

      RETURN 
      END  
 
 
      SUBROUTINE REAL_FOURIER_TRANS_TH (flag) 
use ntypes 
use Domain 
use Grid
use TIME_STEP_VAR
use run_variable 
use   mpi_var
      implicit none
      integer :: i,j,k,N
      logical :: flag             

      
      if (flag) then
      S1X=THX(:,:,:,1)
      CALL MPI_TRANSPOSE_REAL_X_TO_Z(S1X,S1Z)
      varp(:,:,:) = 0.d0
      DO I=0,NXM
       varp(I,:,:)=S1Z(I,:,:)
      ENDDO
      CALL FFT_X_TO_FOURIER_OP(varp,cvarp,0,NY+1,0,NZP)
      DO I=0,NKX
      CS1Z(I,:,:)=cvarp(I,:,:)
      ENDDO
      CS1Z(NKX+1:NX2V-1,:,:)=(0.0,0.0)
      CALL MPI_TRANSPOSE_COMPLEX_Z_TO_X(CS1Z,CS1X)
      
      IF (N_TH .gt. 1) THEN
      S2X=THX(:,:,:,2)
      CALL MPI_TRANSPOSE_REAL_X_TO_Z(S2X,S2Z)
      varp(:,:,:) = 0.d0
      DO I=0,NXM
       varp(I,:,:)=S2Z(I,:,:)
      ENDDO
      CALL FFT_X_TO_FOURIER_OP(varp,cvarp,0,NY+1,0,NZP)
      DO I=0,NKX
      CS2Z(I,:,:)=cvarp(I,:,:)
      ENDDO
      CS2Z(NKX+1:NX2V-1,:,:)=(0.0,0.0)
      CALL MPI_TRANSPOSE_COMPLEX_Z_TO_X(CS2Z,CS2X)
      ENDIF

      call allocation_th(.false.)
      CTHX(:,:,:,1)=CS1X
      IF (N_TH .gt. 1) THEN
      CTHX(:,:,:,2)=CS2X 
      ENDIF

      else 
      
      CS1X=CTHX(:,:,:,1)
      CALL MPI_TRANSPOSE_COMPLEX_X_TO_Z(CS1X,CS1Z)
      cvarp(:,:,:)=(0.d0,0.d0)
      DO I=0,NKX
      cvarp(I,:,:)=CS1Z(I,:,:)
      ENDDO
      CALL FFT_X_TO_PHYSICAL_OP(cvarp,varp,0,NY+1,0,NZP)
      DO I=0,NXM
      S1Z(I,:,:)=varp(I,:,:)
      ENDDO
      S1Z(NXM+1:NXV-1,:,:)=0.0
      CALL MPI_TRANSPOSE_REAL_Z_TO_X(S1Z,S1X)

      IF (N_TH .gt. 1) THEN
      CS2X=CTHX(:,:,:,2)
      CALL MPI_TRANSPOSE_COMPLEX_X_TO_Z(CS2X,CS2Z)
      cvarp(:,:,:)=(0.d0,0.d0)
      DO I=0,NKX
      cvarp(I,:,:)=CS2Z(I,:,:)
      ENDDO
      CALL FFT_X_TO_PHYSICAL_OP(cvarp,varp,0,NY+1,0,NZP)
      DO I=0,NXM
      S2Z(I,:,:)=varp(I,:,:)
      ENDDO
      S2Z(NXM+1:NXV-1,:,:)=0.0
      CALL MPI_TRANSPOSE_REAL_Z_TO_X(S2Z,S2X)
      ENDIF


      call allocation_th(.true.)
      THX(:,:,:,1)=S1X
      IF (N_TH .gt. 1) THEN
      THX(:,:,:,2)=S2X
      ENDIF

      endif

      RETURN 
      END  
 
      SUBROUTINE REAL_FOURIER_TRANS_R1 (flag) 
use ntypes 
use Domain 
use Grid
use TIME_STEP_VAR
use run_variable 
use   mpi_var
      implicit none
      integer :: i,j,k
      logical :: flag             

     
      if (flag) then
      S1X=R1X
      CALL MPI_TRANSPOSE_REAL_X_TO_Z(S1X,S1Z)
      varp(:,:,:) = 0.d0
      DO I=0,NXM
       varp(I,:,:)=S1Z(I,:,:)
      ENDDO
      CALL FFT_X_TO_FOURIER_OP(varp,cvarp,0,NY+1,0,NZP)
      DO I=0,NKX
      CS1Z(I,:,:)=cvarp(I,:,:)
      ENDDO
      CS1Z(NKX+1:NX2V-1,:,:)=(0.0,0.0)
      CALL MPI_TRANSPOSE_COMPLEX_Z_TO_X(CS1Z,CS1X)
      call allocation_R1 (.false.)
      CR1X=CS1X

      else 

      
      CS1X=CR1X
      call allocation_R1(.true.)
      CALL MPI_TRANSPOSE_COMPLEX_X_TO_Z(CS1X,CS1Z)
      cvarp(:,:,:)=(0.d0,0.d0)
      DO I=0,NKX
      cvarp(I,:,:)=CS1Z(I,:,:)
      ENDDO
      CALL FFT_X_TO_PHYSICAL_OP(cvarp,varp,0,NY+1,0,NZP)
      DO I=0,NXM
      S1Z(I,:,:)=varp(I,:,:)
      ENDDO
      S1Z(NXM+1:NXV-1,:,:)=0.0
      CALL MPI_TRANSPOSE_REAL_Z_TO_X(S1Z,S1X)
      R1X=S1X
      endif

      RETURN 
      END  


      SUBROUTINE REAL_FOURIER_TRANS_R2 (flag) 
use ntypes 
use Domain 
use Grid
use TIME_STEP_VAR
use run_variable 
use   mpi_var
      implicit none
      integer :: i,j,k
      logical :: flag             

     
      if (flag) then
      S1X=R2X
      CALL MPI_TRANSPOSE_REAL_X_TO_Z(S1X,S1Z)
      varp(:,:,:) = 0.d0
      DO I=0,NXM
       varp(I,:,:)=S1Z(I,:,:)
      ENDDO
      CALL FFT_X_TO_FOURIER_OP(varp,cvarp,0,NY+1,0,NZP)
      DO I=0,NKX
      CS1Z(I,:,:)=cvarp(I,:,:)
      ENDDO
      CS1Z(NKX+1:NX2V-1,:,:)=(0.0,0.0)
      CALL MPI_TRANSPOSE_COMPLEX_Z_TO_X(CS1Z,CS1X)
      call allocation_R2 (.false.)
      CR2X=CS1X

      else 

      
      CS1X=CR2X
      call allocation_R2(.true.)
      CALL MPI_TRANSPOSE_COMPLEX_X_TO_Z(CS1X,CS1Z)
      cvarp(:,:,:)=(0.d0,0.d0)
      DO I=0,NKX
      cvarp(I,:,:)=CS1Z(I,:,:)
      ENDDO
      CALL FFT_X_TO_PHYSICAL_OP(cvarp,varp,0,NY+1,0,NZP)
      DO I=0,NXM
      S1Z(I,:,:)=varp(I,:,:)
      ENDDO
      S1Z(NXM+1:NXV-1,:,:)=0.0
      CALL MPI_TRANSPOSE_REAL_Z_TO_X(S1Z,S1X)
      R2X=S1X
      endif

      RETURN 
      END  

      SUBROUTINE REAL_FOURIER_TRANS_R3 (flag) 
use ntypes 
use Domain 
use Grid
use TIME_STEP_VAR
use run_variable 
use   mpi_var
      implicit none
      integer :: i,j,k
      logical :: flag             

     
      if (flag) then
      S1X=R3X
      CALL MPI_TRANSPOSE_REAL_X_TO_Z(S1X,S1Z)
      varp(:,:,:) = 0.d0
      DO I=0,NXM
       varp(I,:,:)=S1Z(I,:,:)
      ENDDO
      CALL FFT_X_TO_FOURIER_OP(varp,cvarp,0,NY+1,0,NZP)
      DO I=0,NKX
      CS1Z(I,:,:)=cvarp(I,:,:)
      ENDDO
      CS1Z(NKX+1:NX2V-1,:,:)=(0.0,0.0)
      CALL MPI_TRANSPOSE_COMPLEX_Z_TO_X(CS1Z,CS1X)
      call allocation_R3 (.false.)
      CR3X=CS1X

      else 

      
      CS1X=CR3X
      call allocation_R3(.true.)
      CALL MPI_TRANSPOSE_COMPLEX_X_TO_Z(CS1X,CS1Z)
      cvarp(:,:,:)=(0.d0,0.d0)
      DO I=0,NKX
      cvarp(I,:,:)=CS1Z(I,:,:)
      ENDDO
      CALL FFT_X_TO_PHYSICAL_OP(cvarp,varp,0,NY+1,0,NZP)
      DO I=0,NXM
      S1Z(I,:,:)=varp(I,:,:)
      ENDDO
      S1Z(NXM+1:NXV-1,:,:)=0.0
      CALL MPI_TRANSPOSE_REAL_Z_TO_X(S1Z,S1X)
      R3X=S1X
      endif

      RETURN 
      END  

      SUBROUTINE REAL_FOURIER_TRANS_Rth (flag) 
use ntypes 
use Domain 
use Grid
use TIME_STEP_VAR
use run_variable 
use   mpi_var
      implicit none
      integer :: i,j,k,N
      logical :: flag             

      
      if (flag) then
      S1X=RTHX(:,:,:,1)
      CALL MPI_TRANSPOSE_REAL_X_TO_Z(S1X,S1Z)
      varp(:,:,:) = 0.d0
      DO I=0,NXM
       varp(I,:,:)=S1Z(I,:,:)
      ENDDO
      CALL FFT_X_TO_FOURIER_OP(varp,cvarp,0,NY+1,0,NZP)
      DO I=0,NKX
      CS1Z(I,:,:)=cvarp(I,:,:)
      ENDDO
      CS1Z(NKX+1:NX2V-1,:,:)=(0.0,0.0)
      CALL MPI_TRANSPOSE_COMPLEX_Z_TO_X(CS1Z,CS1X)
      
      IF (N_TH .gt. 1) THEN
      S2X=RTHX(:,:,:,2)
      CALL MPI_TRANSPOSE_REAL_X_TO_Z(S2X,S2Z)
      varp(:,:,:) = 0.d0
      DO I=0,NXM
       varp(I,:,:)=S2Z(I,:,:)
      ENDDO
      CALL FFT_X_TO_FOURIER_OP(varp,cvarp,0,NY+1,0,NZP)
      DO I=0,NKX
      CS2Z(I,:,:)=cvarp(I,:,:)
      ENDDO
      CS2Z(NKX+1:NX2V-1,:,:)=(0.0,0.0)
      CALL MPI_TRANSPOSE_COMPLEX_Z_TO_X(CS2Z,CS2X)
      ENDIF

      call allocation_Rth(.false.)
      CRTHX(:,:,:,1)=CS1X
      IF (N_TH .gt. 1) THEN
       CRTHX(:,:,:,2)=CS2X
      ENDIF

      else 

      
      CS1X=CRTHX(:,:,:,1)
      CALL MPI_TRANSPOSE_COMPLEX_X_TO_Z(CS1X,CS1Z)
      cvarp(:,:,:)=(0.d0,0.d0)
      DO I=0,NKX
      cvarp(I,:,:)=CS1Z(I,:,:)
      ENDDO
      CALL FFT_X_TO_PHYSICAL_OP(cvarp,varp,0,NY+1,0,NZP)
      DO I=0,NXM
      S1Z(I,:,:)=varp(I,:,:)
      ENDDO
      S1Z(NXM+1:NXV-1,:,:)=0.0
      CALL MPI_TRANSPOSE_REAL_Z_TO_X(S1Z,S1X)

      IF (N_TH .gt. 1) THEN
      CS2X=CRTHX(:,:,:,2)
      CALL MPI_TRANSPOSE_COMPLEX_X_TO_Z(CS2X,CS2Z)
      cvarp(:,:,:)=(0.d0,0.d0)
      DO I=0,NKX
      cvarp(I,:,:)=CS2Z(I,:,:)
      ENDDO
      CALL FFT_X_TO_PHYSICAL_OP(cvarp,varp,0,NY+1,0,NZP)
      DO I=0,NXM
      S2Z(I,:,:)=varp(I,:,:)
      ENDDO
      S2Z(NXM+1:NXV-1,:,:)=0.0
      CALL MPI_TRANSPOSE_REAL_Z_TO_X(S2Z,S2X)
      ENDIF

      call allocation_Rth(.true.)

      RTHX(:,:,:,1)=S1X
      IF (N_TH .gt. 1) THEN
      RTHX(:,:,:,2)=S2X
      ENDIF

      endif

      RETURN 
      END  
 
 
      
      SUBROUTINE MPI_BCAST_REAL(IN_R,Z_SIZE,Y_SIZE)

! This subroutine is part of the MPI package for the channel flow
! Diablo package.
! Here, we define a set of ghost cells on each process
! the ghost cells contain information from the neighboring nodes
! and allow us to compute finite differences over the local gridpoints.
! We need to update the contents of the ghost zero cells at the start of
! les process

use ntypes
use Domain
use Grid
use   mpi_var
!INCLUDE 'mpif.h'
implicit none

      integer i,j,k,N,send_rank,Z_SIZE,Y_SIZE
      real(r8),dimension(0:Z_SIZE*Y_SIZE-1) :: IN_R

      send_rank = 0


      CALL MPI_BCAST(IN_R,Z_SIZE*Y_SIZE,MPI_DOUBLE_PRECISION,send_rank,MPI_COMM_WORLD,IERROR) 

      RETURN
      END


      SUBROUTINE MPI_BCAST_COMPLEX(IN_R,Z_SIZE,Y_SIZE)

! This subroutine is part of the MPI package for the channel flow
! Diablo package.
! Here, we define a set of ghost cells on each process
! the ghost cells contain information from the neighboring nodes
! and allow us to compute finite differences over the local gridpoints.
! We need to update the contents of the ghost zero cells at the start of
! les process

use ntypes
use Domain
use Grid
use   mpi_var
!INCLUDE 'mpif.h'
implicit none

      integer i,j,k,N,send_rank,Z_SIZE,Y_SIZE
      complex(r8),dimension(0:Z_SIZE*Y_SIZE-1) :: IN_R

      send_rank = 0


      CALL MPI_BCAST(IN_R,Z_SIZE*Y_SIZE,MPI_DOUBLE_COMPLEX,send_rank,MPI_COMM_WORLD,IERROR)

      RETURN
      END

!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE MPI_COMBINE_STATS(STAT,Z_SIZE,Y_SIZE)
!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
use ntypes
use Domain
use Grid
use   mpi_var
implicit none    
   integer :: I, J, Z_SIZE,Y_SIZE,send_rank
   real(r8),dimension(0:Z_SIZE*Y_SIZE-1)      :: STAT
   real(r8),dimension(0:Z_SIZE*Y_SIZE-1,1:NP) :: STAT_TMP
   

     send_rank = 0
     CALL MPI_GATHER(STAT,Z_SIZE*Y_SIZE,MPI_DOUBLE_PRECISION,STAT_TMP,Z_SIZE*Y_SIZE,MPI_DOUBLE_PRECISION,send_rank,MPI_COMM_WORLD,IERROR)
     IF ( RANK .EQ. 0) THEN
      STAT(:) = 0D0
        DO I = 1,NP
          DO J =0, Z_SIZE*Y_SIZE-1
           STAT(J) = STAT(J) + STAT_TMP(J,I)
          ENDDO
        ENDDO
      ENDIF
 
      RETURN
      END 


!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
 SUBROUTINE MPI_COURANT(DT)
!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
use ntypes
use Domain
use Grid
use mpi_var
implicit none


  REAL(r8),INTENT(INOUT) :: DT
  REAL(r8),DIMENSION(1:NP) :: IPACK
  REAL(r8) :: OPACK

! READ IN LOCAL DT BASED ON CFL
      OPACK = DT

! SEND ALL LOCAL DT'S TO EVERY PROCESS
      CALL MPI_ALLGATHER(OPACK,1,MPI_DOUBLE_PRECISION,IPACK,1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,IERROR)

      DT = MINVAL(IPACK)

      RETURN
      END


!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
 SUBROUTINE FILTER_VAR(filter_type)
!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
use ntypes
use Domain
use Grid
use mpi_var
use run_variable, only: S1X,S1Z
implicit none

     integer :: filter_type

     CALL MPI_TRANSPOSE_REAL_X_TO_Z(S1X,S1Z)
     CALL les_filter_chan (S1Z,0,NZP,0,NY+1,filter_type)
     CALL MPI_TRANSPOSE_REAL_Z_TO_X(S1Z,S1X)


     RETURN
     END

!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE MPI_INST_PROF(OPACK,IPACK)
!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
use ntypes
use Domain
use Grid
use mpi_var
use run_variable
      REAL(8),DIMENSION(1:NY*NP) :: IPACK
      REAL(8),DIMENSION(1:NY) :: OPACK

! READ IN LOCAL DT BASED ON CFL
!      OPACK = DT

! SEND ALL LOCAL DT'S TO EVERY PROCESS
CALL MPI_ALLGATHER(OPACK,NY,MPI_DOUBLE_PRECISION,IPACK,NY,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,IERROR)

!      DT = MINVAL(IPACK)

      RETURN
      END

!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE MPI_COMBINE(STAT,temp,SIZE_1,SIZE_2)
!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
use ntypes
use Domain
use Grid
use   mpi_var
implicit none
   integer :: I, J, SIZE_1, SIZE_2, send_rank,en,st
   real(r8),dimension(1:SIZE_1,1:SIZE_2)      :: STAT
   real(r8),dimension(1:SIZE_1,1:SIZE_2,1:NP) :: STAT_TMP
   real(r8),dimension(1:SIZE_1*NP,1:SIZE_2)      :: temp


     send_rank = 0
     CALL MPI_GATHER(STAT,SIZE_1*SIZE_2,MPI_DOUBLE_PRECISION,STAT_TMP,SIZE_1*SIZE_2, &
                      MPI_DOUBLE_PRECISION,send_rank,MPI_COMM_WORLD,IERROR)




     IF ( RANK .EQ. 0) THEN
        temp(:,:) = 0.0d0
        st=1
        DO I = 1,NP
          en=st+SIZE_1-1
          DO J =1,SIZE_2
           temp(st:en,J) =  STAT_TMP(:,J,I)
          ENDDO
          st=en+1
        ENDDO
      ENDIF



      RETURN
      END
