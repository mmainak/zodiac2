! C******************************************************************************|
! C diablo.f -> DNS In A Box, Laptop Optimized                       VERSION 0.9
! C
! C This Fortran 77 code computes incompressible flow in a box.
! C
! C Primative variables (u,v,w,p) are used, and continuity is enforced with a
! C fractional step algorithm.
! C
! C SPATIAL DERIVATIVES:
! C   0, 1, 2, or 3 directions are taken to be periodic and handled spectrally
! C   (these cases are referred to as the "periodic", "channel", "duct", and
! C    "cavity" cases respectively).
! C   The remaining directions are taken to be bounded by walls and handled with
! C   momentum- and energy-conserving second-order central finite differences.
! C
! C TIME ADVANCEMENT
! C   Two main approaches are implemented:
! C     1. RKW3 on nonlinear terms and CN on viscous terms over each RK substep.
! C     2. RKW3 on y-derivative terms and CN on other terms over each RK substep.
! C
! C The emphasis in this introductory code is on code simplicity:
! C   -> All variables are in core.
! C   -> The code is not explicitly designed for use with either MPI or SMP.
! C   -> Overindexing is not used.
! C A few simple high-performance programming constructs are used:
! C   -> The inner 2 loops are broken out in such a way as to enable out-of-order
! C      execution of these loops as much as possible, thereby leveraging
! C      vector and superscalar CPU architectures.
! C   -> The outer loops are fairly long (including as many operations as
! C      possible inside on a single J plane of data) in order to make effective
! C      use of cache.
! C Multiple time advancement algorithms are implemented for the periodic,
! C channel, duct, and cavity cases in order to compare their efficiency for
! C various flows on different computational architectures.  In a few of the
! C algorithms, the code is broken out fully into a PASS1/PASS2 architecture
! C to maximize the efficient use of cache.
! C
! C This code was developed as a joint project in MAE 223 (CFD), taught by
! C Thomas Bewley, at UC San Diego (spring 2001, spring 2005).
! C Primary contributions follow:
! C Thomas Bewley was the chief software architect
! C John R. Taylor wrote the channel flow solvers
! C******************************************************************************|
! C
! C This code is free software; you can redistribute it and/or modify it
! C under the terms of the GNU General Public License as published by the
! C Free Software Foundation; either version 2 of the License, or (at your
! C option) any later version. This code is distributed in the hope that it
! C will be useful, but WITHOUT ANY WARRANTY; without even the implied
! C warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! C GNU General Public License for more details. You should have received a
! C copy of the GNU General Public License along with this code; if not,
! C write to the Free Software Foundation, Inc., 59 Temple Place - Suite
! C 330, Boston, MA 02111-1307, USA.
! C
! C******************************************************************************|
! 
! C----|--.---------.---------.---------.---------.---------.---------.-|-------|
      PROGRAM DIABLO
!      INCLUDE 'header_duct'
use ntypes
use Domain
use Grid
use Fft_var
use TIME_STEP_VAR
use run_variable
use omp_lib
use mpi_var, only: rank
implicit none

      INTEGER N
      REAL*8  int_time 
      

      WRITE(6,*) 
      WRITE(6,*) '             ****** WELCOME TO DIABLO ******'
      WRITE(6,*)

! c      CALL HELLO

      CALL INITIALIZE
! Initialize START_TIME for run timing



      CALL DATE_AND_TIME (VALUES=TIME_ARRAY)
      START_TIME=TIME_ARRAY(5)*3600+TIME_ARRAY(6)*60 &
        +TIME_ARRAY(7)+TIME_ARRAY(8)*0.001

! C A flag to determine if we are considering the first time-step
      FIRST_TIME=.TRUE. 
     
      INT_TREAT = .FALSE. 
      int_time =  TIME      ! #### chnage has been made
      
      DO TIME_STEP = TIME_STEP+1, TIME_STEP+N_TIME_STEPS
!        rktime = omp_get_wtime ( )

        IF (rank .eq. 0) then
        WRITE(6,*) 'Now beginning TIME_STEP = ',TIME_STEP

        open(99,file='out_screen.txt',form='formatted', &
      status='unknown', position='append')
       WRITE(99,*) 'Now beginning TIME_STEP = ',TIME_STEP
       close(99)
        ENDIF

       DO N=1,N_TH 
        IF (INT_TREAT ) THEN
         If ( (TIME-int_time) .LT. PI/2.0) then
            RI_TAU(N) = RI_FINAL(N)*(TIME-int_time)/(PI/2.0)
          ELSE
            RI_TAU(N) = RI_FINAL(N)
          ENDIF    

          If ( (TIME-int_time) .LT. 10.0/OMEGA0) then
           FACT_AMP = (TIME-int_time)/(10.0/OMEGA0)
           IF (rank .eq. 0) then
            WRITE(6,*) 'RAMPING value', FACT_AMP
           ENDIF
          ELSE
           FACT_AMP = 1.0
          ENDIF
        ELSE
          RI_TAU(N) = RI_FINAL(N)
          FACT_AMP = 1.0
        ENDIF   
    
        ENDDO  
        
   
        DO RK_STEP=1,3
          IF (NUM_PER_DIR.EQ.3) THEN
! c            IF (TIME_AD_METH.EQ.1) CALL RK_PER_1
! c            IF (TIME_AD_METH.EQ.2) CALL RK_PER_2            
          ELSEIF (NUM_PER_DIR.EQ.2) THEN
! c            IF (TIME_AD_METH.EQ.1) CALL RK_CHAN_1
! c            IF (TIME_AD_METH.EQ.2) CALL RK_CHAN_2            
! c            IF (TIME_AD_METH.EQ.3) CALL RK_CHAN_3            
          ELSEIF (NUM_PER_DIR.EQ.1) THEN
            IF (TIME_AD_METH.EQ.1) CALL RK_DUCT_1
            IF (TIME_AD_METH.EQ.2) CALL RK_DUCT_2
            IF (TIME_AD_METH.EQ.3) CALL RK_DUCT_3            
          ELSEIF (NUM_PER_DIR.EQ.0) THEN
            IF (TIME_AD_METH.EQ.1) CALL RK_CAV_1
            IF (TIME_AD_METH.EQ.2) CALL RK_CAV_2            
          END IF
        END DO
        TIME=TIME+DELTA_T
        FIRST_TIME=.FALSE.
! Save statistics to an output file

! c        IF (MOD(TIME_STEP,100).EQ.0) THEN
! c          IF (USE_MPI) THEN
! c           CALL GHOST_CHAN_MPI
! c          END IF
! c          CALL POISSON_P_DUCT
! c        ENDIF        


        IF (MOD(TIME_STEP,SAVE_STATS_INT).EQ.0) THEN
            CALL SAVE_STATS(.FALSE.)
        END IF
! Save the flow to a restart file
        IF (MOD(TIME_STEP,SAVE_FLOW_INT).EQ.0) THEN
          CALL SAVE_FLOW(.FALSE.)
        END IF
! Filter the scalar field
! c        DO N=1,N_TH
! c          IF (FILTER_TH(N)
! c     &       .AND.(MOD(TIME_STEP,FILTER_INT(N)).EQ.0)) THEN
! c          write(*,*) 'Filtering...'
! c          CALL FILTER(N)
! c          END IF 
! c        END DO
!        
! ! If we are considering a Near Wall Model, it may be necessary to
! ! filter the velocity field.  If so, use the following:

          IF (FILTER_VEL           &
             .AND.(MOD(TIME_STEP,FILTER_VEL_INT).EQ.0)) THEN
 
      if (rank == 0)  write(*,*) 'Filtering Velocity'
           call filter_all_real
! c           CALL APPLY_FILTER_VEL(CU1,1,NY)
! c           CALL APPLY_FILTER_VEL(CU2,1,NY)
! c           CALL APPLY_FILTER_VEL(CU3,1,NY)
      END IF
! 
! c        IF (MOD(TIME_STEP,200).EQ.0) THEN
! c          IF (USE_MPI) THEN
! c           CALL GHOST_CHAN_MPI
! c          END IF
! c          CALL POISSON_P_DUCT
! c         ENDIF
 
!      rktime = -(rktime - omp_get_wtime ( ))

      IF (rank .eq. 0) then
      write(6,*) 'Time taken in RK step', rktime

      open(99,file='out_screen.txt',form='formatted', &
      status='unknown', position='append')

      write(99,*) 'Time taken in RK step', rktime
      close(99)
      ENDIF

      END DO

! Calculate and display the runtime for the simulation
      CALL DATE_AND_TIME (VALUES=TIME_ARRAY)
      END_TIME=TIME_ARRAY(5)*3600+TIME_ARRAY(6)*60 &
        +TIME_ARRAY(7)+TIME_ARRAY(8)*0.001
      WRITE(*,*) 'Elapsed Time (sec): ',end_time-start_time
      WRITE(*,*) 'Seconds per Iteration: ' &
          ,(end_time-start_time)/N_TIME_STEPS

      TIME_STEP=TIME_STEP-1
      CALL SAVE_FLOW(.TRUE.)
      CALL SAVE_STATS(.TRUE.)
      WRITE(6,*)
      WRITE(6,*) '        ****** Hello world!  Have a nice day! ******'
      WRITE(6,*)
      END

! C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
SUBROUTINE INITIALIZE
!       INCLUDE 'header_duct'
use ntypes
use Domain
use Grid
use Fft_var
use TIME_STEP_VAR
use run_variable
use mg_vari, only : INIT_FLAG
use les_chan_var , only : C_DYN     
use mpi_var 
implicit none

      REAL    VERSION, CURRENT_VERSION, DAMP_FACT,local_tmp,RNUM1
      CHARACTER*35 file_1  
      logical RESET_TIME, TRIAL_MPI
      INTEGER I, J, K, N
     
      INTEGER, DIMENSION(:), ALLOCATABLE :: seed

! C Initialize the random number generator
      CALL RANDOM_SEED(SIZE = K)
      Allocate (seed(1:K))
      seed(1:K)=10
      CALL RANDOM_SEED(PUT = seed)

      IF (rank .eq. 0) then
      open(99,file='out_screen.txt',form='formatted', &
      status='unknown', position='append') 
 
      OPEN (11,file='input.dat',form='formatted',status='old')      

      WRITE(6,*) 'Note that this code is distributed under the ', &
                'GNU General Public License.'
      WRITE(6,*) 'No warranty is expressed or implied.'
      WRITE(6,*)
 
      WRITE(99,*) 'Note that this code is distributed under the ', &
                'GNU General Public License.'
      WRITE(99,*) 'No warranty is expressed or implied.'
      WRITE(99,*)

      endif 
!ALLOCATE VARIABLES

      call   allocate_var
     
! C Read input file.
! C   (Note - if you change the following section of code, update the
! C    CURRENT_VERSION number to make obsolete previous input files!)
      


      CURRENT_VERSION=0.9
      READ(11,*)
      READ(11,*)
      READ(11,*)
      READ(11,*)
      READ(11,*) VERSION
      IF (VERSION .NE. CURRENT_VERSION) STOP 'Wrong input data format. '
      READ(11,*)
      READ(11,*) USE_MPI
      READ(11,*)
      READ(11,*) NU, LX, LY, LZ
      READ(11,*)
      READ(11,*) N_TIME_STEPS, DELTA_T, RESET_TIME, VARIABLE_DT, CFL &
                 ,DFN, UPDATE_DT
      READ(11,*)
      READ(11,*) NUM_PER_DIR, TIME_AD_METH, LES, LES_MODEL_TYPE, LES_MODEL_TYPE_TH
      READ(11,*)
      READ(11,*) VERBOSITY, SAVE_FLOW_INT, SAVE_STATS_INT, MOVIE
      READ(11,*)
      READ(11,*) CREATE_NEW_FLOW, IC_TYPE, KICK, INT_TREAT
      READ(11,*)
      READ(11,*) F_TYPE, UBULK0, PX0, OMEGA0, AMP_OMEGA0, ANG_BETA,f_0
      IF (NUM_PER_DIR.eq.3) THEN
      READ(11,*)
      READ(11,*) U_BC_XMIN, U_BC_XMIN_C1, U_BC_XMIN_C2, U_BC_XMIN_C3
      READ(11,*) 
      READ(11,*) V_BC_XMIN, V_BC_XMIN_C1, V_BC_XMIN_C2, V_BC_XMIN_C3
      READ(11,*)
      READ(11,*) W_BC_XMIN, W_BC_XMIN_C1, W_BC_XMIN_C2, W_BC_XMIN_C3
      READ(11,*)
      READ(11,*) U_BC_XMAX, U_BC_XMAX_C1, U_BC_XMAX_C2, U_BC_XMAX_C3
      READ(11,*)
      READ(11,*) V_BC_XMAX, V_BC_XMAX_C1, V_BC_XMAX_C2, V_BC_XMAX_C3
      READ(11,*)
      READ(11,*) W_BC_XMAX, W_BC_XMAX_C1, W_BC_XMAX_C2, W_BC_XMAX_C3
      END IF
      IF (NUM_PER_DIR.gt.0) THEN
      READ(11,*)
      READ(11,*) U_BC_YMIN, U_BC_YMIN_C1, U_BC_YMIN_C2, U_BC_YMIN_C3
      READ(11,*) 
      READ(11,*) V_BC_YMIN, V_BC_YMIN_C1, V_BC_YMIN_C2, V_BC_YMIN_C3
      READ(11,*)
      READ(11,*) W_BC_YMIN, W_BC_YMIN_C1, W_BC_YMIN_C2, W_BC_YMIN_C3
      READ(11,*)
      READ(11,*) U_BC_YMAX, U_BC_YMAX_C1, U_BC_YMAX_C2, U_BC_YMAX_C3
      READ(11,*)
      READ(11,*) V_BC_YMAX, V_BC_YMAX_C1, V_BC_YMAX_C2, V_BC_YMAX_C3
      READ(11,*)
      READ(11,*) W_BC_YMAX, W_BC_YMAX_C1, W_BC_YMAX_C2, W_BC_YMAX_C3
      END IF
      IF (NUM_PER_DIR.lt.2) THEN
      READ(11,*)
      READ(11,*) U_BC_ZMIN, U_BC_ZMIN_C1, U_BC_ZMIN_C2, U_BC_ZMIN_C3
      READ(11,*) 
      READ(11,*) V_BC_ZMIN, V_BC_ZMIN_C1, V_BC_ZMIN_C2, V_BC_ZMIN_C3
      READ(11,*)
      READ(11,*) W_BC_ZMIN, W_BC_ZMIN_C1, W_BC_ZMIN_C2, W_BC_ZMIN_C3
      READ(11,*)
      READ(11,*) U_BC_ZMAX, U_BC_ZMAX_C1, U_BC_ZMAX_C2, U_BC_ZMAX_C3
      READ(11,*)
      READ(11,*) V_BC_ZMAX, V_BC_ZMAX_C1, V_BC_ZMAX_C2, V_BC_ZMAX_C3
      READ(11,*)
      READ(11,*) W_BC_ZMAX, W_BC_ZMAX_C1, W_BC_ZMAX_C2, W_BC_ZMAX_C3
      END IF
! Read Stochastic forcing parameters
      READ(11,*)
      READ(11,*) STOCHASTIC_FORCING
! Filter for the velocity field
      READ(11,*)
      READ(11,*) FILTER_VEL, FILTER_VEL_INT, FILTER_TYPE
      
      READ(11,*)
! Read in the parameters for the N_TH scalars
      DO N=1,N_TH
        READ(11,*)
        READ(11,*) CREATE_NEW_TH(N)
        READ(11,*)
        READ(11,*) FILTER_TH(N), FILTER_INT(N)
        READ(11,*)
        READ(11,*) RI_FINAL(N), PR(N), BACKGROUND_GRAD(N), CONT_STRAT
        READ(11,*)
        READ(11,*) TH_BC_YMIN(N),TH_BC_YMIN_C1(N),TH_BC_YMIN_C2(N) &
                  ,TH_BC_YMIN_C3(N)
        READ(11,*)
        READ(11,*) TH_BC_YMAX(N),TH_BC_YMAX_C1(N),TH_BC_YMAX_C2(N) &
                  ,TH_BC_YMAX_C3(N)
        IF (NUM_PER_DIR.lt.2) THEN
         READ(11,*)
         READ(11,*) TH_BC_ZMIN(N),TH_BC_ZMIN_C1(N),TH_BC_ZMIN_C2(N) &
                  ,TH_BC_ZMIN_C3(N)
         READ(11,*)
         READ(11,*) TH_BC_ZMAX(N),TH_BC_ZMAX_C1(N),TH_BC_ZMAX_C2(N) &
                  ,TH_BC_ZMAX_C3(N)
        ENDIF
!        READ(11,*)
      END DO
                            
      IF (N_TH .eq. 2) THEN
        READ(11,*)
        READ(11,*)  Non_linear_ST
      ENDIF

! C If we are using MPI, then Initialize the MPI Variables
      IF (USE_MPI) THEN
        CALL INT_MPI
      END IF

      TRIAL_MPI = .FALSE. 

      DO N=1,N_TH
       RI_TAU(N) = RI_FINAL(N)
      ENDDO
  
!!!!!!!!!!!!!! Some more input parameter !!!!!!!!!!!!!!!!!!!!!!
      ANG_BETA = ANG_BETA*3.14159265/180.0 
      Q_H0 = 1.d0
      H0    = 10.0*LY     
      In_H0 = 1.d0/H0 
      count_data = 0
      Gravity_g = 10.0d0 
      rho_0 = 1000.0d0
      theta_0 = 2.3d0 ;
      Sal_0 = 35.0d0 ; 
      dSaldz = 0.0d0 ; ! 100.0d0 ;

      MELTING_MODEL = .true.
      m_FP   = -0.060d0 ; 
      L_heat = 3.35*10.0d0**5.0
      C_sp_heat = 4.184*10.0d0**3.0


      alpha_w = 6.0d0*10.0d0**(-5)
      gamma_w = 8.0d0*10.0d0**(-4)

      Ratio_gr = Gravity_g/rho_0 ;
      Ratio_gr_a = Gravity_g*alpha_w
      Ratio_gr_g = Gravity_g*gamma_w

      DELTA_T_in = DELTA_T 
! C Initialize grid
      IF (rank .eq. 0) then

      WRITE(6,*) 'Grid size: NX =',NX,', NY =',NY,', NZ =',NZ,'.'
      WRITE(6,*) 'Domain size: LX =',LX,', LY =',LY,', LZ =',LZ,'.' 
      WRITE(6,*)'NU',NU,'DELTA_T',DELTA_T,'Delta_s',sqrt(2.0*NU/OMEGA0) 
      WRITE(6,*)'CFL',CFL, 'Diffu No' , DFN
      WRITE(6,*) 'Variable delta_t applited',  VARIABLE_DT,U_BC_ZMAX_C1

      WRITE(6,*) 'Inclined angle beta', ANG_BETA
! C      WRITE(6,*) 'Height of slope', H0 
       
      WRITE(99,*) 'Grid size: NX =',NX,', NY =',NY,', NZ =',NZ,'.'
      WRITE(99,*) 'Actual Grid size: NXV =',NXV,', NY =',NY,', NZV =',NZV,'.'
      WRITE(99,*) 'Actual Grid Fourier directions: NKXP = ', NKXP, ', NX2P = ',NX2P   
      WRITE(99,*) 'Domain size: LX =',LX,', LY =',LY,', LZ =',LZ,'.'
      WRITE(99,*)'NU',NU,'DELTA_T',DELTA_T,'Delta_s',sqrt(2.0*NU/OMEGA0)
      WRITE(99,*) 'Amp of tide',AMP_OMEGA0,'Frequency of tide',OMEGA0, 'Frequency of Rot',F_0
      WRITE(99,*)'CFL',CFL, 'Diffu No' , DFN
      WRITE(99,*) 'Variable delta_t applited',  VARIABLE_DT,U_BC_ZMAX_C1 

      IF (LES)THEN
      WRITE(99,*) 'Model type for les vel is', LES_MODEL_TYPE, 'and for temperature is ', &
      LES_MODEL_TYPE_TH  
      ENDIF

      DO N=1,N_TH
        WRITE(6,*) 'Scalar number: ',N
        WRITE(6,*) ' Final Richardson number: ',RI_FINAL(N)
        WRITE(6,*) '  Prandlt number: ',PR(N)
        WRITE(6,*) '  Ratio r/rho_0: ', Ratio_gr 
        WRITE(6,*) ' Gravity*gamma_a : ', Ratio_gr_a
        WRITE(6,*) ' Gravity*gamma_w : ', Ratio_gr_g       
        WRITE(6,*) ' Specific heat : ', C_sp_heat
        WRITE(6,*) ' Latent heat : ', L_heat
        WRITE(6,*) ' Factor freezing point : ', m_FP
        WRITE(6,*) ' Magnification of gradient', (PR(2)/PR(1))*(C_sp_heat/L_heat)*theta_0/m_FP 


        WRITE(99,*) 'Scalar number: ',N
        WRITE(99,*) ' Final Richardson number: ',RI_FINAL(N)
        WRITE(99,*) ' Prandlt number: ',PR(N)
        WRITE(99,*) ' Ratio(r/rho_0): ', Ratio_gr
        WRITE(99,*) ' Gravity*gamma_a : ', Ratio_gr_a
        WRITE(99,*) ' Gravity*gamma_w : ', Ratio_gr_g
      END DO
      close(99)
      endif
   
      NXM=NX-1
      NYM=NY-1
      NZM=NZ-1
      NXP_L=NX-(NXP+1)*rank-1
      NX2P_L=NKX+1-(NX2P+1)*rank-1
      write(6,*)'Lower bound rank', rank, min(NXP_L,NXP),min(NX2P_L,NX2P)

      IF (NUM_PER_DIR .GT. 0) THEN
         WRITE (6,*) 'Fourier in X'
         DO I=0,NX
           GX(I)=(I*LX)/NX
           DX(I)=LX/NX
           IF (VERBOSITY .GT. 3) WRITE(6,*) 'GX(',I,') = ',GX(I)
         END DO
      ELSE
         WRITE (6,*) 'Finite-difference in X'
         OPEN (30,file='xgrid.txt',form='formatted',status='old')
         READ (30,*) NX_T
! C Check to make sure that grid file is the correct dimensions
         IF (NX_T.ne.NX) THEN
           WRITE(6,*) 'NX, NX_T',NX,NX_T
           STOP 'Error: xgrid.txt wrong dimensions'
         END IF
         DO I=1,NX+1
           READ(30,*) GX(I)
           IF (VERBOSITY .GT. 3) WRITE(6,*) 'GX(',I,') = ',GX(I)
         END DO 
         DO I=1,NX
           READ(30,*) GXF(I)
           IF (VERBOSITY .GT. 3) WRITE(6,*) 'GXF(',I,') = ',GXF(I)
         END DO
! C Define ghost cells, if needed for this grid...
         GXF(0)=2.d0*GXF(1)-GXF(2)
         GXF(NX+1)=2.d0*GXF(NX)-GXF(NXM)
         GX(0)=2.d0*GX(1)-GX(2)
! C Define the grid spacings 
         DO I=1,NX+1
           DX(I)=(GXF(I)-GXF(I-1))
         END DO
         DO I=1,NX
           DXF(I)=(GX(I+1)-GX(I))
         END DO
         CLOSE(30)
      END IF

      IF (IC_TYPE.EQ.7)THEN
      ELSE

      IF (NUM_PER_DIR .GT. 1) THEN
         WRITE (6,*) 'Fourier in Z'
         DO K=0,NZ
           GZ(K)=(K*LZ)/NZ
           DZ(K)=LZ/NZ
           IF (VERBOSITY .GT. 3) WRITE(6,*) 'GZ(',K,') = ',GZ(K)
         END DO
      ELSE
         WRITE (6,*) 'Finite-difference in Z'
         OPEN (30,file='zgrid.txt',form='formatted',status='old')
         READ (30,*) NZ_T
! C Check to make sure that grid file is the correct dimensions
         IF (NZ_T.ne.NZ) THEN
           WRITE(6,*) 'NZ, NZ_T',NZ,NZ_T
           STOP 'Error: zgrid.txt wrong dimensions'
         END IF
         DO K=1,NZ+1
           READ(30,*) GZ(k)
           IF (VERBOSITY .GT. 3) WRITE(6,*) 'GZ(',K,') = ',GZ(K)
         END DO 
         DO K=1,NZ
           READ(30,*) GZF(k)
           IF (VERBOSITY .GT. 3) WRITE(6,*) 'GZF(',K,') = ',GZF(K)
         END DO
! C Define ghost cells, if needed for this grid...
         GZF(0)=2.d0*GZF(1)-GZF(2)
         GZF(NZ+1)=2.d0*GZF(NZ)-GZF(NZM)
         GZ(0)=2.d0*GZ(1)-GZ(2)
! C Define grid spacing 
         DO K=1,NZ+1
           DZ(K)=(GZF(K)-GZF(K-1))
         END DO
         DO K=1,NZ
           DZF(K)=(GZ(K+1)-GZ(K))
           write(400,*) DZF(K)
         END DO
         DZ(0)=DZ(1)
         DZF(NZ+1)=DZF(NZ) 
         CLOSE(30)
      END IF
      
      IF (NUM_PER_DIR .GT. 2) THEN
         WRITE (6,*) 'Fourier in Y'
         DO J=0,NY
           GY(J)=(J*LY)/NY
           DY(J)=LY/NY
           IF (VERBOSITY .GT. 3) WRITE(6,*) 'GY(',J,') = ',GY(J)
         END DO
      ELSE
         WRITE (6,*) 'Finite-difference in Y'
         OPEN (30,file='./ygrid.txt',form='formatted',status='old')
         READ (30,*) NY_T
! C Check to make sure that grid file is the correct dimensions
         IF (NY_T.ne.NY) THEN
           WRITE(6,*) 'NY, NY_T',NY,NY_T
           STOP 'Error: ygrid.txt wrong dimensions'
         END IF
         DO J=1,NY+1
           READ(30,*) GY(j)
           IF (VERBOSITY .GT. 3) WRITE(6,*) 'GY(',J,') = ',GY(J)
         END DO
         DO J=1,NY
           READ(30,*) GYF(j)
           IF (VERBOSITY .GT. 3) WRITE(6,*) 'GYF(',J,') = ',GYF(J)
         END DO
         CLOSE(30)
 
! C Define ghost cells
           GYF(0)=2.d0*GYF(1)-GYF(2)
           GYF(NY+1)=2.d0*GYF(NY)-GYF(NYM)
           GY(0)=2.d0*GY(1)-GY(2)

! C Define grid spacing
         DO J=1,NY+1
           DY(J)=(GYF(J)-GYF(J-1))
         END DO
         DO J=1,NY
           DYF(J)=(GY(J+1)-GY(J))
         END DO
         DY(0)=DY(1)
         DYF(NY+1)=DYF(NY)
      END IF
      ENDIF



 
      DO J=0,NY+2
       DO K=0,NZ+2
        U2_BAR(K,J)=0.d0
        U1_BAR(J)=0.d0
        U3_BAR(K,J)=0.d0
       ENDDO 
      ENDDO
! C   Boundary condition prescribed by the Blasius for IC_TYPE = 5
! c      IF ( IC_TYPE .EQ. 5) THEN
! c      OPEN (80,file='./int_vel_xy.txt',form='formatted',status='old')     
! c       DO J=1,NY+2
! c        DO K=1,NZ+2
! c        read(80,*) U3_BAR(K,J), U2_BAR(K,J)
! c       ENDDO
! c       ENDDO
! 
! c       DO J=1,NY+1 
! c        DO K=2,NZ+1
! c          S1(0,K,J) = 0.5*(U3_BAR(K,J) + U3_BAR(K-1,J))
! c        ENDDO
! c       ENDDO
! 
! c       DO J=1,NY
! c        DO K=2,NZ
! c          U3_BAR(K,J) =S1(0,K,J)
! c        ENDDO
! c       ENDDO 
!        
! c       DO K=1,NZ
! c          U3_BAR(K,NY+1) =2.0*U3_BAR(K,NY)-U3_BAR(K,NY-1)
! c	  U2_BAR(K,NY+1) =2.0*U2_BAR(K,NY)-U2_BAR(K,NY-1)
! c       ENDDO
       
       DO J=0,NY+1
        DO K=0,NY+1
           U3_BAR(K,J) = 0.0
	   U2_BAR(K,J) = 0.0
	ENDDO   
        ENDDO

!       CALL VEL_PROFILE              
        
        call allocation_ub (.true.)

        IF (LES) THEN
         call allocation_les_var
        ENDIF

        IF ( IC_TYPE .EQ. 7) THEN 
         CALL   JACOBI_TRANS
	 CALL VEL_PROFILE
	ENDIF 
! c      ENDIF	 
	  


     call allocation_u (.true.)
     call allocation_v (.true.)
     call allocation_w (.true.)  
     call allocation_p (.true.)

     call allocation_R1 (.true.) 
     call allocation_R2 (.true.)
     call allocation_R3 (.true.) 
     call allocation_F1 (.true.)
     call allocation_F2 (.true.)
     call allocation_F3 (.true.)
 


      DO N=1,N_TH
       call allocation_th (.true.)
       call allocation_Rth (.true.)
       call allocation_Fth (.true.)
      ENDDO

! C Initialize storage arrays.
!      DO K=0,NZ+1
!        DO I=0,NX+1 
!          DO J=0,NY+1
!            U1(I,K,J)=0.
!            U3(I,K,J)=0.
!            U2(I,K,J)=0.
!            P (I,K,J)=0.
!            R1(I,K,J)=0.
!            R2(I,K,J)=0.
!            R3(I,K,J)=0.
!            F1(I,K,J)=0.
!            F2(I,K,J)=0.
!            F3(I,K,J)=0.
! Array for LES subgrid model
!            IF (LES) THEN
!             NU_T(I,K,J)=0.
!             DO N=1,N_TH
!               KAPPA_T(I,K,J,N)=0.
!             END DO
!            ENDIF
!          END DO
!        END DO
!      END DO


! C Initialize storage arrays.
      DO K=0,NZV-1
        DO I=0,NXP
          DO J=0,NY+1
            U1X(I,K,J)=0.
            U3X(I,K,J)=0.
            U2X(I,K,J)=0.
            PX (I,K,J)=0.
            R1X(I,K,J)=0.
            R2X(I,K,J)=0.
            R3X(I,K,J)=0.
            F1X(I,K,J)=0.
            F2X(I,K,J)=0.
            F3X(I,K,J)=0.
            DO N=1,N_TH
             THX(I,K,J,N)=0.
             RTHX(I,K,J,N)=0. 
             FTHX(I,K,J,N)=0.
            END DO
            
! Array for LES subgrid model
            IF (LES) THEN
             NU_T(I,K,J)=0.
             DO N=1,N_TH
               KAPPA_T(I,K,J,N)=0.
             END DO
            ENDIF
          END DO
        END DO
      END DO


      IF (LES) THEN
       DO K=0,NZ+1
        DO J=0,NY+1
          C_DYN(K,J)=0.d0
        END DO
       END DO
      ENDIF

      write(6,*) 'I am here 1'
! C Initialize FFT package (includes defining the wavenumber vectors).

!$OMP PARALLEL
      CALL INIT_FFT
!$OMP END PARALLEL

 
! C Initialize RKW3 parameters.
      H_BAR(1)=DELTA_T*(8.0/15.0)
      H_BAR(2)=DELTA_T*(2.0/15.0)
      H_BAR(3)=DELTA_T*(5.0/15.0)
      BETA_BAR(1)=1.0
      BETA_BAR(2)=25.0/8.0
      BETA_BAR(3)=9.0/4.0
      ZETA_BAR(1)=0.0
      ZETA_BAR(2)=-17.0/8.0
      ZETA_BAR(3)=-5.0/4.0
      
! C Initialize case-specific packages.
      IF (NUM_PER_DIR.EQ.3) THEN
! C        CALL INIT_PER
      ELSEIF (NUM_PER_DIR.EQ.2) THEN
          CALL INIT_CHAN
      ELSEIF (NUM_PER_DIR.EQ.1) THEN
        CALL INIT_DUCT
      ELSEIF (NUM_PER_DIR.EQ.0) THEN
        CALL INIT_CAV
      END IF

!ccccccccccccccccccccccccccccccccccccccccccccccccc
call open_MP_initialization
!ccccccccccccccccccccccccccccccccccccccccccccccccc

      IF(TRIAL_MPI) THEN  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   Trial section for MPI 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ENDIF 
       

! C Initialize values for reading of scalars
!      NUM_READ_TH=0
!      DO N=1,N_TH
!        IF (CREATE_NEW_TH(N)) THEN
!          NUM_READ_TH=NUM_READ_TH 
!        ELSE
!          NUM_READ_TH=NUM_READ_TH + 1
!          READ_TH_INDEX(NUM_READ_TH)=N
!        END IF
      
                
!      IF (NUM_PER_DIR.EQ.1) THEN
!        CALL CREATE_TH_DUCT
!      ELSE IF (NUM_PER_DIR.EQ.2) THEN
!        CALL CREATE_TH_CHAN 
!      ELSE IF (NUM_PER_DIR.EQ.3) THEN
!         CALL CREATE_TH_PER
!      END IF 
!     
!      ENDDO

      IF (CREATE_NEW_TH(1)) THEN
       IF (NUM_PER_DIR.EQ.1) THEN
         CALL CREATE_TH_DUCT
       ELSE IF (NUM_PER_DIR.EQ.2) THEN
        CALL CREATE_TH_CHAN
       ELSE IF (NUM_PER_DIR.EQ.3) THEN
! C        CALL CREATE_TH_PER
       END IF
      ENDIF 

      
!      write(6,*) 'I am here 1'


      IF (CREATE_NEW_FLOW) THEN
        IF (NUM_PER_DIR.EQ.3) THEN
! C          CALL CREATE_FLOW_PER
        ELSEIF (NUM_PER_DIR.EQ.2) THEN
          CALL CREATE_FLOW_CHAN
        ELSEIF (NUM_PER_DIR.EQ.1) THEN
          CALL CREATE_FLOW_DUCT
        ELSEIF (NUM_PER_DIR.EQ.0) THEN
          CALL CREATE_FLOW_CAV
        END IF
        write(*,*) 'A new flowfield has been created'
! c        CALL SAVE_STATS(.FALSE.)
        CALL SAVE_FLOW(.FALSE.)
        write(*,*) 5000
      ELSE

! C     Convert velocity back to Fourier space
                CALL REAL_FOURIER_TRANS_U1 (.true.)
                CALL REAL_FOURIER_TRANS_U2 (.true.)
                CALL REAL_FOURIER_TRANS_U3 (.true.)
                CALL REAL_FOURIER_TRANS_P (.true.)
        
                CALL REAL_FOURIER_TRANS_R1 (.true.)
                CALL REAL_FOURIER_TRANS_R2 (.true.)
                CALL REAL_FOURIER_TRANS_R3 (.true.)
        
                call allocation_F1 (.false.)
                call allocation_F2 (.false.)
                call allocation_F3 (.false.)

        CU1X(:,:,:) = (0.d0,0.d0)
        CU2X(:,:,:) = (0.d0,0.d0)
        CU3X(:,:,:) = (0.d0,0.d0)
        CPX(:,:,:) = (0.d0,0.d0)

         CALL REAL_FOURIER_TRANS_TH (.true.)
         CALL REAL_FOURIER_TRANS_Rth (.true.)
         call allocation_Fth(.false.)

        
        DO N=1,N_TH
!                  CALL REAL_FOURIER_TRANS_TH (.true.)
!                  CALL REAL_FOURIER_TRANS_Rth (.true.)
!                  call allocation_Fth(.false.)

          CTHX(:,:,:,N) = (0.d0,0.d0)
        ENDDO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        if(rank.EQ.0) write(*,*) 'Reading flow...'
        CALL READ_FLOW
        if(rank.EQ.0) write(*,*) 'Done reading flow'

        IF ( INT_TREAT ) THEN

! TEMPORARY!! SCALE THE VELOCITY FLUCTIONS
       
        ENDIF

! C Initialize flow.
      IF (RESET_TIME .OR. CREATE_NEW_FLOW) THEN
        PREVIOUS_TIME_STEP=0
        TIME_STEP=0
        TIME=0
      END IF

        CALL SAVE_STATS(.FALSE.)
        INIT_FLAG = .TRUE.       
        
        IF (NUM_PER_DIR.EQ.3) THEN
! C          CALL POISSON_P_PER
        ELSEIF (NUM_PER_DIR.EQ.2) THEN
          CALL POISSON_P_CHAN
        ELSEIF (NUM_PER_DIR.EQ.1) THEN
! c          CALL POISSON_P_DUCT
        ELSEIF (NUM_PER_DIR.EQ.0) THEN
          CALL POISSON_P_CAV
        END IF

          CALL SAVE_STATS(.FALSE.)   

      END IF


      RETURN
      END





! C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE SAVE_STATS(FINAL)
! C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
use ntypes
use Domain
use Grid
use Fft_var
use TIME_STEP_VAR
use run_variable
use mpi_var, only: rank      
implicit none

      LOGICAL FINAL

      IF (NUM_PER_DIR.EQ.3) THEN
! C        CALL SAVE_STATS_PER(FINAL)          
      ELSEIF (NUM_PER_DIR.EQ.2) THEN
! C        CALL SAVE_STATS_CHAN(FINAL)          
      ELSEIF (NUM_PER_DIR.EQ.1) THEN
        IF ( IC_TYPE .EQ. 7) THEN
         CALL SAVE_STATS_CURVI(FINAL)
        ELSE
         CALL SAVE_STATS_DUCT(FINAL)
	ENDIF           
      ELSEIF (NUM_PER_DIR.EQ.0) THEN
        CALL SAVE_STATS_CAV(FINAL)          
      END IF

      if (rank ==0) then
      write(*,*) 'done save_stats diablo'
      endif

! c      IF (FINAL) THEN
! c        IF (NUM_PER_DIR.EQ.3) THEN
! C          CALL VIS_FLOW_PER         
! c        ELSEIF (NUM_PER_DIR.EQ.2) THEN
! !          CALL VIS_FLOW_CHAN         
! c        ELSEIF (NUM_PER_DIR.EQ.1) THEN
! c          CALL VIS_FLOW_DUCT          
! c        ELSEIF (NUM_PER_DIR.EQ.0) THEN
! c          CALL VIS_FLOW_CAV         
! c        END IF
! c      END IF

      RETURN
      END

! C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE READ_FLOW
! C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
!       INCLUDE 'header_duct'
use ntypes
use Domain
use Grid
use Fft_var
use TIME_STEP_VAR
use run_variable
use mpi_var, only: rank      
implicit none

      CHARACTER*35 FNAME
      CHARACTER*35 FNAME_TH(N_TH)
      INTEGER I, J, K, N, NUM_PER_DIR_T


      IF (USE_MPI) THEN
       I = RANK+1

        FNAME ='last_saved/diablo_'                 &
              //CHAR(MOD(I,100)/10+48)              &
              //CHAR(MOD(I,10)+48) // '.start'
        DO N=1,N_TH
           FNAME_TH(N)='last_saved/diablo_th'       & 
             //CHAR(MOD(N,100)/10+48)               & 
             //CHAR(MOD(N,10)+48) //'_'             & 
             //CHAR(MOD(I,100)/10+48)               &
             //CHAR(MOD(I,10)+48) // '.start'
        ENDDO
      ELSE
       FNAME='diablo.start'
       DO N=1,N_TH
        FNAME_TH(N)='diablo_th' &
             //CHAR(MOD(N,100)/10+48) &
             //CHAR(MOD(N,10)+48) // '.start'
       END DO
      ENDIF
      WRITE(6,*)   'Reading flow from velocity data ', FNAME
       WRITE(6,*)   'Reading flow from theta data ', FNAME_TH(1)

      OPEN(UNIT=10,FILE=FNAME,STATUS="OLD",FORM="UNFORMATTED")
      READ (10) NX_T, NY_T, NZ_T, NUM_PER_DIR_T, TIME, TIME_STEP


      write(*,*) 'NX_T, NY_T, NZ_T: ',NX_T,NY_T,NZ_T,NUM_PER_DIR_T, &
                  NUM_PER_DIR

      IF ((NX .NE. NX_T) .OR. (NY .NE. NY_T) .OR. (NZ .NE. NZ_T))THEN
           WRITE(*,*) NX_T, NY_T, NZ_T
           STOP 'Error: old flowfield wrong dimensions. '
      ENDIF
         
      IF (NUM_PER_DIR .NE. NUM_PER_DIR_T) &
         STOP 'Error: old flowfield wrong NUM_PER_DIR. '

      write(*,*) 'READING FLOW'
      IF (NUM_PER_DIR.EQ.3) THEN
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ELSEIF (NUM_PER_DIR.EQ.2) THEN
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ELSEIF (NUM_PER_DIR.EQ.1) THEN
        READ (10) (((CU1X(I,K,J),I=0,NX2P),K=0,NZ+1),J=0,NY+1), &
                 (((CU2X(I,K,J),I=0,NX2P),K=0,NZ+1),J=0,NY+1),  &
                 (((CU3X(I,K,J),I=0,NX2P),K=0,NZ+1),J=0,NY+1),  &
                 (((CPX(I,K,J),I=0,NX2P),K=0,NZ+1),J=0,NY+1) 
         Write(6,*) 'Done with velocity field'
        DO N=1,N_TH
! Specify in input.dat which scalars are to be read
          OPEN(UNIT=11,FILE=FNAME_TH(N),STATUS="OLD" &
                ,FORM="UNFORMATTED")
          READ (11) NX_T, NY_T, NZ_T, NUM_PER_DIR_T, TIME, TIME_STEP
          READ (11) (((CTHX(I,K,J,N) &
                ,I=0,NX2P),K=0,NZ+1),J=0,NY+1)
         CLOSE(11)
         Write(6,*) 'Done with theta field', N
        END DO

      ELSEIF (NUM_PER_DIR.EQ.0) THEN
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      END IF
      CLOSE(10)
      CLOSE(11)

! C Apply initial boundary conditions, set ghost cells
      IF (IC_TYPE .EQ. 7)THEN
      ELSE 
        call APPLY_BC_VEL_LOWER
        call APPLY_BC_VEL_UPPER
        call APPLY_BC_VEL_LEFT
        call APPLY_BC_VEL_RIGHT
      ENDIF

      RETURN
      END

! C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE SAVE_FLOW(FINAL)
! C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
use ntypes
use Domain
use Grid
use Fft_var
use TIME_STEP_VAR
use run_variable
use mpi_var, only : rank      
implicit none

      CHARACTER*35 FNAME, frame(4),frame_th(4)
      CHARACTER*35 FNAME_TH(N_TH)
      INTEGER      I, J, K, N, m, NY_min, NY_max,no_p,ny_p
      LOGICAL      FINAL
      PARAMETER  (no_p=4, ny_p=66)
      


      IF (USE_MPI) THEN
       I=1+RANK
       IF (FINAL) THEN
         FNAME='last_saved/diablo_'      &
             //CHAR(MOD(I,100)/10+48)    &
             //CHAR(MOD(I,10)+48) // '.res'


         DO N=1,N_TH
            FNAME_TH(N)='last_saved/diablo_th'  &
              //CHAR(MOD(N,100)/10+48)          &
              //CHAR(MOD(N,10)+48) //'_'        &
              //CHAR(MOD(I,100)/10+48)          &
              //CHAR(MOD(I,10)+48) // '.res'
         END DO

       ELSE
         FNAME='last_saved/diablo_'       &
              //CHAR(MOD(I,100)/10+48)    &
              //CHAR(MOD(I,10)+48) // '.saved'

         DO N=1,N_TH
            FNAME_TH(N)='last_saved/diablo_th'  &
              //CHAR(MOD(N,100)/10+48)          &
              //CHAR(MOD(N,10)+48) //'_'        &
              //CHAR(MOD(I,100)/10+48)          & 
              //CHAR(MOD(I,10)+48) // '.saved'
         END DO
        ENDIF   
      ELSE
       IF (FINAL) THEN
          FNAME='diablo.res'
          DO N=1,N_TH
            FNAME_TH(N)='diablo_th' &
              //CHAR(MOD(N,100)/10+48) &
              //CHAR(MOD(N,10)+48) // '.res'
          END DO
       ELSE
          FNAME='diablo.saved'
          DO N=1,N_TH
            FNAME_TH(N)='diablo_th' &
              //CHAR(MOD(N,100)/10+48) &
              //CHAR(MOD(N,10)+48) // '.saved'
          END DO

        write(*,*) 'bishakh 000000000'
       END IF
      ENDIF
      WRITE(6,*) 'Writing flow to ',FNAME

      OPEN(UNIT=10,FILE=FNAME,STATUS="UNKNOWN",FORM="UNFORMATTED")
      WRITE(10) NX, NY, NZ, NUM_PER_DIR, TIME, TIME_STEP


      IF (NUM_PER_DIR.EQ.3) THEN
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ELSEIF (NUM_PER_DIR.EQ.2) THEN
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ELSEIF (NUM_PER_DIR.EQ.1) THEN
        WRITE(10) (((CU1X(I,K,J),I=0,NX2P),K=0,NZ+1),J=0,NY+1), &
                 (((CU2X(I,K,J),I=0,NX2P),K=0,NZ+1),J=0,NY+1),  &
                 (((CU3X(I,K,J),I=0,NX2P),K=0,NZ+1),J=0,NY+1), &
                 (((CPX(I,K,J),I=0,NX2P),K=0,NZ+1),J=0,NY+1)
        DO N=1,N_TH
          OPEN(UNIT=11,FILE=FNAME_TH(N),STATUS="UNKNOWN" &
            ,FORM="UNFORMATTED")
         WRITE(11) NX, NY, NZ, NUM_PER_DIR, TIME, TIME_STEP
         WRITE(11) (((CTHX(I,K,J,N),I=0,NX2P),K=0,NZ+1),J=0,NY+1)
         CLOSE(11)
        END DO

      ELSEIF (NUM_PER_DIR.EQ.0) THEN
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!       
      END IF
      CLOSE(10)
      CLOSE(11)

      RETURN
      END


! C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE FILTER(n)
! C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
use TIME_STEP_VAR

      integer n

      IF (NUM_PER_DIR.EQ.3) THEN
! C        CALL FILTER_PER
      ELSEIF (NUM_PER_DIR.EQ.2) THEN
! C        CALL FILTER_CHAN(n)
      END IF

      RETURN
      END

! C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE JACOBI_TRANS
! C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
use ntypes
use Domain
use Grid
use Fft_var
use TIME_STEP_VAR
use run_variable
use les_chan_var
use mpi_var, only : rank     
implicit none

      INTEGER  I,J,K,N
     
!       REAL*8   x_zeta(NZ+1,NY+1,2), y_zeta(NZ+1,NY+1,2), &
!               x_eta(NZ+1,NY+1,2), y_eta(NZ+1,NY+1,2),  &
!               zeta_x(NZ+1,NY+1,2), zeta_y(NZ+1,NY+1,2), &
!               eta_x(NZ+1,NY+1,2), eta_y(NZ+1,NY+1,2),  &
!               IN_JACO(NZ+1,NY+1,2),  &
!               xpoint(0:NZ+1,0:NY+1), ypoint(0:NZ+1,0:NY+1) 



      real*8  TH_TOP,b_min,b_max,co_b
      integer NZ_C,NY_C,NZ_min,NZ_max,NY_min,NY_max

      real(r8),allocatable,dimension(:,:,:) :: x_zeta, y_zeta, &
         x_eta, y_eta, zeta_x, zeta_y, eta_x, eta_y, IN_JACO
!      real(r8),allocatable,dimension(:,:)   :: xpoint, ypoint

      allocate (x_zeta(NZ+1,NY+1,2), y_zeta(NZ+1,NY+1,2))
      allocate (x_eta(NZ+1,NY+1,2), y_eta(NZ+1,NY+1,2))
      allocate (zeta_x(NZ+1,NY+1,2), zeta_y(NZ+1,NY+1,2))
      allocate (eta_x(NZ+1,NY+1,2), eta_y(NZ+1,NY+1,2))
      allocate (IN_JACO(NZ+1,NY+1,2))
!      allocate (xpoint(0:NZ+1,0:NY+1), ypoint(0:NZ+1,0:NY+1))
      


!      open(31,file='density_data_input.dat',form='formatted', &
!      status='old')
  
 
      write(6,*) 'Jacobian tranformation is starting'
      open(21,file='GRID_XZ.dat',form='formatted',status='old') 
      xpoint(:,:) =0.0d0 
      ypoint(:,:) =0.0d0 
      DO J=0,NY+1
       DO K=0,NZ+1
         READ(21,*)xpoint(K,J),ypoint(K,J)
         ypoint(K,J) = 1.8*ypoint(K,J)
       ENDDO
      ENDDO
      close(21)
       
111   format(2f16.8)

     DO N=1,N_TH 
      DO J=0,NY+1
        DO K=0,NZ+1
          IF (CONT_STRAT) THEN
             IF (N .eq. 1) THEN
!             background condition for temperature 
              TH_BAR(K,J,N)= theta_0 !(-ypoint(k,j)+ypoint(0,NY+1))
             ELSEIF (N .eq. 2) THEN
              TH_BAR(K,J,N)= Sal_0 + dSaldz*(-ypoint(k,j)+ypoint(0,NY+1))
             ENDIF
          ELSE
             READ(31,*)TH_BAR(K,J,N)
          ENDIF
        ENDDO
       ENDDO
      ENDDO

     close(31)


     IF (N_TH .eq. 2) THEN
     IF (Non_linear_ST) THEN
          call density_TC(.true.,.false.)
     ELSE
      DO J=0,NY+1
      DO K=0,NZ+1
        TH_BAR(K,J,3)= 1.0d0 + gamma_w*(TH_BAR(K,J,2)-Sal_0) - alpha_w*(TH_BAR(K,J,1)-theta_0)
      ENDDO
      ENDDO
     ENDIF
     ENDIF
 


      
      DO J=1,NY
       DO K=1,NZ
         x_zeta(K,J,1) =xpoint(K+1,J)-xpoint(K,J)
         x_zeta(K,J,2) =0.250d0*(xpoint(K+1,J+1)+xpoint(K+1,J) &
                           - xpoint(K-1,J+1)- xpoint(K-1,J)) 

         y_zeta(K,J,1) = ypoint(K+1,J)-ypoint(K,J)
         y_zeta(K,J,2) = 0.250d0*(ypoint(K+1,J+1)+ypoint(K+1,J) &
                           - ypoint(K-1,J+1) - ypoint(K-1,J))


         x_eta(K,J,1) = 0.250d0*(xpoint(K,J+1)+xpoint(K+1,J+1) &
                           - xpoint(K,J-1)- xpoint(K+1,J-1))
         x_eta(K,J,2) = xpoint(K,J+1)-xpoint(K,J)

         y_eta(K,J,1) = 0.250d0*(ypoint(K,J+1)+ypoint(K+1,J+1)  &
                           - ypoint(K,J-1)- ypoint(K+1,J-1))
         y_eta(K,J,2) = (ypoint(K,J+1)- ypoint(K,J)) 
       ENDDO
      ENDDO


      DO K=1,NZ+1
       x_eta(K,1,1)  = 2.0d0*x_eta(K,2,1)- x_eta(K,3,1)
       DO I=1,2
        x_eta(K,NY+1,I)  = 2.0d0*x_eta(K,NY,I) - x_eta(K,NY-1,I)
	x_zeta(K,NY+1,I) = 2.0d0*x_zeta(K,NY,I)- x_zeta(K,NY-1,I) 
       ENDDO
       y_eta(K,1,1)  = 2.0d0*y_eta(K,2,1)- y_eta(K,3,1)
       DO I=1,2
        y_eta(K,NY+1,I)  = 2.0d0*y_eta(K,NY,I)- y_eta(K,NY-1,I)
	y_zeta(K,NY+1,I) = 2.0d0*y_zeta(K,NY,I)- y_zeta(K,NY-1,I)
       ENDDO   
      ENDDO  

      DO J=1,NY+1
        x_zeta(1,J,2) = 2.0d0*x_zeta(2,J,2)- x_zeta(3,J,2)
	DO I=1,2
         x_zeta(NZ+1,J,I) = 2.0d0*x_zeta(NZ,J,I)- x_zeta(NZ-1,J,I)
	 x_eta(NZ+1,J,I) = 2.0d0*x_eta(NZ,J,I)- x_eta(NZ-1,J,I)
        ENDDO
        y_zeta(1,J,2)    = 2.0d0*y_zeta(2,J,2) - y_zeta(3,J,2)
	DO I=1,2
         y_zeta(NZ+1,J,I) = 2.0d0*y_zeta(NZ,J,I) - y_zeta(NZ-1,J,I)
	 y_eta(NZ+1,J,I) = 2.0d0*y_eta(NZ,J,I) - y_eta(NZ-1,J,I)
	ENDDO 	 
      ENDDO

! C make smooth transformation
      DO J=1,NY+1
       DO K=1,NZ+1
         DO I=1,2
	   x_zeta(K,J,I) = dble(AINT(x_zeta(K,J,I)*10**8)) &
     	/dble(10**8)
	   x_eta(K,J,I)  = dble(AINT(x_eta(K,J,I)*10**8)) &
     	/dble(10**8)
	   
	   y_zeta(K,J,I) = dble(AINT(y_zeta(K,J,I)*10**8)) &
     	/dble(10**8)
	   y_eta(K,J,I)  = dble(AINT(y_eta(K,J,I)*10**8)) &
     	/dble(10**8)
	 ENDDO
       ENDDO 	  
      ENDDO 	  
      
      
      DO J=1,NY+1
       DO K=1,NZ+1
        DO I = 1,2
         IN_JACO(K,J,I) = x_zeta(K,J,I)* y_eta(K,J,I) &
                    - x_eta(K,J,I)*y_zeta(K,J,I)
  

         zeta_x(K,J,I) =   y_eta(K,J,I)/IN_JACO(K,J,I)  
         eta_x(K,J,I)  =  -y_zeta(K,J,I)/IN_JACO(K,J,I) 
 
         zeta_y(K,J,I) =  -x_eta(K,J,I)/IN_JACO(K,J,I) 
         eta_y(K,J,I)  =   x_zeta(K,J,I)/IN_JACO(K,J,I)
        ENDDO

        DO I = 1,2
         GMAT_11(K,J,I)   =  (zeta_x(K,J,I)*zeta_x(K,J,I) &
                     +  zeta_y(K,J,I)*zeta_y(K,J,I))*IN_JACO(K,J,I)   
         GMAT_12(K,J,I)   =  (zeta_x(K,J,I)*eta_x(K,J,I) &
                     +  zeta_y(K,J,I)*eta_y(K,J,I))*IN_JACO(K,J,I)
         GMAT_22(K,J,I)   =  (eta_x(K,J,I)*eta_x(K,J,I) &
                     +  eta_y(K,J,I)*eta_y(K,J,I))*IN_JACO(K,J,I) 

         CJOB_11(K,J,I)   = zeta_x(K,J,I)*IN_JACO(K,J,I)  
         CJOB_12(K,J,I)   = zeta_y(K,J,I)*IN_JACO(K,J,I)
         CJOB_21(K,J,I)   = eta_x(K,J,I)*IN_JACO(K,J,I)
         CJOB_22(K,J,I)   = eta_y(K,J,I)*IN_JACO(K,J,I)  
        ENDDO

       ENDDO
      ENDDO
     
      
      
      DO J=1,NY
       DO K=1,NZ
         INT_JACOB(K,J)=0.250d0*(IN_JACO(K,J,1)+ IN_JACO(K,J,2) &
                  + IN_JACO(K,J+1,1) + IN_JACO(K+1,J,2))
       ENDDO
      ENDDO
      
      DO J=0,NY+1
! c       INT_JACOB(NZ,J)   = 2.0*INT_JACOB(NZ-1,J) - INT_JACOB(NZ-2,J)
       INT_JACOB(NZ+1,J) = 2.0d0*INT_JACOB(NZ,J) - INT_JACOB(NZ-1,J)
       INT_JACOB(0,J) = 2.0d0*INT_JACOB(1,J) - INT_JACOB(2,J)
      ENDDO
      
      DO K=0,NZ+1
       INT_JACOB(K,NY+1) = 2.0d0*INT_JACOB(K,NY) - INT_JACOB(K,NY-1)
       INT_JACOB(K,0) = 2.0d0*INT_JACOB(K,1) - INT_JACOB(K,2)
      ENDDO 
      
! c    Correction for 0.0
      INT_JACOB(0,0)  = 0.5*( 2.0d0*INT_JACOB(1,0) - INT_JACOB(2,0) &
            +  2.0d0*INT_JACOB(0,1) - INT_JACOB(0,2) ) 
     
      INT_JACOB(NZ+1,NY+1)  = 0.5*( 2.0d0*INT_JACOB(NZ,NY+1)  &
            - INT_JACOB(NZ-1,NY+1) +  2.0d0*INT_JACOB(NZ+1,NY) &
            - INT_JACOB(NZ+1,NY-1) )   
     
      DO J=0,NY+1
       DO I=1,2
        CJOB_11(0,J,I) = 2.d0*CJOB_11(1,J,I) - CJOB_11(2,J,I)
        CJOB_12(0,J,I) = 2.d0*CJOB_12(1,J,I) - CJOB_12(2,J,I)
        CJOB_21(0,J,I) = 2.d0*CJOB_21(1,J,I) - CJOB_21(2,J,I)
        CJOB_22(0,J,I) = 2.d0*CJOB_22(1,J,I) - CJOB_22(2,J,I)
	
	GMAT_11(0,J,I) = 2.d0*GMAT_11(1,J,I) - GMAT_11(2,J,I)
	GMAT_12(0,J,I) = 2.d0*GMAT_12(1,J,I) - GMAT_12(2,J,I)
	GMAT_22(0,J,I) = 2.d0*GMAT_22(1,J,I) - GMAT_22(2,J,I)
       ENDDO
      ENDDO
      
      DO K=0,NZ+1
       DO I=1,2
        CJOB_11(K,0,I) = 2.d0*CJOB_11(K,1,I) - CJOB_11(K,2,I)
        CJOB_12(K,0,I) = 2.d0*CJOB_12(K,1,I) - CJOB_12(K,2,I)
	CJOB_21(K,0,I) = 2.d0*CJOB_21(K,1,I) - CJOB_21(K,2,I)
	CJOB_22(K,0,I) = 2.d0*CJOB_22(K,1,I) - CJOB_22(K,2,I)
	
	GMAT_11(K,0,I) = 2.0d0*GMAT_11(K,1,I) - GMAT_11(K,2,I)
	GMAT_12(K,0,I) = 2.0d0*GMAT_12(K,1,I) - GMAT_12(K,2,I) 
	GMAT_22(K,0,I) = 2.0d0*GMAT_22(K,1,I) - GMAT_22(K,2,I)
       ENDDO
      ENDDO
      
      DO I = 1,2
      GMAT_11(0,0,I) = 0.5*( 2.0d0*GMAT_11(0,1,I) - GMAT_11(0,2,I) &
                 +  2.d0*GMAT_11(1,0,I) - GMAT_11(2,0,I) )
      GMAT_22(0,0,I) = 0.5*( 2.0d0*GMAT_12(0,1,I) - GMAT_12(0,2,I) &
                 +  2.d0*GMAT_12(1,0,I) - GMAT_12(2,0,I) )
      GMAT_22(0,0,I) = 0.5*( 2.0d0*GMAT_22(0,1,I) - GMAT_22(0,2,I) &
                 +  2.d0*GMAT_22(1,0,I) - GMAT_22(2,0,I) )
          
     
       CJOB_11(0,0,I) = 0.50d0*(2.d0*CJOB_11(0,1,I)-CJOB_11(0,2,I) &
                 +  2.d0*CJOB_11(1,0,I) - CJOB_11(2,0,I) )
       CJOB_12(0,0,I) = 0.50d0*(2.d0*CJOB_12(0,1,I)-CJOB_12(0,2,I) &
                 +  2.d0*CJOB_12(1,0,I) - CJOB_12(2,0,I) )
       CJOB_21(0,0,I) = 0.50d0*(2.d0*CJOB_21(0,1,I)-CJOB_21(0,2,I) &
                 +  2.d0*CJOB_21(1,0,I) - CJOB_21(2,0,I) )
       CJOB_22(0,0,I) = 0.50d0*(2.d0*CJOB_22(0,1,I)-CJOB_22(0,2,I) &
                 +  2.d0*CJOB_22(1,0,I) - CJOB_22(2,0,I) )
      ENDDO
      
! CC Additinal matrix required for MG  
      DO I = 1,2        
      DO K=0,NZ+1
        CJOB_11(K,NY+2,I) = CJOB_11(K,NY+1,I)
	CJOB_12(K,NY+2,I) = CJOB_12(K,NY+1,I)
	CJOB_21(K,NY+2,I) = CJOB_21(K,NY+1,I)
	CJOB_22(K,NY+2,I) = CJOB_22(K,NY+1,I)
      
        GMAT_11(K,NY+2,I) = GMAT_11(K,NY+1,I)
	GMAT_12(K,NY+2,I) = GMAT_12(K,NY+1,I)
	GMAT_22(K,NY+2,I) = GMAT_22(K,NY+1,I)
      ENDDO	 
      DO J=0,NY+1
        CJOB_11(NZ+2,J,I) = CJOB_11(NZ+1,J,I)
	CJOB_12(NZ+2,J,I) = CJOB_12(NZ+1,J,I)
	CJOB_21(NZ+2,J,I) = CJOB_21(NZ+1,J,I)
	CJOB_22(NZ+2,J,I) = CJOB_22(NZ+1,J,I)
      
        GMAT_11(NZ+2,J,I) = GMAT_11(NZ+1,J,I)
	GMAT_12(NZ+2,J,I) = GMAT_12(NZ+1,J,I)
	GMAT_22(NZ+2,J,I) = GMAT_22(NZ+1,J,I)
      ENDDO
      GMAT_11(NZ+2,NY+2,I) = GMAT_11(NZ+1,NY+1,I)
      GMAT_12(NZ+2,NY+2,I) = GMAT_12(NZ+1,NY+1,I)
      GMAT_22(NZ+2,NY+2,I) = GMAT_22(NZ+1,NY+1,I)
      
        CJOB_11(NZ+2,NY+2,I) = CJOB_11(NZ+1,NY+1,I)
	CJOB_12(NZ+2,NY+2,I) = CJOB_12(NZ+1,NY+1,I)
	CJOB_21(NZ+2,NY+2,I) = CJOB_21(NZ+1,NY+1,I)
	CJOB_22(NZ+2,NY+2,I) = CJOB_22(NZ+1,NY+1,I)
      ENDDO 	
     

      IF (LES) THEN

      DO J=1,NY+1
       DO K=1,NZ+1
        DO I=1,2
         GMAT_11_z(K,J,I)   =  (2.0*zeta_x(K,J,I)*zeta_x(K,J,I) &
                     +  zeta_y(K,J,I)*zeta_y(K,J,I))*IN_JACO(K,J,I)
         GMAT_12_z(K,J,I)   =  (2.0*zeta_x(K,J,I)*eta_x(K,J,I) &
                     +  zeta_y(K,J,I)*eta_y(K,J,I))*IN_JACO(K,J,I)
         GMAT_22_z(K,J,I)   =  (2.0*eta_x(K,J,I)*eta_x(K,J,I) &
                     +  eta_y(K,J,I)*eta_y(K,J,I))*IN_JACO(K,J,I)

         GMAT_11_y(K,J,I)   =  (zeta_x(K,J,I)*zeta_x(K,J,I) &
                     +  2.0*zeta_y(K,J,I)*zeta_y(K,J,I))*IN_JACO(K,J,I)
         GMAT_12_y(K,J,I)   =  (zeta_x(K,J,I)*eta_x(K,J,I) &
                     +  2.0*zeta_y(K,J,I)*eta_y(K,J,I))*IN_JACO(K,J,I)
         GMAT_22_y(K,J,I)   =  (eta_x(K,J,I)*eta_x(K,J,I) &
                     +  2.0*eta_y(K,J,I)*eta_y(K,J,I))*IN_JACO(K,J,I)
        ENDDO
       ENDDO
      ENDDO

      DO J=0,NY+1
       DO I=1,2
        GMAT_11_z(0,J,I) = 2.d0*GMAT_11_z(1,J,I) - GMAT_11_z(2,J,I)
        GMAT_12_z(0,J,I) = 2.d0*GMAT_12_z(1,J,I) - GMAT_12_z(2,J,I)
        GMAT_22_z(0,J,I) = 2.d0*GMAT_22_z(1,J,I) - GMAT_22_z(2,J,I)

        GMAT_11_y(0,J,I) = 2.d0*GMAT_11_y(1,J,I) - GMAT_11_y(2,J,I)
        GMAT_12_y(0,J,I) = 2.d0*GMAT_12_y(1,J,I) - GMAT_12_y(2,J,I)
        GMAT_22_y(0,J,I) = 2.d0*GMAT_22_y(1,J,I) - GMAT_22_y(2,J,I)
       ENDDO
      ENDDO

      DO K=0,NZ+1
       DO I=1,2
        GMAT_11_z(K,0,I) = 2.0d0*GMAT_11_z(K,1,I) - GMAT_11_z(K,2,I)
        GMAT_12_z(K,0,I) = 2.0d0*GMAT_12_z(K,1,I) - GMAT_12_z(K,2,I)
        GMAT_22_z(K,0,I) = 2.0d0*GMAT_22_z(K,1,I) - GMAT_22_z(K,2,I)

        GMAT_11_y(K,0,I) = 2.0d0*GMAT_11_y(K,1,I) - GMAT_11_y(K,2,I)
        GMAT_12_y(K,0,I) = 2.0d0*GMAT_12_y(K,1,I) - GMAT_12_y(K,2,I)
        GMAT_22_y(K,0,I) = 2.0d0*GMAT_22_y(K,1,I) - GMAT_22_y(K,2,I)
       ENDDO
      ENDDO

      DO I = 1,2
      GMAT_11_z(0,0,I) = 0.5*( 2.0d0*GMAT_11_z(0,1,I) - GMAT_11_z(0,2,I) &
                 +  2.d0*GMAT_11_z(1,0,I) - GMAT_11_z(2,0,I) )
      GMAT_22_z(0,0,I) = 0.5*( 2.0d0*GMAT_12_z(0,1,I) - GMAT_12_z(0,2,I) &
                 +  2.d0*GMAT_12_z(1,0,I) - GMAT_12_z(2,0,I) )
      GMAT_22_z(0,0,I) = 0.5*( 2.0d0*GMAT_22_z(0,1,I) - GMAT_22_z(0,2,I) &
                 +  2.d0*GMAT_22_z(1,0,I) - GMAT_22_z(2,0,I) )

      GMAT_11_y(0,0,I) = 0.5*( 2.0d0*GMAT_11_y(0,1,I) - GMAT_11_y(0,2,I) &
                 +  2.d0*GMAT_11_y(1,0,I) - GMAT_11_y(2,0,I) )
      GMAT_22_y(0,0,I) = 0.5*( 2.0d0*GMAT_12_y(0,1,I) - GMAT_12_y(0,2,I) &
                 +  2.d0*GMAT_12_y(1,0,I) - GMAT_12_y(2,0,I) )
      GMAT_22_y(0,0,I) = 0.5*( 2.0d0*GMAT_22_y(0,1,I) - GMAT_22_y(0,2,I) &
                 +  2.d0*GMAT_22_y(1,0,I) - GMAT_22_y(2,0,I) )
      ENDDO

      DO I = 1,2
      DO K=0,NZ+1
        GMAT_11_z(K,NY+2,I) = GMAT_11_z(K,NY+1,I)
        GMAT_12_z(K,NY+2,I) = GMAT_12_z(K,NY+1,I)
        GMAT_22_z(K,NY+2,I) = GMAT_22_z(K,NY+1,I)

        GMAT_11_y(K,NY+2,I) = GMAT_11_y(K,NY+1,I)
        GMAT_12_y(K,NY+2,I) = GMAT_12_y(K,NY+1,I)
        GMAT_22_y(K,NY+2,I) = GMAT_22_y(K,NY+1,I)
      ENDDO
      DO J=0,NY+1
        GMAT_11_z(NZ+2,J,I) = GMAT_11_z(NZ+1,J,I)
        GMAT_12_z(NZ+2,J,I) = GMAT_12_z(NZ+1,J,I)
        GMAT_22_z(NZ+2,J,I) = GMAT_22_z(NZ+1,J,I)

        GMAT_11_y(NZ+2,J,I) = GMAT_11_y(NZ+1,J,I)
        GMAT_12_y(NZ+2,J,I) = GMAT_12_y(NZ+1,J,I)
        GMAT_22_y(NZ+2,J,I) = GMAT_22_y(NZ+1,J,I)
      ENDDO
      GMAT_11_z(NZ+2,NY+2,I) = GMAT_11_z(NZ+1,NY+1,I)
      GMAT_12_z(NZ+2,NY+2,I) = GMAT_12_z(NZ+1,NY+1,I)
      GMAT_22_z(NZ+2,NY+2,I) = GMAT_22_z(NZ+1,NY+1,I)

      GMAT_11_y(NZ+2,NY+2,I) = GMAT_11_y(NZ+1,NY+1,I)
      GMAT_12_y(NZ+2,NY+2,I) = GMAT_12_y(NZ+1,NY+1,I)
      GMAT_22_y(NZ+2,NY+2,I) = GMAT_22_y(NZ+1,NY+1,I)
      ENDDO

      ENDIF
 
 



      open(10,file='sample.dat',form='formatted',status='unknown')
      write(10,*) 'title = "sample mesh" '
      write(10,*) 'variables = "x", "y"'
      write(10,*) 'zone i=',NZ+2 ,', j=', NY+2,'DATAPACKING=POINT'
      
      do j=0,NY+1
       do i=0,NZ+1
         write(10,111) xpoint(i,j),ypoint(i,j)
       enddo
      enddo
      close(10)
      
! C111   format(2f12.8)      
	
112   format(f16.10)      
      


! C     Calculating sponge in both direction
      SPONGE_SIGMA_OUT = 0.0d0
      sponge_sigma     = 0.0d0
      sponge_temp      = 0.0d0

      SPONGE_EAST = .TRUE.
      SPONGE_WEST = .FALSE.

      NZ_min = 10
      NZ_max = (NZ+2)-10
      co_b   = 3.0d0
      b_min  = log(co_b)/(REAL(xpoint(NZ_min,NY+1)-xpoint(0,NY+1)))
      b_max  = log(co_b)/(REAL(xpoint(NZ+1,NY+1)-xpoint(NZ_max,NY+1)))

      do j =0,NY+1
      DO K = 0, NZ+1
        IF ( K .GE. NZ_max ) THEN
         IF (SPONGE_EAST) THEN
          SPONGE_SIGMA_OUT(K,J) = exp(b_max*(DABS(xpoint(K,NY+1)  &
                -xpoint(NZ_max,NY+1)))) -1.0
         ELSE                                                  !  this is to eliminate sponge on east boundary 
          SPONGE_SIGMA_OUT(K,J) = SPONGE_SIGMA_OUT(K,J)
          sponge_temp(K,J) = sponge_temp(K,J)
         ENDIF
        ELSEIF( K .LE. NZ_min ) THEN                            !comment this out to eliminate sponge on west boundary
         IF (SPONGE_WEST) THEN
          SPONGE_SIGMA_OUT(K,J) = exp(b_min*(DABS(xpoint(K,NY+1)  &
                -xpoint(NZ_min,NY+1)))) -1.0
         ELSE                                                  !  this is to eliminate sponge on west boundary
          SPONGE_SIGMA_OUT(K,J) = SPONGE_SIGMA_OUT(K,J)
          sponge_temp(K,J) = (exp(b_min*(DABS(xpoint(K,NY+1)  &
                -xpoint(NZ_min,NY+1)))) -1.0d0)/(co_b-1.0d0)
         ENDIF
        ELSE
          SPONGE_SIGMA_OUT(K,J) = 0.0d0
          sponge_temp(K,J) = 0.0d0
        ENDIF
      ENDDO
      enddo


      SPONGE_NORTH = .FALSE.
      SPONGE_SOUTH = .FALSE.
      WAVE_ABS = .FALSE.
      NY_min  = 20
      NY_max = (NY+2)-20
      co_b = 2.0d0 
      b_min    = log(co_b)/(REAL(ypoint(NZ+1,NY_min)-ypoint(NZ+1,0)))
      b_max    = log(co_b)/(REAL(ypoint(NZ+1,NY+1)-ypoint(NZ+1,NY_max)))

      do k=0,NZ+1
      DO J = 0, NY+1
        IF ( J .GT. NY_max ) THEN
          IF (SPONGE_NORTH) THEN
           SPONGE_SIGMA_OUT(K,J) = SPONGE_SIGMA_OUT(K,J)  &
                                   + exp(b_max*(DABS(ypoint(NZ+1,J) -ypoint(NZ+1,NY_max)))) -1.0
!                                  + exp(b_max*(REAL(J-NY_max))) -1.0
          ELSE                                                  !  this is to eliminate sponge on north boundary
           SPONGE_SIGMA_OUT(K,J) = SPONGE_SIGMA_OUT(K,J)
           sponge_sigma(K,J) = sponge_sigma(K,J) &
                               + exp(b_max*(DABS(ypoint(NZ+1,J) -ypoint(NZ+1,NY_max)))) -1.0 
          ENDIF
         ELSEIF( J .LE. NY_min ) THEN                            !comment this out to eliminate sponge on south boundary
          IF (SPONGE_SOUTH) THEN
           SPONGE_SIGMA_OUT(K,J) = SPONGE_SIGMA_OUT(K,J)  &
                                   +  exp(b_min*(DABS(ypoint(NZ+1,J)  &
                                      -ypoint(NZ+1,NY_min)))) -1.0

!                                   + exp(b_min*(REAL(NY_min-J))) -1.0
          ELSE                                                   !  this is to eliminate sponge on south boundary
           SPONGE_SIGMA_OUT(K,J) = SPONGE_SIGMA_OUT(K,J)  
           sponge_sigma(K,J) = sponge_sigma(K,J) &
                                      +  exp(b_min*(DABS(ypoint(NZ+1,J) -ypoint(NZ+1,NY_min)))) -1.0
          ENDIF
         ELSE
         ENDIF
      ENDDO
      enddo

      SPONGE_SIGMA = 0.0d0 ;
      SPONGE_temp = 0.0d0 ;

      if(rank .eq. 0) then
       open(654,file='sponge.dat',form='formatted',status='unknown')
       write(654,*) 'title = "sample mesh" '
       write(654,*) 'variables = "x","y", "sponge", "sponge_sigma","sponge_temp","rho_backg" '
       write(654,*) 'zone i=',NZ+2 ,', j=', NY+2,'DATAPACKING=POINT'

       do j=0,NY+1
        do i=0,NZ+1
          write(654,665) xpoint(i,j),ypoint(i,j), &
          SPONGE_SIGMA_OUT(i,j),sponge_sigma(i,j),sponge_temp(i,j),TH_BAR(i,j,1)
        enddo
       enddo
       close(654)
665    format(6f14.8)
       call  plane_parav_sponge
     endif


      IF (V_BC_ZMIN .EQ. 9) THEN
         call wave_beam_boundary
      ELSEIF (W_BC_ZMIN .EQ. 9) THEN
         call wave_beam_boundary
      ENDIF


      deallocate (x_zeta, y_zeta, &
         x_eta, y_eta, zeta_x, zeta_y, eta_x, eta_y, IN_JACO)
!      deallocate (xpoint, ypoint)
      write(6,*) 'Jacobian tranformation is done'

      
      RETURN
      END 

!------------------------------------
      SUBROUTINE wave_beam_boundary
!-----------------------------------
use Domain
use forcing
use mpi_var, only : RANK
      implicit none
!      real*8 xpoint(0:NZ+1,0:NY+1), ypoint(0:NZ+1,0:NY+1)
      real*8 z !,y(0:NY+1)
!      real*8 zeta(0:NY+1),eta(0:NY+1)
!      real*8 A(0:NY+1), B(0:NY+1), C(0:NY+1)
      real*8 beam_loc, theta, omega, grav, nu, eta_scale, N_infty, alpha, rho_0
      real*8 k_bot, k_top, amp_l !,K_X(0:NY+1)
      real*8 l, u_0,norm_factor,pi,H_W,temp_force
      integer kk,j,i,kkk   
      logical PARTIAL_FORCING
      Parameter (Kk = 250)

      real(r8),allocatable,dimension(:,:)   :: xpoint, ypoint
      real(r8),allocatable,dimension(:)   :: K_X,A,B,C,zeta,eta,y
      complex*16,allocatable, dimension(:) ::  u_ac,v_ac,phi_prime

      allocate (xpoint(0:NZ+1,0:NY+1), ypoint(0:NZ+1,0:NY+1))   
      allocate (zeta(0:NY+1),eta(0:NY+1),y(0:NY+1),A(0:NY+1), &
               B(0:NY+1), C(0:NY+1))     
      allocate ( K_X(0:Kk) )
      allocate (u_ac(0:NY+1),v_ac(0:NY+1),phi_prime(0:NY+1))

!     y,v
!      ^
!      |  -. zeta,v'
!  th-->-/
!      |/ 
!      .---->z,u
!       \
!        \
!         -' eta,u'
!default wavelength is 0.2

      PARTIAL_FORCING = .FALSE.
      write(6,*) 'Calculating beam'

      open(21,file='GRID_XZ.dat',form='formatted',status='old')
      xpoint(:,:) =0.0d0
      ypoint(:,:) =0.0d0
      DO J=0,NY+1
       DO Kkk=0,NZ+1
         READ(21,*)xpoint(kkK,J),ypoint(kkK,J)
       ENDDO
      ENDDO
      close(21) 

      u_0 = 1.0d0*10.0d0**(-03)
      pi = 4.d0*ATAN(1.d0)
      beam_loc = 0.23  ! distance from floor to generation zone (% of total)
      theta = 30.0*pi/180.0
      omega = 2.0d0*pi/(12.4*60.0d0**2.0)
      
      nu = 0.04*10.0**(-9.0)

!scale to get desired wave number
      eta_scale = 1.0

!Kk is the resolution for E_k and k_bot/top are the limits of integration for E_k
!      Kk = 250

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     IF (RANK .eq. 0) THEN
      write(6,*) '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
!      write(6,*) '%%%%%%% Beam Type = thomas_steven',thomas_steven
      write(6,*) '%%%%%%%        Beam parameter   %%%%%%%%%%%%%%%%%%%'
      write(6,*) 'U_0 = ', u_0, 'Location = ', beam_loc, 'Angle = ', theta*180.0d0/pi 
      write(6,*) 'Omega = ', omega 
      write(6,*) '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
     ENDIF

      N_infty = omega/dsin(theta)
      grav=10.0d0 ! N_infty*N_infty
      rho_0 = 1000.0d0 

      k_bot = 0.0 
      k_top = 2.0

      l = 160.0
      im = (0.0,1.0)

!coordinate transformations
! beam at 1st Quadrant  + l*dsin(theta) 
! beam at 4th Quadrant  - l*dsin(theta)

      do j=0,NY+1
!         y(j) = (ypoint(0,j)/abs(ypoint(0,0))) &
!           - 1.0 + l*dsin(theta) + (1.0-beam_loc)*2.0
          y(j) = ((ypoint(0,j)-ypoint(NY+1,j)+ypoint(0,j))  &
           /abs(ypoint(0,0)-ypoint(NY+1,j)+ypoint(0,j)))   &
           - 1.0   + l*dsin(theta) + (1.0-beam_loc)*2.0
      enddo
     
      z = l*dcos(theta)
      alpha = (2.0*N_infty*dcos(theta)/nu)**(1.0/3.0)

      do i = 0,Kk
        K_X(i) = real(i)*k_top/real(Kk) + k_bot
      enddo

      write(6,*) 'Inputs assigned'

!to get beam at 1st quadrant put: zeta(i,k) = z(i,k)*cos(theta) + y(i,k)*sin(theta) + l;
!                                % eta(i,k) = z(i,k)*sin(theta) - y(i,k)*cos(theta);
!to get beam at 4th quadrant put: eta(i,k)  = z(i,k)*sin(theta) - y(i,k)*cos(theta) + l;
!                                % eta(i,k) = z(i,k)*sin(theta) + y(i,k)*cos(theta);

      do j=0,NY+1
        zeta(j) = z*dcos(theta) + y(j)*dsin(theta) + l
        eta(j) = (z*dsin(theta) - y(j)*dcos(theta))/eta_scale
        u_prime(j)=0.0
        v_prime(j)=0.0
        beam_th(j)=0.0
        u_ac(j)=0.0
        v_ac(j)=0.0
      enddo
  
!      min_l_q1 = max(max(z))*dsin(theta)
!      max_l_q1 = -max(max(x))*dcos(theta)

      write(6,*) 'Grid set'
      write(6,*) 'Calculating Variables'

      do j = 0,NY+1
        A(j) = (zeta(j)*(N_infty**2.0)*dsin(theta)/grav)**(-2.0/3.0)
        B(j) = (zeta(j)*(N_infty**2.0)*dsin(theta)*zeta(j)*alpha**(3.0/2.0)/grav)**(-2.0/3.0)
        C(j) = -(N_infty*rho_0/grav)*(zeta(j)*(N_infty**2.0)* & 
                 dsin(theta)/grav)**(-2.0/3.0)
      enddo

      IF (PARTIAL_FORCING) THEN
      B = 0.d0
      C = 0.d0
      ENDIF

!      C = 0.d0
 
      do j = 0,NY+1
        do i = 0,Kk
          u_prime(j) = u_prime(j) + A(j)*(1.0    &
            *K_X(i)*exp(-K_X(i)**3.0)*exp(im*K_X(i)*alpha*eta(j)/ &
             (zeta(j)**(1.0/3.0))))

          v_prime(j) = v_prime(j) + B(j)*(-im    &
            *K_X(i)**3.0*exp(-K_X(i)**3.0)*exp(im*K_X(i)*alpha*eta(j)/ &
             (zeta(j)**(1.0/3.0))))

          beam_th(j) = beam_th(j) + C(j)*(-im    &
            *K_X(i)*exp(-K_X(i)**3.0)*exp(im*K_X(i)*alpha*eta(j)/ & 
             (zeta(j)**(1.0/3.0))))
        enddo
!%      if real(u_prime(i,k)) > 1000.0
!%         imp_part = imag(u_prime(i,k)) ;
!%         u_prime(i,k)=1000.0 + j*imp_part ;
!%      end
!%      if real(v_prime(i,k)) > 1000.0
!%         imp_part = imag(v_prime(i,k)) ;
!%         v_prime(i,k)=1000.0 + j*imp_part ;
!%      end 
!%       if real(u_prime(i,k)) < -1000.0
!%         imp_part = imag(u_prime(i,k)) ;
!%         u_prime(i,k)=-1000.0 + j*imp_part ;
!%      end
!%      if real(v_prime(i,k)) < -1000.0
!%         imp_part = imag(v_prime(i,k)) ;
!%         v_prime(i,k)=-1000.0 + j*imp_part ;
!%      end
!%
      enddo

      IF ( PARTIAL_FORCING ) THEN
      ELSE      
       do j = 0,NY+1
        u_ac(j)    = u_prime(j)*dcos(theta) !+ v_prime(j)*dsin(theta)
        v_ac(j)    = u_prime(j)*dsin(theta) !- v_prime(j)*dcos(theta)
       enddo
       do j = 0,NY+1
        u_prime(j) =  u_ac(j)
        v_prime(j) =  v_ac(j)
       enddo

      ENDIF  

      norm_factor = 0
      do j = 0,NY+1
        if (abs(real(u_prime(j))) .GT. norm_factor) THEN
           norm_factor = abs(real(u_prime(j)))
        endif
      enddo

!scaling
      do j = 0,NY+1
        H_W       = (1.d0-(2.d0*(ypoint(0,j)-ypoint(0,0))/ &
                    (ypoint(0,NY+1)-ypoint(0,0))-1.0)**30.0)**6.0
        u_prime(j)=u_0*u_prime(j)*H_W/norm_factor
        v_prime(j)=u_0*v_prime(j)*H_W/norm_factor
        beam_th(j) = u_0*beam_th(j)*H_W/norm_factor
      IF (RANK == 0) THEN
!        write(6,*) 'H_W', H_w,grav
      ENDIF
      enddo

      IF (RANK == 0) THEN
       open(66,file='beam_profile.dat',form='formatted',status='unknown')
       do j = 0,NY+1
        write(66,666) ypoint(0,j),real(u_prime(j)), real(v_prime(j)),real(beam_th(j))
       enddo
       close(66)
      ENDIF
 
666   format(4f14.8) 
    
      deallocate (xpoint, ypoint)
      deallocate (A,B,C,Y,zeta,eta)
      deallocate (K_X)
      deallocate (u_ac,v_ac,phi_prime)

      write(6,*) 'Beam Complete'
      
                   
      RETURN
      END




!------------------------------------
      SUBROUTINE wave_beam_boundary_old
!-----------------------------------
use Domain
use forcing
use mpi_var, only : RANK
      implicit none
!      real*8 xpoint(0:NZ+1,0:NY+1), ypoint(0:NZ+1,0:NY+1)
      real*8 z !,y(0:NY+1)
!      real*8 zeta(0:NY+1),eta(0:NY+1)
!      real*8 A(0:NY+1), B(0:NY+1), C(0:NY+1)
      real*8 beam_loc, theta, omega, grav, nu, eta_scale, N_infty, alpha
      real*8 k_bot, k_top !,K_X(0:NY+1)
      real*8 l, u_0,norm_factor,pi,H_W,temp_force
      integer kk,j,i,kkk   
      logical PARTIAL_FORCING
      Parameter (Kk = 250)

      real(r8),allocatable,dimension(:,:)   :: xpoint, ypoint
      real(r8),allocatable,dimension(:)   :: K_X,A,B,C,zeta,eta,y

      allocate (xpoint(0:NZ+1,0:NY+1), ypoint(0:NZ+1,0:NY+1))   
      allocate (zeta(0:NY+1),eta(0:NY+1),y(0:NY+1),A(0:NY+1), &
               B(0:NY+1), C(0:NY+1))     
      allocate ( K_X(0:Kk) )
!     y,v
!      ^
!      |  -. zeta,v'
!  th-->-/
!      |/ 
!      .---->z,u
!       \
!        \
!         -' eta,u'
!default wavelength is 0.2

      PARTIAL_FORCING = .TRUE.
      write(6,*) 'Calculating beam'

      open(21,file='GRID_XZ.dat',form='formatted',status='old')
      xpoint(:,:) =0.0d0
      ypoint(:,:) =0.0d0
      DO J=0,NY+1
       DO Kkk=0,NZ+1
         READ(21,*)xpoint(kkK,J),ypoint(kkK,J)
       ENDDO
      ENDDO
      close(21) 


      pi = 3.1415
      beam_loc = 0.3  ! distance from floor to generation zone (% of total)
      theta = 30.0*pi/180.0
      omega = 4.0d0*pi
      
      nu = 0.6*10.0**(-5.0)

!scale to get desired wave number
      eta_scale = 1.0

!Kk is the resolution for E_k and k_bot/top are the limits of integration for E_k
!      Kk = 250

 

      k_bot = 0.0 
      k_top = 2.0

      l = 100.0
      u_0 = 0.6
      im = (0.0,1.0)

!coordinate transformations
! beam at 1st Quadrant  + l*dsin(theta) 
! beam at 4th Quadrant  - l*dsin(theta)

      do j=0,NY+1
         y(j) = (ypoint(0,j)/abs(ypoint(0,0))) &
           - 1.0 + l*dsin(theta) + (1.0-beam_loc)*2.0
      enddo
     
      z = l*dcos(theta)

      N_infty = omega/dsin(theta)
      grav=N_infty*N_infty

      alpha = (2.0*N_infty*dcos(theta)/nu)**(1.0/3.0)

      do i = 0,Kk
        K_X(i) = real(i)*k_top/real(Kk) + k_bot
      enddo

      write(6,*) 'Inputs assigned'

!to get beam at 1st quadrant put: zeta(i,k) = z(i,k)*cos(theta) + y(i,k)*sin(theta) + l;
!                                % eta(i,k) = z(i,k)*sin(theta) - y(i,k)*cos(theta);
!to get beam at 4th quadrant put: eta(i,k)  = z(i,k)*sin(theta) - y(i,k)*cos(theta) + l;
!                                % eta(i,k) = z(i,k)*sin(theta) + y(i,k)*cos(theta);

      do j=0,NY+1
        zeta(j) = z*dcos(theta) + y(j)*dsin(theta) + l
        eta(j) = (z*dsin(theta) - y(j)*dcos(theta))/eta_scale
        u_prime(j)=0.0
        v_prime(j)=0.0
        beam_th(j)=0.0
      enddo
  
!      min_l_q1 = max(max(z))*dsin(theta)
!      max_l_q1 = -max(max(x))*dcos(theta)

      write(6,*) 'Grid set'
      write(6,*) 'Calculating Variables'

      do j = 0,NY+1
        A(j) = (zeta(j)*(N_infty**2.0)*dsin(theta)/grav)**(-2.0/3.0)
        B(j) = (zeta(j)*(N_infty**2.0)*dsin(theta)*zeta(j)*alpha**(3.0/2.0)/grav)**(-2.0/3.0)
        C(j) = -(N_infty/grav)*(zeta(j)*(N_infty**2.0)* & 
                 dsin(theta)/grav)**(-2.0/3.0)
      enddo

      IF (PARTIAL_FORCING) THEN
      B = 0.d0
      C = 0.d0
      ENDIF

!      C = 0.d0
 
      do j = 0,NY+1
        do i = 0,Kk
          u_prime(j) = u_prime(j) + A(j)*(1.0    &
            *K_X(i)*exp(-K_X(i)**3.0)*exp(im*K_X(i)*alpha*eta(j)/ &
             (zeta(j)**(1.0/3.0))))

          v_prime(j) = v_prime(j) + B(j)*(-im    &
            *K_X(i)**3.0*exp(-K_X(i)**3.0)*exp(im*K_X(i)*alpha*eta(j)/ &
             (zeta(j)**(1.0/3.0))))

          beam_th(j) = beam_th(j) + C(j)*(-im    &
            *K_X(i)*exp(-K_X(i)**3.0)*exp(im*K_X(i)*alpha*eta(j)/ & 
             (zeta(j)**(1.0/3.0))))
        enddo
!%      if real(u_prime(i,k)) > 1000.0
!%         imp_part = imag(u_prime(i,k)) ;
!%         u_prime(i,k)=1000.0 + j*imp_part ;
!%      end
!%      if real(v_prime(i,k)) > 1000.0
!%         imp_part = imag(v_prime(i,k)) ;
!%         v_prime(i,k)=1000.0 + j*imp_part ;
!%      end 
!%       if real(u_prime(i,k)) < -1000.0
!%         imp_part = imag(u_prime(i,k)) ;
!%         u_prime(i,k)=-1000.0 + j*imp_part ;
!%      end
!%      if real(v_prime(i,k)) < -1000.0
!%         imp_part = imag(v_prime(i,k)) ;
!%         v_prime(i,k)=-1000.0 + j*imp_part ;
!%      end
!%
      enddo

      IF ( PARTIAL_FORCING ) THEN
      ELSE      
       do j = 0,NY+1
        temp_force = u_prime(j)*dcos(theta) + v_prime(j)*dsin(theta)
        v_prime(j) = u_prime(j)*dsin(theta) - v_prime(j)*dcos(theta)
        u_prime(j) = temp_force
       enddo
      ENDIF  

      norm_factor = 0
      do j = 0,NY+1
        if (abs(real(u_prime(j))) .GT. norm_factor) THEN
           norm_factor = abs(real(u_prime(j)))
        endif
      enddo

!scaling
      do j = 0,NY+1
        H_W       = (1.d0-(2.d0*(ypoint(0,j)-ypoint(0,0))/ &
                    (ypoint(0,NY+1)-ypoint(0,0))-1.0)**30.0)**6.0
        u_prime(j)=u_0*u_prime(j)*H_W/norm_factor
        v_prime(j)=u_0*v_prime(j)*H_W/norm_factor
        beam_th(j) = u_0*beam_th(j)*H_W/norm_factor
      IF (RANK == 0) THEN
!        write(6,*) 'H_W', H_w,grav
!       write(6,*) u_prime(j), beam_th(j)
      ENDIF
      enddo
     
      deallocate (xpoint, ypoint)
      deallocate (A,B,C,Y,zeta,eta)
      deallocate (K_X)

      write(6,*) 'Beam Complete'
     
                   
      RETURN
      END



!--------------------------------------------------------------
            SUBROUTINE open_MP_initialization
!--------------------------------------------------------------
use ntypes
use Domain
use Grid
use Fft_var
use TIME_STEP_VAR
use run_variable
use omp_lib
use mpi_var, only : rank
implicit none

!--.---------.---------.---------.---------.---------.-|-------|
!$OMP PARALLEL PRIVATE(iam, nthreads, chunk)
! Compute the subset of iterations
! executed by each thread

!      nthreads = omp_get_num_threads()
!      iam = omp_get_thread_num()
!      chunk = (NKX+1)/nthreads

!      kstart = iam * chunk
!      kend = min((iam + 1) * chunk-1, NKX)

!      if(iam == nthreads-1)then
!      kend = NKX
!      endif

!      chunk = (NXM+1)/nthreads

!      k_start = iam * chunk
!      k_end = min((iam + 1) * chunk-1, NXM)

!      if(iam == nthreads-1)then
!      k_end = NXM
!      endif

!      chunk = (JEND-JSTART+1)/nthreads

!      j_start = iam * chunk +1 
!      j_end = min((iam + 1) * chunk, JEND)

!      if(iam == nthreads-1)then
!      j_end = JEND
!      endif

!      chunk = (NY+2)/nthreads
         
!      jj_start = iam * chunk  
!      jj_end = min((iam + 1) * chunk-1, NY+1)
            
!      if(iam == nthreads-1)then
!      jj_end = NY+1
!      endif 

!       write(6,*)'LOOPs', iam,'kstart',kstart,'kend',kend
!      write(6,*)'LOOPs',iam,'k_start',k_start,'k_end',k_end
!      write(6,*)'LOOPs',iam,'j_start',j_start,'j_end',j_end
      if (rank .eq. 0) then 
       write(6,*)'LOOPs',iam,'jj_start',jj_start,'jj_end',jj_end
      endif
!$omp end parallel

      write(6,*)'Threads are intialized'

      RETURN
      END
 
