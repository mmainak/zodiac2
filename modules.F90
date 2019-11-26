subroutine global_allocation  

return
end

!@c
!DATA TYPES
module ntypes
 integer, parameter :: r4=4
 integer, parameter :: r8=8 
 integer, parameter :: i4=4 
end module ntypes 

!DOMAIN
module Domain
 use ntypes

 integer(i4)       :: nx, ny, nz, N_TH, NXM, NYM, NZM, TNKZ, TNKY 
 integer(i4)       :: NKX
 integer(i4)       :: NP,NXP,NZP,NXV,NZV, NKXV, NKXP,NX2V,NX2P,NXP_L,NX2P_L

 INCLUDE   'grid_def'
 INCLUDE 'mpif.h'

end module Domain

!GRID
module Grid
 use ntypes


 real(r8)          :: LX, LY, LZ, CSX, CSY, CSZ         !Length
 INTEGER           :: jstart, jend , ZSTART,  ZEND                       ! 

 real(r8),allocatable,dimension(:) :: GX,  GY,  GZ,  DX,  DY,  DZ, &
 GXF, GYF, GZF, DXF, DYF, DZF

 real(r8),allocatable,dimension(:,:)   :: xpoint, ypoint
end module Grid

!TIME_STEP_VAR
module TIME_STEP_VAR
 use ntypes
real(r8) DELTA_T, DELTA_T_in
integer  N_TIME_STEPS, NUM_PER_DIR, TIME_AD_METH

end module TIME_STEP_VAR

!Fft
module fft_var
use ntypes, only : r8

!FFT plans
 integer(8) :: FFTW_X_TO_P_PLAN, FFTW_X_TO_F_PLAN, &
               FFTW_Y_TO_P_PLAN, FFTW_Y_TO_F_PLAN, &
               FFTW_Z_TO_P_PLAN, FFTW_Z_TO_F_PLAN
 
 integer    NKY, NKZ
 real(r8)   PI, EPS, RNX, RNY, RNZ
! PARAMETER  (NKX=NX/3) 

!FFT VARIABLES
 !REAL FIELDS
 real(r8),allocatable,dimension(:,:,:)    :: FIELD_IN
 real(r8),allocatable,dimension(:,:)      :: PLANE_IN

 !COMPLEX FIELDS
 complex(r8)                              :: CI  
 complex(r8),allocatable,dimension(:,:)   :: CZX_PLANE, CYZ_PLANE
 complex(r8),allocatable,dimension(:)     :: CIKX, CIKY, CIKZ,CIKXP

!WAVE INDICES  
 real(r8),allocatable,dimension(:)    :: kx, ky, kz,kxp

!WAVE NUMBERS 
 real(r8),allocatable,dimension(:)    :: rkx, rky, rkz

!MODIFIED WAVE NUMBERS
 real(r8),allocatable,dimension(:)    ::  kx2, ky2, kz2,kx2p

      common /fft_pln/ FFTW_X_TO_P_PLAN, FFTW_X_TO_F_PLAN
!$omp threadprivate(/fft_pln/)
end module fft_var

!Run_vari
module run_variable
use ntypes, only : r8
use Domain
!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
! Input parameters and runtime variables
!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      LOGICAL  ::USE_MPI,CONT_STRAT, SPONGE_NORTH, SPONGE_SOUTH, SPONGE_EAST,SPONGE_WEST, MELTING_MODEL
      real(r8) :: NU, KICK, UBULK0, PX0,U_BC_XMIN_C1, U_BC_XMIN_C2, U_BC_XMIN_C3,       &
                V_BC_XMIN_C1,V_BC_XMIN_C2,V_BC_XMIN_C3, W_BC_XMIN_C1, W_BC_XMIN_C2,   &
                W_BC_XMIN_C3,U_BC_YMIN_C1,U_BC_YMIN_C2,U_BC_YMIN_C3, V_BC_YMIN_C1,    &
                V_BC_YMIN_C2, V_BC_YMIN_C3, W_BC_YMIN_C1, W_BC_YMIN_C2, W_BC_YMIN_C3, &
                U_BC_ZMIN_C1, U_BC_ZMIN_C2, U_BC_ZMIN_C3, V_BC_ZMIN_C1, V_BC_ZMIN_C2, &
                V_BC_ZMIN_C3, W_BC_ZMIN_C1, W_BC_ZMIN_C2, W_BC_ZMIN_C3, TH_BC_XMIN_C1(1:N_TH), &
                TH_BC_XMIN_C2(1:N_TH), TH_BC_XMIN_C3(1:N_TH), TH_BC_YMIN_C1(1:N_TH),  &
                TH_BC_YMIN_C2(1:N_TH), TH_BC_YMIN_C3(1:N_TH), TH_BC_ZMIN_C1(1:N_TH),  &
                TH_BC_ZMIN_C2(1:N_TH), TH_BC_ZMIN_C3(1:N_TH)

      real(r8)  :: U_BC_XMAX_C1, U_BC_XMAX_C2, U_BC_XMAX_C3, &
                V_BC_XMAX_C1, V_BC_XMAX_C2, V_BC_XMAX_C3, &
                W_BC_XMAX_C1, W_BC_XMAX_C2, W_BC_XMAX_C3, &
                U_BC_YMAX_C1, U_BC_YMAX_C2, U_BC_YMAX_C3, &
                V_BC_YMAX_C1, V_BC_YMAX_C2, V_BC_YMAX_C3, &
                W_BC_YMAX_C1, W_BC_YMAX_C2, W_BC_YMAX_C3, &
                U_BC_ZMAX_C1, U_BC_ZMAX_C2, U_BC_ZMAX_C3, &
                V_BC_ZMAX_C1, V_BC_ZMAX_C2, V_BC_ZMAX_C3, &
                W_BC_ZMAX_C1, W_BC_ZMAX_C2, W_BC_ZMAX_C3
 
      real(r8)  :: TH_BC_XMAX_C1(1:N_TH),TH_BC_XMAX_C2(1:N_TH), &
                TH_BC_XMAX_C3(1:N_TH),TH_BC_YMAX_C1(1:N_TH), &
                TH_BC_YMAX_C2(1:N_TH),TH_BC_YMAX_C3(1:N_TH), &
                TH_BC_ZMAX_C1(1:N_TH), TH_BC_ZMAX_C2(1:N_TH), &
                TH_BC_ZMAX_C3(1:N_TH), CFL,DFN,INT_PI,RI_FINAL(1:N_TH), &
                dtc,C_avg


      real(r8),allocatable,dimension(:,:)    ::   SPONGE_SIGMA,SPONGE_SIGMA_OUT, &
                                                  SPONGE_TEMP  
      real(r8),allocatable,dimension(:)      ::   C_int, C_int_le
      real(r8)                               ::   dthdz_mean(0:NY+1,1:N_TH), m_FP, C_sp_heat, L_heat          

      integer :: NY_S, NX_T,NY_T,NZ_T

      integer :: VERBOSITY, &  
              SAVE_FLOW_INT, SAVE_STATS_INT, IC_TYPE, F_TYPE,COUNT_DATA, &
              U_BC_XMIN, V_BC_XMIN, W_BC_XMIN, TH_BC_XMIN(1:N_TH), &
              U_BC_XMAX, V_BC_XMAX, W_BC_XMAX, TH_BC_XMAX(1:N_TH), &
              U_BC_YMIN, V_BC_YMIN, W_BC_YMIN, TH_BC_YMIN(1:N_TH), &
              U_BC_YMAX, V_BC_YMAX, W_BC_YMAX, TH_BC_YMAX(1:N_TH), &
              U_BC_ZMIN, V_BC_ZMIN, W_BC_ZMIN, TH_BC_ZMIN(1:N_TH), &
              U_BC_ZMAX, V_BC_ZMAX, W_BC_ZMAX, TH_BC_ZMAX(1:N_TH)

      integer :: PREVIOUS_TIME_STEP,SIGNAL_BC,UPDATE_DT

      logical :: VARIABLE_DT,FIRST_TIME, INT_TREAT, WAVE_ABS, STOCHASTIC_FORCING
      logical :: MOVIE,CREATE_NEW_FLOW, Non_linear_ST

      real(r8) :: TIME, START_TIME,END_TIME
      integer  :: TIME_STEP, RK_STEP
      integer  :: TIME_ARRAY(8)
      
      real(r8) :: RI_TAU(1:N_TH), PR(1:N_TH), FACT_AMP,Ratio_gr,Gravity_g,rho_0, &
                  alpha_w,gamma_w,Ratio_gr_a,Ratio_gr_g, Sal_0, theta_0, dSaldz 
      logical  :: CREATE_NEW_TH(1:N_TH), BACKGROUND_GRAD(1:N_TH)
      integer  :: NUM_READ_TH
      integer  :: READ_TH_INDEX(1:N_TH)
      logical  :: FILTER_TH(1:N_TH)
      integer  :: FILTER_INT(1:N_TH)
      integer  :: JSTART_TH(1:N_TH),JEND_TH(1:N_TH),ZSTART_TH(1:N_TH), &
                  ZEND_TH(1:N_TH) 
      real(r8) :: OMEGA0, AMP_OMEGA0, ANG_BETA,In_H0,Q_H0,H0,f_0

      real(r8),allocatable,dimension(:,:) :: U_BC_LOWER_NWM, &
               W_BC_LOWER_NWM,U_BC_UPPER_NWM, &
               W_BC_UPPER_NWM
      real(r8),allocatable,dimension(:) ::  U1_bar
      real(r8),allocatable,dimension(:,:) :: U2_bar, &
              U3_bar
      real(r8),allocatable,dimension(:,:,:) :: TH_BAR 
      complex(r8),allocatable,dimension(:,:) :: temp_mean 

      real(r8) :: UTAU_MEAN_LOWER,UTAU_MEAN_UPPER,TAUWALL_MEAN
      real(r8) :: UTAU_AVE
      
      complex(r8),allocatable,dimension(:,:) :: CU_BC_LOWER_NWM, &
                CU_BC_UPPER_NWM, CW_BC_LOWER_NWM, &
                CW_BC_UPPER_NWM

      LOGICAL :: FILTER_VEL
      INTEGER :: FILTER_VEL_INT, filter_type

!----*|--.---------.---------.---------.---------.---------.---------.-|-------|! Parameters for Large Eddy Simulation
!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      LOGICAL ::LES
      INTEGER :: LES_MODEL_TYPE,LES_MODEL_TYPE_TH
      real(r8),allocatable,dimension(:,:,:) :: NU_T
      real(r8),allocatable,dimension(:,:,:,:) :: KAPPA_T
      integer ::  J1,J2

!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
! RKW3 parameters
!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      real(r8)  H_BAR(3), BETA_BAR(3), ZETA_BAR(3)

!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
! Global variables
!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
!     ARRAYS FOR WHEN DATA IS SPLITED IN X DIRECTION

!      real(r8),target ::  U1 (0:NX+1,0:NZ+1,0:NY+1)
!      real(r8),target ::  U2 (0:NX+1,0:NZ+1,0:NY+1)
!      real(r8),target ::  U3 (0:NX+1,0:NZ+1,0:NY+1)
!      real(r8),target ::  P  (0:NX+1,0:NZ+1,0:NY+1)
!      real(r8),target ::  TH (0:NX+1,0:NZ+1,0:NY+1,1:N_TH)

!      real(r8)        ::  U1b(0:NX+1,0:NZ+1,0:NY+1)
!      real(r8)        ::  U2b(0:NX+1,0:NZ+1,0:NY+1)
!      real(r8)        ::  U3b(0:NX+1,0:NZ+1,0:NY+1)

!     ARRAYS FOR WHEN DATA IS SPLITED IN X DIRECTION

!      complex(r8),target :: CU1(0:NX/2,0:NZ+1,0:NY+1)
!      complex(r8),target :: CU2(0:NX/2,0:NZ+1,0:NY+1)
!      complex(r8),target :: CU3(0:NX/2,0:NZ+1,0:NY+1)
!      complex(r8),target :: CP (0:NX/2,0:NZ+1,0:NY+1)
!      complex(r8),target :: CTH(0:NX/2,0:NZ+1,0:NY+1,1:N_TH)
  
!      complex(r8)        ::  CU1b(0:NX/2,0:NZ+1,0:NY+1)
!      complex(r8)        ::  CU2b(0:NX/2,0:NZ+1,0:NY+1) 
!      complex(r8)        ::  CU3b(0:NX/2,0:NZ+1,0:NY+1) 

!      EQUIVALENCE (U1,CU1), (U2,CU2), (U3,CU3), &
!        (U1b,CU1b), (U2b,CU2b), (U3b,CU3b),     &
!        (TH, CTH), (P,CP) 


!     real(r8) R1 (0:NX+1,0:NZ+1,0:NY+1), R2 (0:NX+1,0:NZ+1,0:NY+1),    &
!             R3 (0:NX+1,0:NZ+1,0:NY+1), F1 (0:NX+1,0:NZ+1,0:NY+1),     &
!             F2 (0:NX+1,0:NZ+1,0:NY+1), F3 (0:NX+1,0:NZ+1,0:NY+1),     &
!             S1 (0:NX+1,0:NZ+1,0:NY+1), S2 (0:NX+1,0:NZ+1,0:NY+1)

!     real(r8)FTH (0:NX+1,0:NZ+1,0:NY+1,1:N_TH),                      &
!             RTH (0:NX+1,0:NZ+1,0:NY+1,1:N_TH)



!     complex(r8) CR1(0:NX/2,0:NZ+1,0:NY+1), CR2(0:NX/2,0:NZ+1,0:NY+1),  &
!                 CR3(0:NX/2,0:NZ+1,0:NY+1), CF1(0:NX/2,0:NZ+1,0:NY+1),  &
!                 CF2(0:NX/2,0:NZ+1,0:NY+1), CF3(0:NX/2,0:NZ+1,0:NY+1),  &
!                 CS1(0:NX/2,0:NZ+1,0:NY+1), CS2(0:NX/2,0:NZ+1,0:NY+1)

!     complex(r8) CFTH(0:NX/2,0:NZ+1,0:NY+1,1:N_TH), &
!                CRTH(0:NX/2,0:NZ+1,0:NY+1,1:N_TH) 


!     EQUIVALENCE (R1,CR1), (R2,CR2) &
!        , (R3,CR3) , (F1,CF1), (F2,CF2), (F3,CF3), (S1,CS1) &
!        , (S2,CS2),  (RTH,CRTH) , (FTH, CFTH)


!      real(r8),dimension(:,:,:),pointer      :: u_p
!      complex(r8),dimension(:,:,:),pointer    :: cu_p

!       EQUIVALENCE (temp,ctemp)

      real(r8),allocatable,dimension(:,:) :: TH_BACK

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!       TRAIL VARIABLE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      real(r8),DIMENSION(0:NXP,0:NZV-1,0:NY+1)    ::  S1X,S2X
      complex(r8),DIMENSION(0:NX2P,0:NZV-1,0:NY+1)::  CS1X, CS2X
      
      real(r8),DIMENSION(0:NXV-1,0:NZP,0:NY+1)    ::  S1Z,S2Z 
      complex(r8),DIMENSION(0:NX2V-1,0:NZP,0:NY+1)::  CS1Z,CS2Z
      
      real(r8)    ::  VARP (0:NX+1,0:NZP,0:NY+1)
      complex(r8) ::  CVARP(0:NX/2,0:NZP,0:NY+1)

      EQUIVALENCE                                                                                &
                  (S1X,CS1X,S1Z,CS1Z), (S2X,CS2X,S2Z,CS2Z)


      EQUIVALENCE (VARP,CVARP)

     real(r8),allocatable,dimension(:,:,:), target ::  U1X,U2X,U3X,PX,U2bX,U3bX,R1X,R2X,R3X,F1X,F2X,F3X
     complex(r8),allocatable,dimension(:,:,:), target :: CU1X,CU2X,CU3X,CPX,CU2bX,CU3bX,CR1X,CR2X,CR3X,CF1X,CF2X,CF3X

     real(r8),allocatable,dimension(:,:,:,:), target ::        THX,RTHX,FTHX
     complex(r8),allocatable,dimension(:,:,:,:), target ::     CTHX,CRTHX,CFTHX
 
!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
! ADI variables
!----*|--.---------.---------.---------.---------.---------.---------.-|-------|

!     real(r8),allocatable,dimension(:,:) ::  MATL, MATD, &
!               MATU, VEC
!     real(r8),allocatable,dimension(:,:) ::  MATL_Z, MATD_Z, MATU_Z, VEC_Z

!      REAL(r8)     MATL (0:NX-1,0:NY+1), MATD(0:NX-1,0:NY+1), &
!          MATU(0:NX-1,0:NY+1), VEC(0:NX-1,0:NY+1)
!      REAL(r8)     MATL_Z (0:NX-1,0:NZ+1), MATD_Z(0:NX-1,0:NZ+1), &
!          MATU_Z(0:NX-1,0:NZ+1), VEC_Z(0:NX-1,0:NZ+1)

      real*8 TIME_OLD,TIME_OLD_ENG
     
      INTEGER NSAMPLES

      real(r8),allocatable,dimension(:,:) :: INT_JACOB
      real(r8),allocatable,dimension(:,:,:) :: GMAT_11,&
         GMAT_12, GMAT_22, &
         CJOB_11, CJOB_12, & 
         CJOB_21, CJOB_22

!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
! openMP variables - ALL ARRAYS
!----*|--.---------.---------.---------.---------.---------.---------.-|-------|

      integer kstart, kend, chunk,nthreads,iam,& !omp_get_num_threads, &
             k_start, k_end, j_start, j_end, jj_start, jj_end
      real(r8)  wtime, rktime !, omp_get_wtime
      common /bounds/ kstart, kend, k_start, k_end, j_start, j_end, &
      jj_start, jj_end

!      common /THOMAS_TH/ MATL, MATD, MATU, VEC
!!$omp threadprivate(/bounds/)
end module run_variable

module ADI_var
use ntypes, only : r8
use Domain

!      REAL(r8)     MATL (0:NX-1,0:NY+1), MATD(0:NX-1,0:NY+1), &
!          MATU(0:NX-1,0:NY+1), VEC(0:NX-1,0:NY+1)
!      REAL(r8)     MATL_Z (0:NX-1,0:NZ+1), MATD_Z(0:NX-1,0:NZ+1), &
!          MATU_Z(0:NX-1,0:NZ+1), VEC_Z(0:NX-1,0:NZ+1)
!     common /THOMAS_TH_Z/ MATL_Z, MATD_Z, MATU_Z, VEC_Z
!     common /THOMAS_TH/ MATL, MATD, MATU, VEC
!!$omp threadprivate(/THOMAS_TH_Z/, /THOMAS_TH/)
     REAL(r8),allocatable,DIMENSION(:,:) :: MATLX,MATDX,MATUX,VECX
     REAL(r8),allocatable,DIMENSION(:,:) :: MATLY,MATDY,MATUY,VECY
end module ADI_var

module variable_stat
use ntypes, only : r8
use Domain, only : N_TH
! Variables for outputting statistics
     real(r8),allocatable,dimension(:) :: UBAR,VBAR,WBAR
     real(r8),allocatable,dimension(:,:) :: URMS,VRMS, &
            WRMS, UV,UW, WV,PV,PU, &
            DWDY, DUDY,DWDZ, DUDZ, &
            SHEAR, OMEGA_X,OMEGA_Y,OMEGA_Z,TKE
         
     real(r8) ::     URMS_B,VRMS_B,WRMS_B,TKE_B

! Variables needed for SAVE_STATS_TH
     
      real(r8),allocatable,dimension(:,:) ::   THBAR, PE_DISS, Rig,tim

      real(r8) THRMS_B(1:N_TH)

      real(r8),allocatable,dimension(:,:,:) :: THRMS, &
                THV, THW, DTHDY, DTHDZ, DELTA_N2
      real(r8),allocatable,dimension(:,:,:) :: th_wh

! Variables for tkebudget
      real(r8),allocatable,dimension(:,:) :: epsilon, &
              tke_mean,tke_mean_old,energy_mean,energy_mean_old


      real(r8),allocatable,dimension(:,:) :: tke_1,tke_2, &
          tke_2_1, tke_2_2, tke_3, tke_3_1, tke_3_2, tke_3_3, &
          tke_4, tke_6_1,tke_6_1_1, tke_6_1_2, tke_6_1_3, &
          tke_6_2, tke_6_2_1, tke_6_2_2, tke_6_2_3, tke_7, &
          S1_mean, p_mean, transport, tke_5_1 !tke_5_1 is introduced by mainak
!for calculating u'rho'
 

      real(r8),allocatable,dimension(:,:,:) :: tke_5
      
end module variable_stat



module mg_vari
use ntypes
      
!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
! Multigrid variables - ALL ARRAYS 1D
!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      INTEGER ifail, istart, maxit, iprep, ndid, nout
      INTEGER bc(4), iout(6)
      LOGICAL INIT_FLAG, CALL_CHEEK_DIV 

      real(r8),allocatable,dimension(:,:) :: V, VC, &
             VB,VBC, RHS, RHSC, A, &
             AC, LDU, LDUC, WORK, WA, WAC, WB, WBC

      REAL*8 TOL, RESNO

end module mg_vari

module pr_rem
use ntypes
use Domain
REAL*8 P_TMP(0:NZ+1,0:NY+1), P_TMP2(0:NZ+1,0:NY+1), &
      RHS_TMP(0:NZ+1,0:NY+1)

      common /mg_rh/ P_TMP, P_TMP2,RHS_TMP
!$omp threadprivate(/mg_rh/)

end module pr_rem

module IO
 use ntypes, only: i4, r8
 integer(i4)             :: I_OUT,IOUT_MASTER,IOUT_SLAVE !Unit to write output 6 is screen
 character(len=100)      :: resultDIR,tempDIR,runDIR,ext,penDIR,plnDIR,statDIR,flowDIR, gridDIR,MGDIR,relaxDIR
end module IO


module les_chan_var
use ntypes
use Domain

! Variables for dynamic Smagrinsky
      real(r8),allocatable,dimension(:,:,:) :: numerator,denominator
      real(r8),allocatable,dimension(:,:) :: denominator_sum,numerator_sum

      real(r8),allocatable,dimension(:,:) :: NU_T_mean, DELTA_Y,U1_bar_les, &
                                             U3_bar_les,U2_bar_les,C_DYN
      real(r8),allocatable,dimension(:,:,:) :: KAPPA_T_mean,C_DYN_TH
      real(r8),allocatable,dimension(:,:,:) :: Sij_mean, TAU_mean
      real(r8),allocatable,dimension(:,:) :: tke_sgs_diss,tke_sgs_p,tke_sgs_t
      real(r8),allocatable,dimension(:,:,:,:) :: U_BAR_TIL,U_2BAR,U_4BAR
      real(r8),allocatable,dimension(:,:,:)   :: S_2BAR
      real(r8),allocatable,dimension(:,:,:,:)   :: St_rij 
      real(r8),allocatable,dimension(:,:,:)     :: GMAT_11_z,GMAT_11_y,GMAT_12_z,GMAT_12_y,GMAT_22_z,GMAT_22_y


      integer J1i, J2e
      

      REAL*8 C_DYN_H(0:NY+1),C_DYN_V(0:NY+1)
! Variables for plane-averaged momentum budget
      real*8 NU_U1(0:NY+1)
      real*8 NU_U3(0:NY+1)
      real*8 DELTA_YF(0:NY+1)

! For the TKE part
      real*8 tke_sgs_mm(0:NY+1),tke_sgs_evm(0:NY+1)

      REAL*8 TEMP(0:NXP,0:NZV-1,0:NY+1),Mij(0:NXP,0:NZV-1,0:NY+1)
      REAL*8 TEMP_1(0:NXP,0:NZV-1,0:NY+1)
      REAL*8 TEMP_2(0:NXP,0:NZV-1,0:NY+1)
      REAL*8 TEMP_3(0:NXP,0:NZV-1,0:NY+1)
      REAL*8  Sij(0:NXP,0:NZV-1,0:NY+1,1:6)

      REAL*8 cross

      COMPLEX*16 CSij(0:NX2P,0:NZV-1,0:NY+1,1:6)
      COMPLEX*16 CTEMP(0:NX2P,0:NZV-1,0:NY+1)
      COMPLEX*16 CTEMP_1(0:NX2P,0:NZV-1,0:NY+1)
      COMPLEX*16 CTEMP_2(0:NX2P,0:NZV-1,0:NY+1)
      COMPLEX*16 CTEMP_3(0:NX2P,0:NZV-1,0:NY+1)
      COMPLEX*16 CMij(0:NX2P,0:NZV-1,0:NY+1)

      EQUIVALENCE (Sij,CSij)    &
              ,(TEMP,CTEMP)     &
              ,(Mij,CMij)       &
              ,(TEMP_1,CTEMP_1) &
              ,(TEMP_2,CTEMP_2) &
              ,(TEMP_3,CTEMP_3) 

end module les_chan_var


module mpi_var
use ntypes
use Domain
INTEGER :: MPI_COMM_ROW, MPI_COMM_COL,MPI_COMM_CART
INTEGER :: NPROCES, RANK, RANK_ROW, RANK_COL
INTEGER :: status(MPI_STATUS_SIZE), IERROR


 integer(i4),parameter                   :: realtype  =MPI_DOUBLE_PRECISION
 integer(i4),parameter                   :: inttype   =MPI_INTEGER
 integer(i4),parameter                   :: chartype  =MPI_CHARACTER
 integer(i4),parameter                   :: logictype =MPI_LOGICAL
 integer(i4),parameter                   :: cmplxtype =MPI_DOUBLE_COMPLEX
 integer(i4),parameter                   :: commx1x2x3=MPI_COMM_WORLD

! complex(r8),allocatable,DIMENSION(:,:,:)    :: IN_CZ,IN_CX
! complex(r8),allocatable,DIMENSION(:,:,:)    :: OUT_CZ,OUT_CX

   
 real(r8),allocatable,DIMENSION(:)        :: M_IN,M_OUT
 complex(r8),allocatable,DIMENSION(:)     :: M_IN_C,M_OUT_C
end module mpi_var

module forcing
 use Domain
 use ntypes
      real(r8),allocatable,dimension(:,:) :: f, g
      real(r8),allocatable,dimension(:,:) :: f_prime, g_prime
      complex(r8),allocatable, dimension(:) ::  u_prime,v_prime, beam_th
      complex*8 im
!      real(r8) x0,y0,x_local,y_local, D, C 
end module forcing

module mpi_stitch
 use Domain
 use ntypes
     real(r8)    ::  U1_Y (0:NX+1,0:NZ+1), U2_Y(0:NX+1,0:NZ+1),U3_Y(0:NX+1,0:NZ+1)
     real(r8)    ::  TH_Y(0:NX+1,0:NZ+1,3)
 end module mpi_stitch
