! C******************************************************************************|
! C duct.f, the duct-flow solvers for diablo.                        VERSION 0.9
! C These solvers were written by ? and ? (spring 2001).
! C******************************************************************************|

! C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE INIT_DUCT
! C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
! C Initialize any constants here
use ntypes
use Domain
use Grid
use Fft_var
use TIME_STEP_VAR
use run_variable
use mg_vari, only : INIT_FLAG
use mpi_var
      
implicit none

      INTEGER J, K,N
      
      PI=4.D0*ATAN(1.D0)
! 
! ! At YMIN Location
! c        write(*,*) 'U_BC_YMIN: ',U_BC_YMIN
        IF (W_BC_YMIN.EQ.0) THEN
         IF (IC_TYPE == 7 ) THEN
          JSTART=1
         ELSE
          JSTART=2
         ENDIF
        ELSE IF (W_BC_YMIN.EQ.1) THEN
          JSTART=1
	ELSE IF (W_BC_YMIN.EQ.6) THEN
          JSTART=1  
        ELSE IF (W_BC_YMIN.EQ.7) THEN
          JSTART=1
        ELSE
          JSTART=2
        END IF
! ! At ZMIN Location  
! c        write(*,*) 'U_BC_ZMIN: ',U_BC_ZMIN
        IF (U_BC_ZMIN.EQ.0) THEN
	 IF (IC_TYPE == 7 ) THEN
          ZSTART=1
	 ELSE
	  ZSTART=2
	 ENDIF   
        ELSE IF (U_BC_ZMIN.EQ.1) THEN	
          ZSTART=1
	ELSE IF (U_BC_ZMIN.EQ.6) THEN	
          ZSTART=1  
        ELSE
          ZSTART=2
        END IF
! Now, set the indexing for the scalar equations
! At YMIN Location
        DO N=1,N_TH
          IF (TH_BC_YMIN(N).EQ.0) THEN
            IF (IC_TYPE == 7 ) THEN
             JSTART_TH(N)=1
            ELSE
             JSTART_TH(N)=2
            ENDIF 
          ELSE IF (TH_BC_YMIN(N).EQ.1) THEN
            JSTART_TH(N)=1
	  ELSE IF (TH_BC_YMIN(N).EQ.6) THEN
	    JSTART_TH(N)=1
          ELSE IF (TH_BC_YMIN(N).EQ.7) THEN
            JSTART_TH(N)=1
          ELSE IF (TH_BC_YMIN(N).EQ.8) THEN
            JSTART_TH(N)=1
          ELSE
            JSTART_TH(N)=2
          END IF
        END DO
! At ZMIN Location    
        DO N=1,N_TH
          IF (TH_BC_ZMIN(N).EQ.0) THEN
           IF (IC_TYPE == 7 ) THEN
            ZSTART_TH(N)=1
           ELSE
            ZSTART_TH(N)=2
           ENDIF
          ELSE IF (TH_BC_ZMIN(N).EQ.1) THEN
            ZSTART_TH(N)=1
	  ELSE IF (TH_BC_ZMIN(N).EQ.6) THEN
            ZSTART_TH(N)=1    
          ELSE IF (TH_BC_ZMIN(N).EQ.9) THEN
            ZSTART_TH(N)=1
          ELSE
            ZSTART_TH(N)=2
          END IF
        END DO      

! At YMAX Location

        IF (W_BC_YMAX.EQ.0) THEN
	 IF (IC_TYPE==7)THEN
	  JEND=NY-1
	 ELSE
          JEND=NY-1
	 ENDIF 
        ELSE IF (W_BC_YMAX.EQ.1) THEN
          JEND=NY
	ELSE IF (W_BC_YMAX.EQ.6) THEN
          JEND=NY  
        ELSE
          JEND=NY-1
        END IF

! At ZMAX Location
        IF (U_BC_ZMAX.EQ.0) THEN
	 IF (IC_TYPE == 7)THEN 
          ZEND=NZ
	 ELSE
	  ZEND=NZ-1
	 ENDIF 
        ELSE IF (U_BC_ZMAX.EQ.1) THEN
          ZEND=NZ
	ELSE IF (U_BC_ZMAX.EQ.6) THEN
          ZEND=NZ  
        ELSE
          ZEND=NZ-1
        END IF

! Set the upper and lower limits of timestepping of the scalar equations
! At YMAX Location 
        DO N=1,N_TH
        IF (TH_BC_YMAX(N).EQ.0) THEN
         IF (IC_TYPE == 7)THEN
          JEND_TH(N)=NY
         ELSE 
          JEND_TH(N)=NY-1
         ENDIF
        ELSE IF (TH_BC_YMAX(N).EQ.1) THEN
          JEND_TH(N)=NY
        ELSE IF (TH_BC_YMAX(N).EQ.6) THEN
	  JEND_TH(N)=NY
	ELSE IF (TH_BC_YMAX(N).EQ.5) THEN
          JEND_TH(N)=NY
        ELSE
          JEND_TH(N)=NY-1
        END IF
        END DO
! At ZMAX Location
        DO N=1,N_TH
        IF (TH_BC_ZMAX(N).EQ.0) THEN
         IF (IC_TYPE == 7 ) THEN
          ZEND_TH(N)=NZ
         ELSE
          ZEND_TH(N)=NZ-1
         ENDIF
        ELSE IF (TH_BC_ZMAX(N).EQ.1) THEN
          ZEND_TH(N)=NZ
        ELSE IF (TH_BC_ZMAX(N).EQ.5) THEN
          ZEND_TH(N)=NZ
	ELSE IF (TH_BC_ZMAX(N).EQ.6) THEN
          ZEND_TH(N)=NZ  
        ELSE
          ZEND_TH(N)=NZ-1
        END IF
        END DO

       IF (RANK .eq. 0 ) THEN  
       WRITE(6,*) 'CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC'
                  
       WRITE(6,*)'Boundary condition for U'
       WRITE(6,*)'U_BC_YMIN',U_BC_YMIN,'U_BC_YMAX',U_BC_YMAX
       WRITE(6,*)'U_BC_ZMIN',U_BC_ZMIN,'U_BC_ZMAX',U_BC_ZMAX
       WRITE(6,*) 'CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC'
                  
                  
       WRITE(6,*)'U_BC_YMIN_C1',U_BC_YMIN_C1,'U_BC_YMAX_C1',U_BC_YMAX_C1
       WRITE(6,*)'U_BC_ZMIN_C1',U_BC_ZMIN_C1,'U_BC_ZMAX_C1',U_BC_ZMAX_C1

       WRITE(6,*)'Boundary condition for W'
       WRITE(6,*)'W_BC_YMIN',W_BC_YMIN,'W_BC_YMAX', W_BC_YMAX
       WRITE(6,*)'W_BC_ZMIN',W_BC_ZMIN,'W_BC_ZMAX', W_BC_ZMAX
       WRITE(6,*) 'CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC'
                  
                  
       WRITE(6,*)'W_BC_YMIN_C1',W_BC_YMIN_C1,'W_BC_YMAX_C1',W_BC_YMAX_C1
       WRITE(6,*)'W_BC_ZMIN_C1',W_BC_ZMIN_C1,'W_BC_ZMAX_C1',W_BC_ZMAX_C1


       WRITE(6,*)'Boundary condition for V'
       WRITE(6,*)'V_BC_YMIN',V_BC_YMIN,'V_BC_YMAX',V_BC_YMAX
       WRITE(6,*)'V_BC_ZMIN',V_BC_ZMIN,'V_BC_ZMAX',V_BC_ZMAX
       WRITE(6,*) 'CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC'
                  
                   
       WRITE(6,*)'V_BC_YMIN_C1',V_BC_YMIN_C1,'V_BC_YMAX_C1',V_BC_YMAX_C1
       WRITE(6,*)'V_BC_ZMIN_C1',V_BC_ZMIN_C1,'V_BC_ZMAX_C1',V_BC_ZMAX_C1


       WRITE(6,*)'CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC'
                  
                 
       WRITE(6,*) 'JSATRT', JSTART, 'JEND', JEND
       WRITE(6,*) 'ZSATRT', ZSTART, 'ZEND', ZEND 
       WRITE(6,*)'CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC'
                  
                 
       DO N=1,N_TH
       WRITE(6,*)'Boundary condition for TH', N
       WRITE(6,*)'TH_BC_YMIN',TH_BC_YMIN(N),'TH_BC_YMAX', TH_BC_YMAX(N)
       WRITE(6,*)'TH_BC_ZMIN',TH_BC_ZMIN(N),'TH_BC_ZMAX', TH_BC_ZMAX(N)
       WRITE(6,*) 'CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC'
                  
                  
      WRITE(6,*)'TH_BC_YMIN_C1',TH_BC_YMIN_C1(N),   &
                'TH_BC_YMAX_C1',TH_BC_YMAX_C1(N)
      WRITE(6,*)'TH_BC_ZMIN_C1',TH_BC_ZMIN_C1(N),   &
                 'TH_BC_ZMAX_C1',TH_BC_ZMAX_C1(N)
      WRITE(6,*)'CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC'
                 
                 
       WRITE(6,*) 'JSATRT_TH', JSTART_TH(N), 'JEND_TH', JEND_TH(N)
       WRITE(6,*) 'ZSATRT_TH', ZSTART_TH(N), 'ZEND_TH', ZEND_TH(N) 
       WRITE(6,*)'CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC'
                  
                 
       ENDDO
       ENDIF
 

      RETURN
      END
! 
! C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE RK_DUCT_1
! C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
! C Main time-stepping algorithm for the duct-flow case.
! C INPUTS  (in Fourier space):  CUi, P, and (if k>1) CFi at (k-1)  (for i=1,2,3)
! C OUTPUTS (in Fourier space):  CUi, P, and (if k<3) CFi at (k)
! C----*|--.---------.---------.---------.---------.---------.---------.-|-------|

      RETURN
      END

! C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE RK_DUCT_2
! C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
! C Alternative time-stepping algorithm for the duct-flow case with ADI.
! C INPUTS  (in Fourier space):  CUi, P, and (if k>1) CFi at (k-1)  (for i=1,2,3)
! C OUTPUTS (in Fourier space):  CUi, P, and (if k<3) CFi at (k)
! C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
use ntypes
use Domain
use Grid
use Fft_var
use TIME_STEP_VAR
use run_variable
use ADI_var
use mg_vari, only : INIT_FLAG
use omp_lib      
use les_chan_var
use mpi_var
implicit none     

      INTEGER I,J,K,N      
      REAL*8 TEMP1, TEMP2, TEMP3, TEMP4, TEMP5, UBULK, &
            MEAN_U1,D
! C Communicate the information between ghost cells 
! 
! c      CALL GHOST_CHAN_MPI
! 
! C Define the constants that are used in the time-stepping
! C For reference, see Numerical Renaissance
      TEMP1=NU * H_BAR(RK_STEP) / 2.0
      TEMP2=H_BAR(RK_STEP) / 2.0
      TEMP3=ZETA_BAR(RK_STEP) * H_BAR(RK_STEP)
      TEMP4=H_BAR(RK_STEP)
      TEMP5=BETA_BAR(RK_STEP) * H_BAR(RK_STEP)

      
! C this is required convective boundary condition 
      dtc = TEMP4      
!      wtime =  omp_get_wtime ( )
!       
! C First, we will compute the explicit RHS terms and store in Ri
! C Note, Momentum equation and hence the RHS is evaluated at the
! C corresponding velocity points.

     
     CR1X(:,:,:) = (0.d0,0.d0)
     CR2X(:,:,:) = (0.d0,0.d0)
     CR3X(:,:,:) = (0.d0,0.d0)    
    

      DO J=JSTART,JEND
        DO K=ZSTART,ZEND
          DO I=0,NX2P
            CR1X(I,K,J)=INT_JACOB(K,J)*CU1X(I,K,J)
          END DO
        END DO
      END DO

      DO J= JSTART,JEND
        DO K=ZSTART,ZEND
          DO I=0,NX2P
            CR2X(I,K,J)=INT_JACOB(K,J)*CU2X(I,K,J)
          END DO
        END DO
      END DO

      DO J=JSTART,JEND
        DO K=ZSTART,ZEND
          DO I=0,NX2P
            CR3X(I,K,J)=INT_JACOB(K,J)*CU3X(I,K,J)
          END DO
        END DO
      END DO


! C Add the R-K term from the rk-1 step

      IF (RK_STEP .GT. 1) THEN
        DO J=JSTART,JEND
          DO K=ZSTART,ZEND
            DO I=0,NX2P
              CR1X(I,K,J)=CR1X(I,K,J)+TEMP3*CF1X(I,K,J)
            END DO
          END DO
        END DO
        DO J=JSTART,JEND
          DO K=ZSTART,ZEND
            DO I=0,NX2P
              CR2X(I,K,J)=CR2X(I,K,J)+TEMP3*CF2X(I,K,J)
            END DO
          END DO
        END DO
        DO J=JSTART,JEND
          DO K=ZSTART,ZEND
            DO I=0,NX2P
              CR3X(I,K,J)=CR3X(I,K,J)+TEMP3*CF3X(I,K,J)
            END DO
          END DO
        END DO
      END IF

      
! c       DO J=JSTART,JEND
! c         DO K=ZSTART,ZEND 
! c          CR1X(0,K,J)=CR1X(0,K,J)-TEMP4*INT_JACOB(K,J)*PX0
! c	  CR2X(0,K,J)=CR2X(0,K,J)-TEMP4*INT_JACOB(K,J)*PX0/2.0
! C	  CR3X(0,K,J)=CR3X(0,K,J)-TEMP4*INT_JACOB(K,J)*PX0
! c         END DO
! c        END DO   
! c      IF (F_TYPE == 1) THEN
! C Add the pressure gradient to the RHS as explicit Euler	
      DO J=JSTART,JEND
        DO K=ZSTART,ZEND
	 DO I=0,NX2P  
! C  dp/dz => d[J-1*\zeta_x*p]/d\zeta + d[J-1*\eta_x*p]/d\eta	 
	  CS1X(I,K,J) = 0.5*CJOB_11(K+1,J,2)*(CPX(I,K,J) + CPX(I,K+1,J)) &
     	               - 0.5*CJOB_11(K,J,2)*(CPX(I,K,J) + CPX(I,K-1,J)) &
                    + 0.5*CJOB_21(K,J+1,1)*(CPX(I,K,J) + CPX(I,K,J+1))  &
     	               - 0.5*CJOB_21(K,J,1)*(CPX(I,K,J) + CPX(I,K,J-1))  
! C  dp/dy => d[J-1*\zeta_y*p]/d\zeta + d[J-1*\eta_y*p]/d\eta          
          CS2X(I,K,J) = 0.5*CJOB_12(K+1,J,2)*(CPX(I,K,J) + CPX(I,K+1,J)) &
     	               - 0.5*CJOB_12(K,J,2)*(CPX(I,K,J) + CPX(I,K-1,J)) &
                     + 0.5*CJOB_22(K,J+1,1)*(CPX(I,K,J) + CPX(I,K,J+1)) &
      	               - 0.5*CJOB_22(K,J,1)*(CPX(I,K,J) + CPX(I,K,J-1))
        ENDDO
       END DO
      END DO	 
	 
! C Add the pressure gradient to the RHS as explicit Euler
      DO J=JSTART,JEND
        DO K=ZSTART,ZEND
          DO I=0,NX2P
! C  dp/dx => iKX(I)*J-1*{P}	  
            CR1X(I,K,J)=CR1X(I,K,J)-TEMP4*CIKXP(I) &
       	    *INT_JACOB(K,J)*CPX(I,K,J)
          END DO
        END DO
      END DO

      DO J=JSTART,JEND
        DO K=ZSTART,ZEND
          DO I=0,NX2P
            CR2X(I,K,J)=CR2X(I,K,J)-TEMP4*CS2X(I,K,J)
          END DO
        END DO
      END DO

      DO J=JSTART,JEND
        DO K=ZSTART,ZEND
          DO I=0,NX2P
            CR3X(I,K,J)=CR3X(I,K,J)-TEMP4*CS1X(I,K,J)
          END DO
        END DO
      END DO      
           	  	   
      IF (RANK .EQ. 0) THEN
      If ( F_TYPE == 1) THEN
       DO J=JSTART,JEND
         DO K=ZSTART,ZEND
             CR3X(0,K,J)=CR3X(0,K,J)-TEMP4*PX0*INT_JACOB(K,J)
         END DO
       END DO
      ELSE IF (F_TYPE.EQ.2) THEN
! C If oscillatory pressure gradient
        DO J=JSTART,JEND
         DO K=ZSTART,ZEND
          CR3X(0,K,J)=CR3X(0,K,J)-TEMP4*(PX0+AMP_OMEGA0  &
            *cos(OMEGA0*TIME))*INT_JACOB(K,J)
         ENDDO
        END DO 
      ENDIF
      ENDIF 
! C Now compute the term R-K term Ai
! C Compile terms of Ai in CFi which will be saved for next time step
! C First, store the horizontal viscous terms in CFi
      DO J=JSTART,JEND
        DO K=ZSTART,ZEND
          DO I=0,NX2P
            CF1X(I,K,J)=-NU * KX2P(I)*INT_JACOB(K,J)*CU1X(I,K,J) 
          END DO
        END DO
      END DO

      DO J=JSTART,JEND
        DO K=ZSTART,ZEND
          DO I=0,NX2P
            CF2X(I,K,J)=-NU * KX2P(I)*INT_JACOB(K,J)*CU2X(I,K,J) 
          END DO 
        END DO
      END DO

      DO J=JSTART,JEND
        DO K=ZSTART,ZEND
          DO I=0,NX2P
            CF3X(I,K,J)=-NU * KX2P(I)*INT_JACOB(K,J)* CU3X(I,K,J)
          END DO
        END DO
      END DO      
      
!  Adding rotation when velocities are in fourier space
      DO J=JSTART,JEND
        DO K=ZSTART,ZEND
          DO I=0,NX2P
            CF3X(I,K,J)=CF3X(I,K,J) + f_0*CU1X(I,K,J)*INT_JACOB(K,J)
            CF1X(I,K,J)=CF1X(I,K,J) - f_0*CU3X(I,K,J)*INT_JACOB(K,J)
          END DO
        END DO
      END DO     
 
! Do for each scalar
     IF (N_TH .eq. 1) THEN
! If a scalar contributes to the denisty, RI_TAU is not equal to zero and
! add the buoyancy term as explicit R-K.  Don't add the 0,0 mode which 
! corresponds to a plane average.  The plane averaged density balances
! the hydrostratic pressure component.

      DO J=JSTART,JEND
        DO K=ZSTART,ZEND
          DO I=0,NX2P
             CF2X(I,K,J)=CF2X(I,K,J) - INT_JACOB(K,J)*Ratio_gr* &
           CTHX(I,K,J,1) 
          END DO
        END DO
      END DO

      ELSEIF (N_TH .gt. 1) THEN

      IF (Non_linear_ST) THEN
       ! This term is added during the calculation of nonlinear terms bcz multiplication will be performed  
      ELSE
       DO J=JSTART,JEND
        DO K=ZSTART,ZEND
          DO I=0,NX2P
!             CF2X(I,K,J)=CF2X(I,K,J) + INT_JACOB(K,J)*Ratio_gr_a* &
!           CTHX(I,K,J,1) - INT_JACOB(K,J)*Ratio_gr_g*CTHX(I,K,J,2)
            CF2X(I,K,J)=CF2X(I,K,J) + (INT_JACOB(K,J)*Ratio_gr_a* &
           CTHX(I,K,J,1) - INT_JACOB(K,J)*Ratio_gr_g*CTHX(I,K,J,2))* &
           SIN(ANG_BETA)

           CF3X(I,K,J)=CF3X(I,K,J) - (INT_JACOB(K,J)*Ratio_gr_a* &
           CTHX(I,K,J,1) - INT_JACOB(K,J)*Ratio_gr_g*CTHX(I,K,J,2))* &
           COS(ANG_BETA)

          END DO
         END DO
       END DO
      ENDIF

      ENDIF      

      
! Do for each scalar
      DO N=1,N_TH
! Now, compute the RHS vector for the scalar equations
! Since TH is defined at horizontal velocity points, the
! scalar update equation will be very similar to the horizontal
! velocity update.

! We will store the RHS scalar terms in CRTH, RTHX
! The k-1 term for the R-K stepping is saved in FTHX, CFTHX

      CRTHX(:,:,:,N) = (0.d0,0.d0)

! First, build the RHS vector, use CRTH
      DO J=JSTART_TH(N),JEND_TH(N)
        DO K=ZSTART_TH(N),ZEND_TH(N)
          DO I=0,NX2P
          CRTHX(I,K,J,N)=INT_JACOB(K,J)*CTHX(I,K,J,N)
         ENDDO
       END DO
      END DO
! Add term from k-2 step to free up CFTHX variable
      IF (RK_STEP .GT. 1) THEN
        DO J=JSTART_TH(N),JEND_TH(N)
          DO K=ZSTART_TH(N),ZEND_TH(N)
            DO I=0,NX2P
              CRTHX(I,K,J,N)=CRTHX(I,K,J,N)+TEMP3*CFTHX(I,K,J,N)
            END DO
          END DO
        END DO
       END IF
   
       
! Now compute the explicit R-K term Ai
! Compile terms of Ai in CFi which will be saved for next time step
      DO J=JSTART_TH(N),JEND_TH(N)
        DO K=ZSTART_TH(N),ZEND_TH(N)
          DO I=0,NX2P
           CFTHX(I,K,J,N)=-(NU/PR(N))*INT_JACOB(K,J)*KX2P(I)*CTHX(I,K,J,N)
          END DO
        END DO
      END DO
      
!      IF ( CONT_STRAT ) THEN
!       DO J=JSTART_TH(N),JEND_TH(N)
!        DO K=ZSTART_TH(N),ZEND_TH(N)
!          DO I=0,NX2P
!            CFTHX(I,K,J,N)=CFTHX(I,K,J,N) + INT_JACOB(K,J)*  &
!                    CU2X(I,K,J)!/(1.d0 + SPONGE_SIGMA_OUT(K))
!          END DO
!        END DO
!       END DO
!      ELSE
       DO J=JSTART_TH(N),JEND_TH(N)
         DO K=ZSTART_TH(N),ZEND_TH(N)
           DO I=0,NX2P
              CFTHX(I,K,J,N)=CFTHX(I,K,J,N) -                               &
                       CU2X(I,K,J)*(                                        &
                   0.5*CJOB_12(K+1,J,2)*(TH_BAR(K,J,N) + TH_BAR(K+1,J,N))   &
                 - 0.5*CJOB_12(K,J,2)*(TH_BAR(K,J,N) + TH_BAR(K-1,J,N))     &
                 + 0.5*CJOB_22(K,J+1,1)*(TH_BAR(K,J,N) + TH_BAR(K,J+1,N))   &
                 - 0.5*CJOB_22(K,J,1)*(TH_BAR(K,J,N) + TH_BAR(K,J-1,N))  )

!!!!   Addition of kappa*d^2(rho_bar(x,z))/dz^2 + kappa*d^2(rho_bar(x,z))/dx^2
!!!!   at right hand side of density equation

               CFTHX(I,K,J,N)=CFTHX(I,K,J,N) +                                &
                    0.25*(NU/PR(N))*(GMAT_12(K+1,J,2)*                        &
                   (TH_BAR(K+1,J+1,N) + TH_BAR(K,J+1,N)                       &
                 - TH_BAR(K,J-1,N) - TH_BAR(K+1,J-1,N))                       &
                 - GMAT_12(K,J,2)*(TH_BAR(K-1,J+1,N) + TH_BAR(K,J+1,N)        &
                 - TH_BAR(K,J-1,N) - TH_BAR(K-1,J-1,N)))                      &
                 + 0.25*(NU/PR(N))*(GMAT_12(K,J+1,1)*                         &
                   (TH_BAR(K+1,J+1,N) + TH_BAR(K+1,J,N)                       &
                 - TH_BAR(K-1,J+1,N) - TH_BAR(K-1,J,N))                       &
                 -   GMAT_12(K,J,1)*(TH_BAR(K+1,J,N) + TH_BAR(K+1,J-1,N)      &
                 - TH_BAR(K-1,J,N) - TH_BAR(K-1,J-1,N)))                      &
                 + (NU/PR(N))                                                 &
                   *( GMAT_11(K+1,J,2)*(TH_BAR(K+1,J,N) - TH_BAR(K,J,N))      &
                 - GMAT_11(K,J,2)*(TH_BAR(K,J,N)   - TH_BAR(K-1,J,N)) )       &
                 + (NU/PR(N))                                                 &
                   *( GMAT_22(K,J+1,1)*(TH_BAR(K,J+1,N) - TH_BAR(K,J,N))      &
                 - GMAT_22(K,J,1)*(TH_BAR(K,J,N)- TH_BAR(K,J-1,N))  )
           END DO
         END DO
       END DO
!      ENDIF
! End of  loop for passive scalars (N_TH)

      CALL sponge_th(N)   

      IF (MELTING_MODEL)THEN
       dthdz_mean(:,n) = 0.0d0 ;
       IF (RANK ==0) THEN
        k=1
        do j=1,NY
         dthdz_mean(j,n)= dble(CTHX(0,k+1,j,n)-CTHX(0,k,j,n))
!                     (0.125*(CJOB_21(K,J,1)+CJOB_21(K,J+1,1) &
!                     + CJOB_21(K,J,2)+CJOB_21(K+1,J,2) )*    &
!                 dble(CTHX(0,k,j+1,n)-CTHX(0,k,j-1,n))       &
!               +  0.125*(CJOB_11(K,J,1)+CJOB_11(K,J+1,1)     &
!                     + CJOB_11(K,J,2)+CJOB_11(K+1,J,2) )*    &
!                  dble(CTHX(0,k+1,j,n)-CTHX(0,k-1,j,n)))/    &
!                  INT_JACOB(K,J)
         write(666,166) ypoint(1,j), dble(CU3X(0,1,j)),dble(CTHX(0,1,j,1)), dble(CTHX(0,1,j,2)), & 
                        dthdz_mean(j,1), dthdz_mean(j,2)
        enddo
        close(666) 
       ENDIF

       do j=0,NY+1
         CALL MPI_BCAST_REAL(dthdz_mean(j,n),1,1)
       enddo
      ENDIF 

      END DO
            
166 format (8f17.9)
      CALL sponge_vel
     


      IF (LES.AND.((.NOT.CREATE_NEW_FLOW).OR.(TIME_STEP.GT.1))) THEN

      call les_chan
! Add the subgrid scale scalar flux to the scalar equations
       DO N=1,N_TH
        call les_chan_th(N)
!         CALL FFT_X_TO_PHYSICAL(CTH(0,0,0,N),TH(0,0,0,N),0,NY+1,0,NZ+1)
       END DO


      ELSE 
! C If the subgrid model hasn't been called, then it is necessary to 
! C convert to physical space.

      CALL REAL_FOURIER_TRANS_U1 (.false.)
      CALL REAL_FOURIER_TRANS_U2 (.false.)
      CALL REAL_FOURIER_TRANS_U3 (.false.)

      
!C        CALL FFT_X_TO_PHYSICAL(CU1,U1,0,NY+1,0,NZ+1)
!C        CALL FFT_X_TO_PHYSICAL(CU2,U2,0,NY+1,0,NZ+1)
!C        CALL FFT_X_TO_PHYSICAL(CU3,U3,0,NY+1,0,NZ+1)

! Transform THETA to physical space for computation of nonlinear terms
! Here pass the first location in memory of the array for scalar n
      IF (N_TH .gt. 0) then
       CALL REAL_FOURIER_TRANS_TH (.false.)
      ENDIF
      
      END IF
! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!  Gravity terms are added here for  Nonlinear equation of state
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      S1X(:,:,:) =0.d0
      IF ( (N_TH .eq.2 ) .and. (Non_linear_ST) ) then
       call density_TC (.false.,.false.)
      ENDIF
      
! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! CCCCC        Non-linear terms calculation   CCCCC
! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! C actual nonlinear terms
! c      IF ( F_TYPE == 4 ) THEN
! C     For x-momentum equation
! C d(U1*u1)/dx
      S1X(:,:,:) =0.d0

      DO J=JSTART,JEND
        DO K=ZSTART,ZEND
          DO I=0,NXP
            S1X(I,K,J)=INT_JACOB(K,J)*U1X(I,K,J)*U1X(I,K,J)
          END DO
        END DO
      END DO

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

!C      CALL FFT_X_TO_FOURIER(S1,CS1,0,NY+1,0,NZ+1)
      
      DO J=JSTART,JEND
        DO K=ZSTART,ZEND
          DO I=0,NX2P
            CF1X(I,K,J)=CF1X(I,K,J) - CIKXP(I) * CS1X(I,K,J) 
          END DO
        END DO
      END DO

! C Calculation of Contravariant velocity

      DO J=JSTART,JEND+1
        DO K=ZSTART,ZEND+1
          DO I=0,NXP
! C      U2 required to define at bottom cell face 	  
            U2bX(I,K,J)=0.5*CJOB_22(K,J,1)*(U2X(I,K,J)+U2X(I,K,J-1))  &
                     + 0.5*CJOB_21(K,J,1)*(U3X(I,K,J)+U3X(I,K,J-1))
! C      U3 required to define at side cell face     	    
	    U3bX(I,K,J)=0.5*CJOB_11(K,J,2)*(U3X(I,K,J)+U3X(I,K-1,J)) &
     	              + 0.5*CJOB_12(K,J,2)*(U2X(I,K,J)+U2X(I,K-1,J))
          END DO
        END DO
      END DO

! C d(U2*u1)/d\eta and d(U3*u1)/d\zeta together       
      DO J=JSTART,JEND
        DO K=ZSTART,ZEND
          DO I=0,NXP
! C d(U2*u1)/d\eta	  
            S1X(I,K,J)=0.5*U2bX(I,K,J+1)*(U1X(I,K,J)+U1X(I,K,J+1)) &
       	    - 0.5*U2bX(I,K,J)*(U1X(I,K,J)+U1X(I,K,J-1))
! C d(U3*u1)/d\zeta   
            S2X(I,K,J)=0.5*U3bX(I,K+1,J)*(U1X(I,K,J)+U1X(I,K+1,J)) &
       	    - 0.5*U3bX(I,K,J)*(U1X(I,K,J)+U1X(I,K-1,J))              
          END DO
        END DO
      END DO
      

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

!C      CALL FFT_X_TO_FOURIER(S1,CS1,0,NY+1,0,NZ+1)
!C      CALL FFT_X_TO_FOURIER(S2,CS2,0,NY+1,0,NZ+1)
      
      
      DO J=JSTART,JEND
        DO K=ZSTART,ZEND
          DO I=0,NX2P
            CF1X(I,K,J)=CF1X(I,K,J) - CS1X(I,K,J) - CS2X(I,K,J)
          END DO
        END DO
      END DO
           
! C     For y-momentum equation
! C d(U1*u2)/dx

      DO J=JSTART,JEND
        DO K=ZSTART,ZEND
          DO I=0,NXP
            S1X(I,K,J)=INT_JACOB(K,J)*U1X(I,K,J)*U2X(I,K,J)
          END DO
        END DO
      END DO

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

!C      CALL FFT_X_TO_FOURIER(S1,CS1,0,NY+1,0,NZ+1)
      
      DO J=JSTART,JEND
        DO K=ZSTART,ZEND
          DO I=0,NX2P
            CF2X(I,K,J)=CF2X(I,K,J) - CIKXP(I) * CS1X(I,K,J) 
          END DO
        END DO
      END DO

! C d(U2*u2)/d\eta and d(U3*u2)/d\zeta together 
! 
! C Now at this time U2 at cell bottom and U3 at cell
! C    side face have been calculated previous step.  
      
      DO J=JSTART,JEND
        DO K=ZSTART,ZEND
          DO I=0,NXP
! C d(U2*u2)/d\eta	  
            S1X(I,K,J)=0.5*U2bX(I,K,J+1)*(U2X(I,K,J)+U2X(I,K,J+1)) &
       	    - 0.5*U2bX(I,K,J)*(U2X(I,K,J)+U2X(I,K,J-1))
! C d(U3*u2)/d\zeta   
            S2X(I,K,J)=0.5*U3bX(I,K+1,J)*(U2X(I,K,J)+U2X(I,K+1,J)) &
       	    - 0.5*U3bX(I,K,J)*(U2X(I,K,J)+U2X(I,K-1,J))              
          END DO
        END DO
      END DO

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
      
!C       CALL FFT_X_TO_FOURIER(S1,CS1,0,NY+1,0,NZ+1)
!C       CALL FFT_X_TO_FOURIER(S2,CS2,0,NY+1,0,NZ+1)
      
      
      DO J=JSTART,JEND
        DO K=ZSTART,ZEND
          DO I=0,NX2P
            CF2X(I,K,J)=CF2X(I,K,J) - CS1X(I,K,J) - CS2X(I,K,J)
          END DO
        END DO
      END DO      

! C     For z-momentum equation
! C d(U1*u3)/dx

      DO J=JSTART,JEND
        DO K=ZSTART,ZEND
          DO I=0,NXP
            S1X(I,K,J)=INT_JACOB(K,J)*U1X(I,K,J)*U3X(I,K,J)
          END DO
        END DO
      END DO
      
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

!C      CALL FFT_X_TO_FOURIER(S1,CS1,0,NY+1,0,NZ+1)
      
      DO J=JSTART,JEND
        DO K=ZSTART,ZEND
          DO I=0,NX2P
            CF3X(I,K,J)=CF3X(I,K,J) - CIKXP(I) * CS1X(I,K,J) 
          END DO
        END DO
      END DO

! C d(U2*u3)/d\eta and d(U3*u3)/d\zeta together 
! 
! C Now at this time U2 at cell bottom and U3 at cell
! C    side face have been calculated previous step.  
      
      DO J=JSTART,JEND
        DO K=ZSTART,ZEND
          DO I=0,NXP
! C d(U2*u3)/d\eta	  
            S1X(I,K,J)=0.5*U2bX(I,K,J+1)*(U3X(I,K,J)+U3X(I,K,J+1)) &
       	    - 0.5*U2bX(I,K,J)*(U3X(I,K,J)+U3X(I,K,J-1))
! C d(U3*u3)/d\zeta   
            S2X(I,K,J)=0.5*U3bX(I,K+1,J)*(U3X(I,K,J)+U3X(I,K+1,J)) &
       	    - 0.5*U3bX(I,K,J)*(U3X(I,K,J)+U3X(I,K-1,J))              
          END DO
        END DO
      END DO
  
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
    
!C      CALL FFT_X_TO_FOURIER(S1,CS1,0,NY+1,0,NZ+1)
!C      CALL FFT_X_TO_FOURIER(S2,CS2,0,NY+1,0,NZ+1)
      
      
      DO J=JSTART,JEND
        DO K=ZSTART,ZEND
          DO I=0,NX2P
            CF3X(I,K,J)=CF3X(I,K,J) - CS1X(I,K,J) - CS2X(I,K,J)
          END DO
        END DO
      END DO      
      
! c      ENDIF      
! ! No-linear terms for U1 comes from the diffusion terms
! 
! ! Kappa*d/d\zeta(GMAT_12(:,:,2)dU1/d\eta

      DO J=JSTART,JEND
        DO K=ZSTART,ZEND
          DO I=0,NXP
         S1X(I,K,J)=  0.25*NU*(GMAT_12(K+1,J,2)*  &
     	          (U1X(I,K+1,J+1) + U1X(I,K,J+1)   &
                - U1X(I,K,J-1) - U1X(I,K+1,J-1))  &
             -  GMAT_12(K,J,2)*(U1X(I,K-1,J+1) + U1X(I,K,J+1) &
                - U1X(I,K,J-1) - U1X(I,K-1,J-1)))
          END DO
        END DO
      END DO

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

      
!C      CALL FFT_X_TO_FOURIER(S1,CS1,0,NY+1,0,NZ+1)

      DO J=JSTART,JEND
        DO K=ZSTART,ZEND
          DO I=0,NX2P
            CF1X(I,K,J)=CF1X(I,K,J) + CS1X(I,K,J)
          END DO
        END DO
      END DO             
         
! Kappa*d/d\eta(GMAT_12(:,:,1))dU1/d\zeta

      DO J=JSTART,JEND
        DO K=ZSTART,ZEND
          DO I=0,NXP
         S1X(I,K,J)=  0.25*NU*(GMAT_12(K,J+1,1)* &
     	          (U1X(I,K+1,J+1) + U1X(I,K+1,J)  &
                - U1X(I,K-1,J+1) - U1X(I,K-1,J)) &
                -   GMAT_12(K,J,1)*(U1X(I,K+1,J) + U1X(I,K+1,J-1) &
                - U1X(I,K-1,J) - U1X(I,K-1,J-1)))
          END DO
        END DO
      END DO

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

!C      CALL FFT_X_TO_FOURIER(S1,CS1,0,NY+1,0,NZ+1)

      DO J=JSTART,JEND
        DO K=ZSTART,ZEND
          DO I=0,NX2P
            CF1X(I,K,J)=CF1X(I,K,J) + CS1X(I,K,J)
          END DO
        END DO
      END DO

      
! No-linear terms for U2 comes from the diffusion terms

! Kappa*d/d\zeta(GMAT_12(:,:,2)dU2/d\eta

      DO J=JSTART,JEND
        DO K=ZSTART,ZEND
          DO I=0,NXP
          S1X(I,K,J)=  0.25*NU*(GMAT_12(K+1,J,2)*  &
     	          (U2X(I,K+1,J+1) + U2X(I,K,J+1)   &
                - U2X(I,K,J-1) - U2X(I,K+1,J-1))   &
             -  GMAT_12(K,J,2)*(U2X(I,K-1,J+1) + U2X(I,K,J+1) &
                - U2X(I,K,J-1) - U2X(I,K-1,J-1)))
          END DO
        END DO
      END DO

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

!C      CALL FFT_X_TO_FOURIER(S1,CS1,0,NY+1,0,NZ+1)

      DO J=JSTART,JEND
        DO K=ZSTART,ZEND
          DO I=0,NX2P
            CF2X(I,K,J)=CF2X(I,K,J) + CS1X(I,K,J)
          END DO
        END DO
      END DO             
         
! Kappa*d/d\eta(GMAT_12(:,:,1))dU2/d\zeta

      DO J=JSTART,JEND
        DO K=ZSTART,ZEND
          DO I=0,NXP
         S1X(I,K,J)=  0.25*NU*(GMAT_12(K,J+1,1)* &
     	          (U2X(I,K+1,J+1) + U2X(I,K+1,J)  &
                - U2X(I,K-1,J+1) - U2X(I,K-1,J)) &
                -   GMAT_12(K,J,1)*(U2X(I,K+1,J) + U2X(I,K+1,J-1) &
                - U2X(I,K-1,J) - U2X(I,K-1,J-1)))
          END DO
        END DO
      END DO

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

!C      CALL FFT_X_TO_FOURIER(S1,CS1,0,NY+1,0,NZ+1)

      DO J=JSTART,JEND
        DO K=ZSTART,ZEND
          DO I=0,NX2P
            CF2X(I,K,J)=CF2X(I,K,J) + CS1X(I,K,J)
          END DO
        END DO
      END DO      
            

! No-linear terms for U3 comes from the diffusion terms

! Kappa*d/d\zeta(GMAT_12(:,:,2)dU3/d\eta

      DO J=JSTART,JEND
        DO K=ZSTART,ZEND
          DO I=0,NXP
         S1X(I,K,J)=  0.25*NU*(GMAT_12(K+1,J,2)*  &
     	          (U3X(I,K+1,J+1) + U3X(I,K,J+1)   &
                - U3X(I,K,J-1) - U3X(I,K+1,J-1))  &
             -  GMAT_12(K,J,2)*(U3X(I,K-1,J+1) + U3X(I,K,J+1) &
                - U3X(I,K,J-1) - U3X(I,K-1,J-1)))
          END DO
        END DO
      END DO

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

!C      CALL FFT_X_TO_FOURIER(S1,CS1,0,NY+1,0,NZ+1)
      DO J=JSTART,JEND
        DO K=ZSTART,ZEND
          DO I=0,NX2P
            CF3X(I,K,J)=CF3X(I,K,J) + CS1X(I,K,J)
          END DO
        END DO
      END DO             
         
! Kappa*d/d\eta(GMAT_12(:,:,1))dU3/d\zeta

      DO J=JSTART,JEND
        DO K=ZSTART,ZEND
          DO I=0,NXP
         S1X(I,K,J)=  0.25*NU*(GMAT_12(K,J+1,1)*  &
     	          (U3X(I,K+1,J+1) + U3X(I,K+1,J)   &
                - U3X(I,K-1,J+1) - U3X(I,K-1,J))  &
                - GMAT_12(K,J,1)*(U3X(I,K+1,J)+U3X(I,K+1,J-1) &
                - U3X(I,K-1,J) - U3X(I,K-1,J-1)))
          END DO
        END DO
      END DO

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

!C      CALL FFT_X_TO_FOURIER(S1,CS1,0,NY+1,0,NZ+1)

      DO J=JSTART,JEND
        DO K=ZSTART,ZEND
          DO I=0,NX2P
            CF3X(I,K,J)=CF3X(I,K,J) + CS1X(I,K,J)
          END DO
        END DO
      END DO       



      IF (LES.AND.((.NOT.CREATE_NEW_FLOW).OR.(TIME_STEP.GT.1))) THEN
       DO J=JSTART,JEND
        DO K=ZSTART,ZEND
          DO I=0,NXP
           S1X(I,K,J)=   0.25*(GMAT_12(K+1,J,2)*0.5d0*(NU_T(I,K+1,J)+NU_T(I,K,J))*   &
                          (U1X(I,K+1,J+1) + U1X(I,K,J+1)- U1X(I,K,J-1) - U1X(I,K+1,J-1))&
                        -  GMAT_12(K,J,2)*0.5d0*(NU_T(I,K,J)+NU_T(I,K-1,J))*        &
                          (U1X(I,K-1,J+1) + U1X(I,K,J+1) -  U1X(I,K,J-1) - U1X(I,K-1,J-1)) ) &
                      + 0.25*(GMAT_12(K,J+1,1)*0.5d0*(NU_T(I,K,J)+NU_T(I,K,J+1))*     &
                          (U1X(I,K+1,J+1) + U1X(I,K+1,J) - U1X(I,K-1,J+1) - U1X(I,K-1,J)) &
                         - GMAT_12(K,J,1)*0.5d0*(NU_T(I,K,J)+NU_T(I,K,J-1))*        &
                          (U1X(I,K+1,J) + U1X(I,K+1,J-1) - U1X(I,K-1,J) - U1X(I,K-1,J-1)))
        END DO  
       END DO
      END DO

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

!C      CALL FFT_X_TO_FOURIER(S1,CS1,0,NY+1,0,NZ+1)

      DO J=JSTART,JEND
        DO K=ZSTART,ZEND
          DO I=0,NX2P
            CF1X(I,K,J)=CF1X(I,K,J) + CS1X(I,K,J)
          END DO
        END DO
      END DO

      DO J=JSTART,JEND
        DO K=ZSTART,ZEND
          DO I=0,NXP
           S1X(I,K,J)=   0.25*(GMAT_12_y(K+1,J,2)*0.5d0*(NU_T(I,K+1,J)+NU_T(I,K,J))*   &
                          (U2X(I,K+1,J+1) + U2X(I,K,J+1)- U2X(I,K,J-1) - U2X(I,K+1,J-1))  &
                        -  GMAT_12_y(K,J,2)*0.5d0*(NU_T(I,K,J)+NU_T(I,K-1,J))*        &
                          (U2X(I,K-1,J+1) + U2X(I,K,J+1) -  U2X(I,K,J-1) - U2X(I,K-1,J-1)) ) &
                       + 0.25*(GMAT_12_y(K,J+1,1)*0.5d0*(NU_T(I,K,J)+NU_T(I,K,J+1))*     &
                          (U2X(I,K+1,J+1) + U2X(I,K+1,J) - U2X(I,K-1,J+1) - U2X(I,K-1,J)) &
                         - GMAT_12_y(K,J,1)*0.5d0*(NU_T(I,K,J)+NU_T(I,K,J-1))*        &
                          (U2X(I,K+1,J) + U2X(I,K+1,J-1) - U2X(I,K-1,J) - U2X(I,K-1,J-1)))

          END DO
        END DO
       END DO

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

!C       CALL FFT_X_TO_FOURIER(S1,CS1,0,NY+1,0,NZ+1)

      DO J=JSTART,JEND
        DO K=ZSTART,ZEND
          DO I=0,NX2P
            CF2X(I,K,J)=CF2X(I,K,J) + CS1X(I,K,J)
          END DO
        END DO
      END DO
                
      DO J=JSTART,JEND
        DO K=ZSTART,ZEND
          DO I=0,NXP
           S1X(I,K,J)=    0.25*(GMAT_12_z(K+1,J,2)*0.5d0*(NU_T(I,K+1,J)+NU_T(I,K,J))*   &
                          (U3X(I,K+1,J+1) + U3X(I,K,J+1)- U3X(I,K,J-1) - U3X(I,K+1,J-1))  &
                        -  GMAT_12_z(K,J,2)*0.5d0*(NU_T(I,K,J)+NU_T(I,K-1,J))*        &
                          (U3X(I,K-1,J+1) + U3X(I,K,J+1) -  U3X(I,K,J-1) - U3X(I,K-1,J-1)) ) &
                       + 0.25*(GMAT_12_z(K,J+1,1)*0.5d0*(NU_T(I,K,J)+NU_T(I,K,J+1))*     &
                          (U3X(I,K+1,J+1) + U3X(I,K+1,J) - U3X(I,K-1,J+1) - U3X(I,K-1,J)) &
                         - GMAT_12_z(K,J,1)*0.5d0*(NU_T(I,K,J)+NU_T(I,K,J-1))*        &
                          (U3X(I,K+1,J) + U3X(I,K+1,J-1) - U3X(I,K-1,J) - U3X(I,K-1,J-1)))
            END DO
          END DO
        END DO

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

!C       CALL FFT_X_TO_FOURIER(S1,CS1,0,NY+1,0,NZ+1)

      DO J=JSTART,JEND
        DO K=ZSTART,ZEND
          DO I=0,NX2P
            CF3X(I,K,J)=CF3X(I,K,J) + CS1X(I,K,J)
          END DO
        END DO
      END DO


      ENDIF


            
! C -- At this point, we are done computing the nonlinear terms --
! 
!       
! C Finally, Add CFi to CRi

      DO J=JSTART,JEND
        DO K=ZSTART,ZEND
          DO I=0,NX2P
            CR1X(I,K,J)=CR1X(I,K,J) + TEMP5 * CF1X(I,K,J)
          END DO
        END DO
      END DO
      DO J=JSTART,JEND
        DO K=ZSTART,ZEND
          DO I=0,NX2P
            CR2X(I,K,J)=CR2X(I,K,J) + TEMP5 * CF2X(I,K,J)
          END DO
        END DO
      END DO
 
      DO J=JSTART,JEND
        DO K=ZSTART,ZEND
          DO I=0,NX2P
            CR3X(I,K,J)=CR3X(I,K,J) + TEMP5 * CF3X(I,K,J)
          END DO
        END DO
      END DO

! C Convert RHS terms to physical space
! C made a change in FFT_X_TO_PHYSICAL(CR2X,R2X,0,NY+1,0,NZ+1)
! C
! 
! c      IF (MOD(TIME_STEP,5).EQ.0) THEN 
! c       CALL sponge_vel
! c      endif

      CALL REAL_FOURIER_TRANS_R1 (.false.)
      CALL REAL_FOURIER_TRANS_R2 (.false.)
      CALL REAL_FOURIER_TRANS_R3 (.false.)

!      CALL FFT_X_TO_PHYSICAL(CR1X,R1X,0,NY+1,0,NZ+1)                 
!      CALL FFT_X_TO_PHYSICAL(CR2X,R2X,0,NY+1,0,NZ+1)                 
!      CALL FFT_X_TO_PHYSICAL(CR3X,R3X,0,NY+1,0,NZ+1)      
      
          
      DO N=1,N_TH
! ! Do for each scalar:
! 
! c      IF ( F_TYPE .EQ. 4 ) THEN
! ! Compute the nonlinear terms that are present in the explicit term A
! U1*TH
      DO J=JSTART_TH(N),JEND_TH(N)
        DO K=ZSTART_TH(N),ZEND_TH(N)
          DO I=0,NXP
            S1X(I,K,J)=INT_JACOB(K,J)*THX(I,K,J,N)*U1X(I,K,J)
          END DO
        END DO
      END DO

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

!C      CALL FFT_X_TO_FOURIER(S1,CS1,0,NY+1,0,NZ+1)

      DO J=JSTART_TH(N),JEND_TH(N)
        DO K=ZSTART_TH(N),ZEND_TH(N)
          DO I=0,NX2P
            CFTHX(I,K,J,N)=CFTHX(I,K,J,N) - CIKXP(I) * CS1X(I,K,J)
          END DO
        END DO
      END DO
! 
! ! U3*TH and U2*TH together
! C d(U2*TH)/d\eta and d(U3*TH)/d\zeta together       
      DO J=JSTART_TH(N),JEND_TH(N)
        DO K=ZSTART_TH(N),ZEND_TH(N)
          DO I=0,NXP
! C d(U2*TH)/d\eta	  
            S1X(I,K,J)=0.5*U2bX(I,K,J+1)*(THX(I,K,J,N)+THX(I,K,J+1,N)) &
       	    - 0.5*U2bX(I,K,J)*(THX(I,K,J,N)+THX(I,K,J-1,N))
! C d(U3*TH)/d\zeta   
            S2X(I,K,J)=0.5*U3bX(I,K+1,J)*(THX(I,K,J,N)+THX(I,K+1,J,N)) &
       	    - 0.5*U3bX(I,K,J)*(THX(I,K,J,N)+THX(I,K-1,J,N))              
          END DO
        END DO
      END DO
  

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
    
!C      CALL FFT_X_TO_FOURIER(S1,CS1,0,NY+1,0,NZ+1)
!C      CALL FFT_X_TO_FOURIER(S2,CS2,0,NY+1,0,NZ+1)
      
      
      DO J=JSTART_TH(N),JEND_TH(N)
        DO K=ZSTART_TH(N),ZEND_TH(N)
          DO I=0,NX2P
            CFTHX(I,K,J,N)=CFTHX(I,K,J,N) - CS1X(I,K,J) - CS2X(I,K,J)
          END DO
        END DO
      END DO
! Done with non-linear terms for theta

! No-linear terms comes from the diffusion terms

! Kappa*d/d\zeta(GMAT_12(:,:,2)dTH/d\eta

      DO J=JSTART_TH(N),JEND_TH(N)
        DO K=ZSTART_TH(N),ZEND_TH(N)
          DO I=0,NXP
         S1X(I,K,J)=  0.25*(NU/PR(N))*(GMAT_12(K+1,J,2)* &
     	          (THX(I,K+1,J+1,N) + THX(I,K,J+1,N)      &
                - THX(I,K,J-1,N) - THX(I,K+1,J-1,N))     &
                -   GMAT_12(K,J,2)*(THX(I,K-1,J+1,N) + THX(I,K,J+1,N) &
                - THX(I,K,J-1,N) - THX(I,K-1,J-1,N)))
          END DO
        END DO
      END DO

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

      

!C      CALL FFT_X_TO_FOURIER(S1,CS1,0,NY+1,0,NZ+1)

      DO J=JSTART_TH(N),JEND_TH(N)
        DO K=ZSTART_TH(N),ZEND_TH(N)
          DO I=0,NX2P
            CFTHX(I,K,J,N)=CFTHX(I,K,J,N) + CS1X(I,K,J)
          END DO
        END DO
      END DO             
         
! Kappa*d/d\eta(GMAT_12(:,:,1))dTH/d\zeta

      DO J=JSTART_TH(N),JEND_TH(N)
        DO K=ZSTART_TH(N),ZEND_TH(N)
          DO I=0,NXP
           S1X(I,K,J)=  0.25*(NU/PR(N))*(GMAT_12(K,J+1,1)* &
     	             (THX(I,K+1,J+1,N) + THX(I,K+1,J,N)      &
                   - THX(I,K-1,J+1,N) - THX(I,K-1,J,N))     &
                   -   GMAT_12(K,J,1)*(THX(I,K+1,J,N) + THX(I,K+1,J-1,N) &
                   - THX(I,K-1,J,N) - THX(I,K-1,J-1,N)))
          END DO
        END DO
      END DO

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
 
!C      CALL FFT_X_TO_FOURIER(S1,CS1,0,NY+1,0,NZ+1)

      DO J=JSTART_TH(N),JEND_TH(N)
        DO K=ZSTART_TH(N),ZEND_TH(N)
          DO I=0,NX2P
            CFTHX(I,K,J,N)=CFTHX(I,K,J,N) + CS1X(I,K,J)
          END DO
        END DO
      END DO      
      
! We are done with the horizontal derivatives of the nonlinear terms

      IF (LES.AND.((.NOT.CREATE_NEW_FLOW).OR.(TIME_STEP.GT.1))) THEN

       DO J=JSTART_TH(N),JEND_TH(N)
        DO K=ZSTART_TH(N),ZEND_TH(N)
          DO I=0,NXP
           S1X(I,K,J)=   0.25*(GMAT_12(K+1,J,2)*0.5d0*(KAPPA_T(I,K+1,J,N)+KAPPA_T(I,K,J,N))*   &
                          (THX(I,K+1,J+1,N) + THX(I,K,J+1,N)- THX(I,K,J-1,N) - THX(I,K+1,J-1,N))&
                        -  GMAT_12(K,J,2)*0.5d0*(KAPPA_T(I,K,J,N)+KAPPA_T(I,K-1,J,N))*        &
                          (THX(I,K-1,J+1,N) + THX(I,K,J+1,N) -  THX(I,K,J-1,N) - THX(I,K-1,J-1,N)) ) &
                      + 0.25*(GMAT_12(K,J+1,1)*0.5d0*(KAPPA_T(I,K,J,N)+KAPPA_T(I,K,J+1,N))*     &
                          (THX(I,K+1,J+1,N) + THX(I,K+1,J,N) - THX(I,K-1,J+1,N) - THX(I,K-1,J,N)) &
                         - GMAT_12(K,J,1)*0.5d0*(KAPPA_T(I,K,J,N)+KAPPA_T(I,K,J-1,N))*        &
                          (THX(I,K+1,J,N) + THX(I,K+1,J-1,N) - THX(I,K-1,J,N) - THX(I,K-1,J-1,N)))
        END DO
       END DO
      END DO

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
 

!C      CALL FFT_X_TO_FOURIER(S1,CS1,0,NY+1,0,NZ+1)

      DO J=JSTART_TH(N),JEND_TH(N)
        DO K=ZSTART_TH(N),ZEND_TH(N)
          DO I=0,NX2P
            CFTHX(I,K,J,N)=CFTHX(I,K,J,N) + CS1X(I,K,J)
          END DO
        END DO
      END DO

      ENDIF



! Add CFTHX to the RHS vector CRTH
      DO J=JSTART_TH(N),JEND_TH(N)
        DO K=ZSTART_TH(N),ZEND_TH(N)
          DO I=0,NX2P
            CRTHX(I,K,J,N)=CRTHX(I,K,J,N) + TEMP5 * CFTHX(I,K,J,N)
          END DO
        END DO
      END DO
! Done with computation of RHS, explicit terms for the THETA equation
! Transform back to physical space
      END DO   ! end of the scalar variables 


        CALL REAL_FOURIER_TRANS_Rth (.false.)

!C      CALL FFT_X_TO_PHYSICAL(CRTH(0,0,0,N),RTHX(0,0,0,N),0,NY+1,0,NZ+1)

        DO N=1,N_TH
! ! Do for each scalar:
     
! 
! C Before going to  ADI1 we have to store RHS after subtracting previous TH_n 
! C at previous time step to bulid new RHS for  ADI2 based on TH_n+1/2 at intermidiate 



      DO J=JSTART_TH(N),JEND_TH(N)
        DO K=ZSTART_TH(N),ZEND_TH(N)
          DO I=0,NXP
            S1X(I,K,J)= RTHX(I,K,J,N)-INT_JACOB(K,J)*THX(I,K,J,N)
          END DO
        END DO
      END DO


! C Compute the vertical viscous term in physical space and add to RHS
! C  at  ADI1 step .... Important in this step we have to add TH_n
! C at previous time step so that after multiplying a factor 1/2 into the 
! C RHS side during ADI1 operation RHS will be taken care TH_n not 1/2*TH_n


      IF (WAVE_ABS) THEN
      DO J=JSTART_TH(N),JEND_TH(N)
        DO K=ZSTART_TH(N),ZEND_TH(N)
          DO I=0,NXP
            RTHX(I,K,J,N)=0.5*S1X(I,K,J) + INT_JACOB(K,J)*THX(I,K,J,N) &
            +  (TEMP1/PR(N))*(1.d0 + 5.d0*SPONGE_SIGMA_OUT(K,J))     &
              *( GMAT_11(K+1,J,2)*(THX(I,K+1,J,N) - THX(I,K,J,N))    &
                -GMAT_11(K,J,2)*(THX(I,K,J,N)   - THX(I,K-1,J,N)) )  
          END DO
        END DO
      END DO
      ELSE 

      DO J=JSTART_TH(N),JEND_TH(N)
        DO K=ZSTART_TH(N),ZEND_TH(N)
          DO I=0,NXP
            RTHX(I,K,J,N)=0.5*S1X(I,K,J) + INT_JACOB(K,J)*THX(I,K,J,N)  &
            +  (TEMP1/PR(N)) &
              *( GMAT_11(K+1,J,2)*(THX(I,K+1,J,N) - THX(I,K,J,N))    &
                -GMAT_11(K,J,2)*(THX(I,K,J,N)   - THX(I,K-1,J,N)) )
          END DO
        END DO
      END DO
      ENDIF

      IF (LES.AND.((.NOT.CREATE_NEW_FLOW).OR.(TIME_STEP.GT.1))) THEN
       DO J=JSTART_TH(N),JEND_TH(N)
        DO K=ZSTART_TH(N),ZEND_TH(N)
          DO I=0,NXP
            RTHX(I,K,J,N)=RTHX(I,K,J,N) + TEMP2*(0.5d0*(KAPPA_T(I,K+1,J,N)+KAPPA_T(I,K,J,N)) &
                    *(GMAT_11(K+1,J,2)*(THX(I,K+1,J,N) - THX(I,K,J,N)))                      &
                    -0.5d0*(KAPPA_T(I,K,J,N) + KAPPA_T(I,K-1,J,N))                         &
                    *(GMAT_11(K,J,2)*(THX(I,K,J,N)  - THX(I,K-1,J,N))))
         END DO
        END DO
       END DO
      ENDIF


! c      IF ( TH_BC_ZMAX(N)  .EQ. 1 ) THEN
! c       DO J=JSTART_TH(N),JEND_TH(N)
! c	   DO I=0,NXM
! c             RTHX(I,NZ+1,J,N)=RTHX(I,NZ,J,N)
! c       	   END DO
! c       END DO  
! c      ENDIF
!  
! c      IF ( TH_BC_ZMIN(N)  .EQ. 1 ) THEN
! c       DO J=JSTART_TH(N),JEND_TH(N)
! c	   DO I=0,NXM
! c             RTHX(I,0,J,N)=RTHX(I,1,J,N)
! c       	   END DO
! c       END DO  
! c      ENDIF

      
! C Solve for TH
! C Note, here the matrix will be indexed from 1...NY+1 corresponding to U1(0:NY)

! Initialize the matrix used to store implicit coefficients
      DO J=0,NY+1
        DO I=0,NXP
          MATLY(I,J)=0.
          MATDY(I,J)=1.
          MATUY(I,J)=0.
          VECY(I,J)=0.
        END DO
      END DO 

! C Build the implicit system of equations for U1 
!$OMP PARALLEL DO PRIVATE(K,J,I)
      DO K=ZSTART_TH(N),ZEND_TH(N)

        DO J=JSTART_TH(N),JEND_TH(N)
          DO I=0,NXP
            MATLY(I,J)=-(TEMP1/PR(N))*GMAT_22(K,J,1)
            MATDY(I,J)= INT_JACOB(K,J)-(TEMP1/PR(N))* &
              (-GMAT_22(K,J+1,1) - GMAT_22(K,J,1))   
            MATUY(I,J)=-(TEMP1/PR(N))*GMAT_22(K,J+1,1)
            VECY(I,J)=RTHX(I,K,J,N)
          END DO
        END DO

! IF using a subgrid model (LES) then add the eddy viscosity part implicitly

         IF (LES.AND.((.NOT.CREATE_NEW_FLOW).OR.(TIME_STEP.GT.1))) THEN

          DO J=JSTART_TH(N),JEND_TH(N)
          DO I=0,NXP
              MATLY(I,J) = MATLY(I,J)                                                  &
           - TEMP2 * 0.5d0*(KAPPA_T(I,K,J,N)+KAPPA_T(I,K,J-1,N))*GMAT_22(K,J,1)
              MATDY(I,J) = MATDY(I,J)                                                  &
           + TEMP2 * 0.5d0*(KAPPA_T(I,K,J,N) + KAPPA_T(I,K,J+1,N))*GMAT_22(K,J+1,1)  &
           + TEMP2 * 0.5d0*(KAPPA_T(I,K,J,N) + KAPPA_T(I,K,J-1,N))*GMAT_22(K,J,1)
              MATUY(I,J) = MATUY(I,J)                                                  &
           - TEMP2 * 0.5d0*(KAPPA_T(I,K,J,N) + KAPPA_T(I,K,J+1,N))*GMAT_22(K,J+1,1)
          END DO
        END DO

        END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! C Set the boundary conditions for TH
          CALL APPLY_BC_TH_LOWER(N,K)
          CALL APPLY_BC_TH_UPPER(N,K)
! C Now, solve the tridiagonal system for TH(:,k,:)
          CALL THOMAS_REAL(MATLY,MATDY,MATUY,VECY,NY+1,NXP)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



        DO J=0,NY+1
          DO I=0,NXP
            THX(I,K,J,N)=VECY(I,J)
          END DO
        END DO
! End do k
      END DO
!$OMP END PARALLEL DO

      
      
! C Compute the horizontal viscous term in physical space and add to new RHS
! C already store in S1 
! C Important we have to multiply 1/2 with S1 before adding to RHS due
! C half time splitting during AD1

      DO J=JSTART_TH(N),JEND_TH(N)
        DO K=ZSTART_TH(N),ZEND_TH(N)
          DO I=0,NXP
            RTHX(I,K,J,N)=0.5*S1X(I,K,J) + INT_JACOB(K,J)*THX(I,K,J,N)  &
              + (TEMP1/PR(N))                                        &
             *(  GMAT_22(K,J+1,1)*(THX(I,K,J+1,N) - THX(I,K,J,N))      &
              - GMAT_22(K,J,1)*(THX(I,K,J,N)   - THX(I,K,J-1,N))  )

          END DO
        END DO
      END DO


      IF (LES.AND.((.NOT.CREATE_NEW_FLOW).OR.(TIME_STEP.GT.1))) THEN
       DO J=JSTART_TH(N),JEND_TH(N)
        DO K=ZSTART_TH(N),ZEND_TH(N)
          DO I=0,NXP
            RTHX(I,K,J,N)=RTHX(I,K,J,N)+TEMP2*(0.5d0*(KAPPA_T(I,K,J+1,N)+KAPPA_T(I,K,J,N)) &
                    *(GMAT_22(K,J+1,1)*(THX(I,K,J+1,N) - THX(I,K,J,N)))                    &
                    -0.5d0*(KAPPA_T(I,K,J-1,N)+KAPPA_T(I,K,J,N))                         &
                    *(GMAT_22(K,J,1)*(THX(I,K,J,N)   - THX(I,K,J-1,N))))
         END DO
        END DO
       END DO
      ENDIF


! c      IF ( IC_TYPE == 7 ) THEN
! c       DO K=ZSTART_TH(N),ZEND_TH(N)
! c	   DO I=0,NXM
! c             RTHX(I,K,NY+1,N)=RTHX(I,K,NY,N)
! c       	   END DO
! c       END DO  
! c      ENDIF
! Initialize the matrix used to store implicit coefficients
      DO K=0,NZ+1
        DO I=0,NXP
          MATLX(I,K)=0.
          MATDX(I,K)=1.
          MATUX(I,K)=0.
          VECX(I,K)=0.
        END DO
      END DO 

! C Build the implicit system of equations for TH 
!$OMP PARALLEL DO PRIVATE(K,J,I)
      DO J=JSTART_TH(N),JEND_TH(N)
       IF(WAVE_ABS) THEN
        DO K=ZSTART_TH(N),ZEND_TH(N)
          DO I=0,NXP
            MATLX(I,K)=-(TEMP1/PR(N))*GMAT_11(K,J,2)*  &
           (1.d0 + 5.d0*SPONGE_SIGMA_OUT(K,J))           
            MATDX(I,K)=INT_JACOB(K,J)-(TEMP1/PR(N))*    &
           (1.d0 + 5.d0*SPONGE_SIGMA_OUT(K,J))*            &
                       (-GMAT_11(K+1,J,2) - GMAT_11(K,J,2))
            MATUX(I,K)=-(TEMP1/PR(N))*GMAT_11(K+1,J,2)*  &
           (1.d0 + 5.d0*SPONGE_SIGMA_OUT(K,J))
            VECX(I,K)=RTHX(I,K,J,N)
          END DO
        END DO
       ELSE
        DO K=ZSTART_TH(N),ZEND_TH(N)
          DO I=0,NXP
            MATLX(I,K)=-(TEMP1/PR(N))*GMAT_11(K,J,2)
            MATDX(I,K)=INT_JACOB(K,J)-(TEMP1/PR(N))*   &
      	                (-GMAT_11(K+1,J,2) - GMAT_11(K,J,2)) 
            MATUX(I,K)=-(TEMP1/PR(N))*GMAT_11(K+1,J,2)
            VECX(I,K)=RTHX(I,K,J,N)
          END DO
        END DO
       ENDIF 

! IF using a subgrid model (LES) then add the eddy viscosity part implicitly

        IF (LES.AND.((.NOT.CREATE_NEW_FLOW).OR.(TIME_STEP.GT.1))) THEN
          DO K=ZSTART_TH(N),ZEND_TH(N)
            DO I=0,NXP
              MATLX(I,K) = MATLX(I,K) &
           - TEMP2 * 0.5d0*(KAPPA_T(I,K,J,N)+KAPPA_T(I,K-1,J,N))*GMAT_11(K,J,2)
              MATDX(I,K) = MATDX(I,K)                                         &
           + TEMP2 * 0.5d0*(KAPPA_T(I,K+1,J,N)+KAPPA_T(I,K,J,N))*GMAT_11(K+1,J,2)   &
           + TEMP2 * 0.5d0*(KAPPA_T(I,K,J,N)+KAPPA_T(I,K-1,J,N))*GMAT_11(K,J,2)
              MATUX(I,K) = MATUX(I,K)                                         &
           - TEMP2 * 0.5d0*(KAPPA_T(I,K+1,J,N)+KAPPA_T(I,K,J,N))*GMAT_11(K+1,J,2)
            END DO
          END DO
        END IF


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! C Set the boundary conditions for U1
          CALL APPLY_BC_TH_LEFT(N,J)
          CALL APPLY_BC_TH_RIGHT(N,J)
! C Now, solve the tridiagonal system for TH(:,k,:)
        IF ( TH_BC_ZMAX(N)  .EQ. 5 ) THEN
         D=1.0
         CALL THOMAS_REAL_SP(MATLX,MATDX,MATUX,VECX,D,NZ,NZ+1,NXP)
        ELSE   
         CALL THOMAS_REAL(MATLX,MATDX,MATUX,VECX,NZ+1,NXP)
        ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

        DO K=0,NZ+1
          DO I=0,NXP
            THX(I,K,J,N)=VECX(I,K)
          END DO
        END DO

        
! End do J
      END DO
!$OMP END PARALLEL DO

! C NEAD TO UPDATE THE CORNERS


      IF ( TH_BC_ZMAX(N)  .EQ. 1 ) THEN
	   DO I=0,NXP
             THX(I,NZ+1,NY+1,N)=THX(I,NZ,NY+1,N)
             THX(I,NZ+1,0,N)=THX(I,NZ,0,N)
       	   END DO
      ENDIF
 
      IF ( TH_BC_ZMIN(N)  .EQ. 1 ) THEN
	   DO I=0,NXP
             THX(I,0,NY+1,N)=THX(I,1,NY+1,N)
             THX(I,0,0,N)=THX(I,1,0,N)
       	   END DO  
      ENDIF
     
 
      END DO

! c      read(6,*)
! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! CCCCCCCCC                   Now same for velocity          CCCCCCCCCCCCCCCCCC
! 
! 
! C ADI STEP FOR  U1 
! C Before going to  ADI1 we have to store RHS after subtracting previous U1_n 
! C at previous time step to bulid new RHS for  ADI2 based on U1_n+1/2 at intermidiate 



      DO J=JSTART,JEND
        DO K=ZSTART,ZEND
          DO I=0,NXP
            S1X(I,K,J)= R1X(I,K,J) - INT_JACOB(K,J)*U1X(I,K,J)
          END DO
        END DO
      END DO


! C Compute the horizontal viscous term in physical space and add to RHS
! C  at  ADI1 step .... Important in this step we have to add U1_n
! C at previous time step so that after multiplying a factor 1/2 into the 
! C RHS side during ADI1 operation RHS will be taken care U1_n not 1/2*U1_n
 

      DO J=JSTART,JEND
        DO K=ZSTART,ZEND
          DO I=0,NXP
            R1X(I,K,J)=0.5*S1X(I,K,J) + INT_JACOB(K,J)*U1X(I,K,J)   &
            + TEMP1*(GMAT_11(K+1,J,2)*(U1X(I,K+1,J) - U1X(I,K,J)) &
                    -GMAT_11(K,J,2)*(U1X(I,K,J)   - U1X(I,K-1,J)))
          END DO
        END DO
      END DO

      IF (LES.AND.((.NOT.CREATE_NEW_FLOW).OR.(TIME_STEP.GT.1))) THEN
       DO J=JSTART,JEND
        DO K=ZSTART,ZEND
          DO I=0,NXP
            R1X(I,K,J)=R1X(I,K,J)+TEMP2*(0.5d0*(NU_T(I,K+1,J)+NU_T(I,K,J)) &
                    *(GMAT_11(K+1,J,2)*(U1X(I,K+1,J) - U1X(I,K,J)))        &
                    -0.5d0*(NU_T(I,K,J)+NU_T(I,K-1,J))                  &
                    *(GMAT_11(K,J,2)*(U1X(I,K,J)   - U1X(I,K-1,J))))
         END DO
        END DO
       END DO
      ENDIF

       
!       
! C Solve for U1
! C Note, here the matrix will be indexed from 1...NY+1 corresponding to U1(0:NY)

! Initialize the matrix used to store implicit coefficients
      DO J=0,NY+1
        DO I=0,NXP
          MATLY(I,J)=0.
          MATDY(I,J)=1.
          MATUY(I,J)=0.
          VECY(I,J)=0.
        END DO
      END DO 
! 
! C Build the implicit system of equations for U1 
!$OMP PARALLEL DO PRIVATE(K,J,I)
      DO K=ZSTART,ZEND

        DO J=JSTART,JEND
          DO I=0,NXP
            MATLY(I,J)=-TEMP1*GMAT_22(K,J,1)
            MATDY(I,J)= INT_JACOB(K,J)-TEMP1*  &
              (-GMAT_22(K,J+1,1) - GMAT_22(K,J,1)) 
            MATUY(I,J)=-TEMP1*GMAT_22(K,J+1,1)
            VECY(I,J)=R1X(I,K,J)
          END DO
        END DO

! IF using a subgrid model (LES) then add the eddy viscosity part implicitly
       

         IF (LES.AND.((.NOT.CREATE_NEW_FLOW).OR.(TIME_STEP.GT.1))) THEN
          DO J=JSTART,JEND
          DO I=0,NXP
              MATLY(I,J) = MATLY(I,J) &
           - TEMP2 * 0.5d0*(NU_T(I,K,J)+NU_T(I,K,J-1))*GMAT_22(K,J,1)
              MATDY(I,J) = MATDY(I,J)                &
           + TEMP2 * 0.5d0*(NU_T(I,K,J)+NU_T(I,K,J+1))*GMAT_22(K,J+1,1) &
           + TEMP2 * 0.5d0*(NU_T(I,K,J)+NU_T(I,K,J-1))*GMAT_22(K,J,1)
              MATUY(I,J) = MATUY(I,J)     &
           - TEMP2 * 0.5d0*(NU_T(I,K,J)+NU_T(I,K,J+1))*GMAT_22(K,J+1,1)
          END DO
        END DO
        END IF



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! C Else, we are running in serial mode
! C Set the boundary conditions for U1
          CALL APPLY_BC_1_LOWER(K)
          CALL APPLY_BC_1_UPPER(K)
! C Now, solve the tridiagonal system for U1(:,k,:)
          CALL THOMAS_REAL(MATLY,MATDY,MATUY,VECY,NY+1,NXP)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


        DO J=0,NY+1
          DO I=0,NXP
            U1X(I,K,J)=VECY(I,J)
          END DO
        END DO
! End do k
      END DO
!$OMP END PARALLEL DO

         	
      
      
      
! C Compute the horizontal viscous term in physical space and add to new RHS
! C already store in S1 
! C Important we have to multiply 1/2 with S1 before adding to RHS due
! C half time splitting during AD1

      DO J=JSTART,JEND
        DO K=ZSTART,ZEND
	   DO I=0,NXP
            R1X(I,K,J)=0.5*S1X(I,K,J) + INT_JACOB(K,J)*U1X(I,K,J)  &
           + TEMP1*(GMAT_22(K,J+1,1)*(U1X(I,K,J+1) - U1X(I,K,J)) &
                  - GMAT_22(K,J,1)*(U1X(I,K,J)   - U1X(I,K,J-1)))

          END DO
        END DO
      END DO

      IF (LES.AND.((.NOT.CREATE_NEW_FLOW).OR.(TIME_STEP.GT.1))) THEN
       DO J=JSTART,JEND
        DO K=ZSTART,ZEND
          DO I=0,NXP
            R1X(I,K,J)=R1X(I,K,J)+TEMP2*(0.5d0*(NU_T(I,K,J+1)+NU_T(I,K,J)) &
                    *(GMAT_22(K,J+1,1)*(U1X(I,K,J+1) - U1X(I,K,J)))        &
                    -0.5d0*(NU_T(I,K,J-1)+NU_T(I,K,J))                   &
                    *(GMAT_22(K,J,1)*(U1X(I,K,J)   - U1X(I,K,J-1))))
         END DO
        END DO
       END DO
      ENDIF


      IF ( IC_TYPE == 7 ) THEN
       DO K=ZSTART,ZEND
	   DO I=0,NXP
             R1X(I,K,NY+1)=R1X(I,K,NY)
       	   END DO
       END DO  
      ENDIF
      
! Initialize the matrix used to store implicit coefficients
      DO K=0,NZ+1
        DO I=0,NXP
          MATLX(I,K)=0.
          MATDX(I,K)=1.
          MATUX(I,K)=0.
          VECX(I,K)=0.
        END DO
      END DO 

! C Build the implicit system of equations for TH 
!$OMP PARALLEL DO PRIVATE(K,J,I)
      DO J=JSTART,JEND
        DO K=ZSTART,ZEND
          DO I=0,NXP
            MATLX(I,K)=-TEMP1*GMAT_11(K,J,2)
            MATDX(I,K)=INT_JACOB(K,J)  - TEMP1* &
      	                (-GMAT_11(K+1,J,2) - GMAT_11(K,J,2)) 
            MATUX(I,K)=-TEMP1*GMAT_11(K+1,J,2)
            VECX(I,K)=R1X(I,K,J)
          END DO
        END DO
! IF using a subgrid model (LES) then add the eddy viscosity part implicitly
         
        IF (LES.AND.((.NOT.CREATE_NEW_FLOW).OR.(TIME_STEP.GT.1))) THEN
          DO K=ZSTART,ZEND
            DO I=0,NXP
              MATLX(I,K) = MATLX(I,K) &
           - TEMP2 * 0.5d0*(NU_T(I,K,J)+NU_T(I,K-1,J))*GMAT_11(K,J,2)
              MATDX(I,K) = MATDX(I,K)    &
           + TEMP2 * 0.5d0*(NU_T(I,K+1,J)+NU_T(I,K,J))*GMAT_11(K+1,J,2) &
           + TEMP2 * 0.5d0*(NU_T(I,K,J)+NU_T(I,K-1,J))*GMAT_11(K,J,2)
              MATUX(I,K) = MATUX(I,K)       &
           - TEMP2 * 0.5d0*(NU_T(I,K+1,J)+NU_T(I,K,J))*GMAT_11(K+1,J,2)
            END DO
          END DO
        END IF


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! C Set the boundary conditions for U1
          CALL APPLY_BC_1_LEFT(J)
          CALL APPLY_BC_1_RIGHT(J)
! C Now, solve the tridiagonal system for TH(:,k,:)
        IF ( U_BC_ZMAX  .EQ. 5 ) THEN
         D=GMAT_11(NZ,J,2)
         CALL THOMAS_REAL_SP(MATLX,MATDX,MATUX,VECX,D,NZ,NZ+1,NXP)
        ELSE   
         CALL THOMAS_REAL(MATLX,MATDX,MATUX,VECX,NZ+1,NXP)
        ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


        DO K=0,NZ+1
          DO I=0,NXP
            U1X(I,K,J) = VECX(I,K)
          END DO
        END DO
! End do J
      END DO
!$OMP END PARALLEL DO
   
! c      read(6,*) 
      
 
      
! C ADI STEP FOR  U2 
! C Before going to  ADI1 we have to store RHS after subtracting previous U2_n 
! C at previous time step to bulid new RHS for  ADI2 based on U2_n+1/2 at intermidiate 



      DO J=JSTART,JEND
        DO K=ZSTART,ZEND
          DO I=0,NXP
            S1X(I,K,J)= R2X(I,K,J) - INT_JACOB(K,J)*U2X(I,K,J)
          END DO
        END DO
      END DO


! C Compute the horizontal viscous term in physical space and add to RHS
! C  at  ADI1 step .... Important in this step we have to add U2_n
! C at previous time step so that after multiplying a factor 1/2 into the 
! C RHS side during ADI1 operation RHS will be taken care U2_n not 1/2*U2_n
 
      IF (WAVE_ABS) THEN
       DO J=JSTART,JEND
        DO K=ZSTART,ZEND
         DO I=0,NXP
! C    wave absorb layer
          R2X(I,K,J)=0.5*S1X(I,K,J) + INT_JACOB(K,J)*U2X(I,K,J)  &
            + TEMP1*(1.d0 + 5.d0*SPONGE_SIGMA_OUT(K,J))*        &
                    (GMAT_11(K+1,J,2)*(U2X(I,K+1,J) - U2X(I,K,J))  &
                   -GMAT_11(K,J,2)*(U2X(I,K,J)   - U2X(I,K-1,J)))
          END DO
        END DO
       END DO
      ELSE
      DO J=JSTART,JEND
        DO K=ZSTART,ZEND
          DO I=0,NXP
            R2X(I,K,J)=0.5*S1X(I,K,J) + INT_JACOB(K,J)*U2X(I,K,J)  &
            + TEMP1*(GMAT_11(K+1,J,2)*(U2X(I,K+1,J) - U2X(I,K,J)) &
                    -GMAT_11(K,J,2)*(U2X(I,K,J)   - U2X(I,K-1,J)))
          END DO
        END DO
      END DO
      ENDIF

     IF (LES.AND.((.NOT.CREATE_NEW_FLOW).OR.(TIME_STEP.GT.1))) THEN
       DO J=JSTART,JEND
        DO K=ZSTART,ZEND
          DO I=0,NXP
            R2X(I,K,J)=R2X(I,K,J)+TEMP2*(0.5d0*(NU_T(I,K+1,J)+NU_T(I,K,J)) &
                    *(GMAT_11_y(K+1,J,2)*(U2X(I,K+1,J) - U2X(I,K,J)))        &
                    -0.5d0*(NU_T(I,K,J)+NU_T(I,K-1,J))                  &
                    *(GMAT_11_y(K,J,2)*(U2X(I,K,J)   - U2X(I,K-1,J))))
         END DO
        END DO
       END DO
      ENDIF

! C Solve for U2
! C Note, here the matrix will be indexed from 1...NY+1 corresponding to U1(0:NY)

! Initialize the matrix used to store implicit coefficients
      DO J=0,NY+1
        DO I=0,NXP
          MATLY(I,J)=0.
          MATDY(I,J)=1.
          MATUY(I,J)=0.
          VECY(I,J)=0.
        END DO
      END DO 

! C Build the implicit system of equations for U1 
!$OMP PARALLEL DO PRIVATE(K,J,I)
      DO K=ZSTART,ZEND

        DO J=JSTART,JEND
          DO I=0,NXP
            MATLY(I,J)=-TEMP1*GMAT_22(K,J,1)  
            MATDY(I,J)= INT_JACOB(K,J)-TEMP1*  &
              (-GMAT_22(K,J+1,1) - GMAT_22(K,J,1)) 
            MATUY(I,J)=-TEMP1*GMAT_22(K,J+1,1)
            VECY(I,J)=R2X(I,K,J)
          END DO
        END DO

! IF using a subgrid model (LES) then add the eddy viscosity part implicitly

        IF (LES.AND.((.NOT.CREATE_NEW_FLOW).OR.(TIME_STEP.GT.1))) THEN
          DO J=JSTART,JEND
          DO I=0,NXP
              MATLY(I,J) = MATLY(I,J) &
           - TEMP2 * 0.5d0*(NU_T(I,K,J)+NU_T(I,K,J-1))*GMAT_22_y(K,J,1)
              MATDY(I,J) = MATDY(I,J)   &
           + TEMP2 * 0.5d0*(NU_T(I,K,J)+NU_T(I,K,J+1))*GMAT_22_y(K,J+1,1) &
           + TEMP2 * 0.5d0*(NU_T(I,K,J)+NU_T(I,K,J-1))*GMAT_22_y(K,J,1)
              MATUY(I,J) = MATUY(I,J)    &
           - TEMP2 * 0.5d0*(NU_T(I,K,J)+NU_T(I,K,J+1))*GMAT_22_y(K,J+1,1)
          END DO
        END DO
        END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! C Set the boundary conditions for U1
          CALL APPLY_BC_2_LOWER(K)
          CALL APPLY_BC_2_UPPER(K)
! C Now, solve the tridiagonal system for U1(:,k,:)
          CALL THOMAS_REAL(MATLY,MATDY,MATUY,VECY,NY+1,NXP)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


        DO J=0,NY+1
          DO I=0,NXP
            U2X(I,K,J)=VECY(I,J)
          END DO
        END DO
! End do k
      END DO
!$OMP END PARALLEL DO

      
      
!       
! C Compute the horizontal viscous term in physical space and add to new RHS
! C already store in S1 
! C Important we have to multiply 1/2 with S1 before adding to RHS due
! C half time splitting during AD1

      DO J=JSTART,JEND
        DO K=ZSTART,ZEND
	   DO I=0,NXP
            R2X(I,K,J)=0.5*S1X(I,K,J) + INT_JACOB(K,J)*U2X(I,K,J)  &
           + TEMP1*(GMAT_22(K,J+1,1)*(U2X(I,K,J+1) - U2X(I,K,J)) &
                  - GMAT_22(K,J,1)*(U2X(I,K,J) - U2X(I,K,J-1)))

          END DO
        END DO
      END DO

      IF (LES.AND.((.NOT.CREATE_NEW_FLOW).OR.(TIME_STEP.GT.1))) THEN
       DO J=JSTART,JEND
        DO K=ZSTART,ZEND
          DO I=0,NXP
            R2X(I,K,J)=R2X(I,K,J)+TEMP2*(0.5d0*(NU_T(I,K,J+1)+NU_T(I,K,J)) &
                    *(GMAT_22_y(K,J+1,1)*(U2X(I,K,J+1) - U2X(I,K,J)))        &
                    -0.5d0*(NU_T(I,K,J-1)+NU_T(I,K,J))                   &
                    *(GMAT_22_y(K,J,1)*(U2X(I,K,J)   - U2X(I,K,J-1))))
         END DO
        END DO
       END DO
      ENDIF


      IF ( IC_TYPE == 7 ) THEN
       DO K=ZSTART,ZEND
	   DO I=0,NXP
             R2X(I,K,NY+1)=R2X(I,K,NY)
       	   END DO
       END DO  
      ENDIF 

! Initialize the matrix used to store implicit coefficients
      DO K=0,NZ+1
        DO I=0,NXP
          MATLX(I,K)=0.
          MATDX(I,K)=1.
          MATUX(I,K)=0.
          VECX(I,K)=0.
        END DO
      END DO 

! C Build the implicit system of equations for U2 
!$OMP PARALLEL DO PRIVATE(K,J,I)
      DO J=JSTART,JEND
       IF (WAVE_ABS) THEN
        DO K=ZSTART,ZEND
          DO I=0,NXP
! C wave absorver layer
            MATLX(I,K)=-TEMP1*(1.d0 + 5.d0*SPONGE_SIGMA_OUT(K,J)) &
                        *GMAT_11(K,J,2)
            MATDX(I,K)=INT_JACOB(K,J)  - TEMP1*          &
                       (1.d0 + 5.d0*SPONGE_SIGMA_OUT(K,J))* &
                              (-GMAT_11(K+1,J,2) - GMAT_11(K,J,2))
            MATUX(I,K)=-TEMP1*(1.d0 + 5.d0*SPONGE_SIGMA_OUT(K,J))*  &
               GMAT_11(K+1,J,2)
            VECX(I,K)=R2X(I,K,J)

          END DO
        END DO
       ELSE 
        DO K=ZSTART,ZEND
          DO I=0,NXP
            MATLX(I,K)=-TEMP1*GMAT_11(K,J,2)
            MATDX(I,K)=INT_JACOB(K,J)  - TEMP1*  &
      	                (-GMAT_11(K+1,J,2) - GMAT_11(K,J,2)) 
            MATUX(I,K)=-TEMP1*GMAT_11(K+1,J,2)
            VECX(I,K)=R2X(I,K,J)
          END DO
        END DO
       ENDIF
! IF using a subgrid model (LES) then add the eddy viscosity part implicitly
       
       IF (LES.AND.((.NOT.CREATE_NEW_FLOW).OR.(TIME_STEP.GT.1))) THEN
          DO K=ZSTART,ZEND
            DO I=0,NXP
              MATLX(I,K) = MATLX(I,K) &
           - TEMP2 * 0.5d0*(NU_T(I,K,J)+NU_T(I,K-1,J))*GMAT_11_y(K,J,2)
              MATDX(I,K) = MATDX(I,K) &
           + TEMP2 * 0.5d0*(NU_T(I,K+1,J)+NU_T(I,K,J))*GMAT_11_y(K+1,J,2) &
           + TEMP2 * 0.5d0*(NU_T(I,K,J)+NU_T(I,K-1,J))*GMAT_11_y(K,J,2)
              MATUX(I,K) = MATUX(I,K)      &
           - TEMP2 * 0.5d0*(NU_T(I,K+1,J)+NU_T(I,K,J))*GMAT_11_y(K+1,J,2)
            END DO
          END DO
        END IF 


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! C Set the boundary conditions for U1
          CALL APPLY_BC_2_LEFT(J)
          CALL APPLY_BC_2_RIGHT(J)
! C Now, solve the tridiagonal system for TH(:,k,:)
        IF ( V_BC_ZMAX  .EQ. 5 ) THEN
         D=1.0d0
        CALL THOMAS_REAL_SP(MATLX,MATDX,MATUX,VECX,D,NZ,NZ+1,NXP)
        ELSE   
         CALL THOMAS_REAL(MATLX,MATDX,MATUX,VECX,NZ+1,NXP)
        ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        DO K=0,NZ+1
          DO I=0,NXP
            U2X(I,K,J) = VECX(I,K)
          END DO
        END DO
! End do J
      END DO
!$OMP END PARALLEL DO
     
      
! C ADI STEP FOR  U3 
! C Before going to  ADI1 we have to store RHS after subtracting previous U3_n 
! C at previous time step to bulid new RHS for  ADI2 based on U3_n+1/2 at intermidiate 

      

      DO J=JSTART,JEND
        DO K=ZSTART,ZEND
          DO I=0,NXP
            S1X(I,K,J)= R3X(I,K,J) - INT_JACOB(K,J)*U3X(I,K,J)
          END DO
        END DO
      END DO


! C Compute the horizontal viscous term in physical space and add to RHS
! C  at  ADI1 step .... Important in this step we have to add U3_n
! C at previous time step so that after multiplying a factor 1/2 into the 
! C RHS side during ADI1 operation RHS will be taken care U3_n not 1/2*U3_n
 
      IF (WAVE_ABS) THEN
        DO J=JSTART,JEND
         DO K=ZSTART,ZEND
          DO I=0,NXP
           R3X(I,K,J)=0.5*S1X(I,K,J) + INT_JACOB(K,J)*U3X(I,K,J) &
            + TEMP1*(1.d0 + 5.d0*SPONGE_SIGMA_OUT(K,J))*         &
            (GMAT_11(K+1,J,2)*(U3X(I,K+1,J) - U3X(I,K,J))        &
             -GMAT_11(K,J,2)*(U3X(I,K,J)   - U3X(I,K-1,J)))
          END DO
        END DO
      END DO
      ELSE
      DO J=JSTART,JEND
        DO K=ZSTART,ZEND
          DO I=0,NXP
            R3X(I,K,J)=0.5*S1X(I,K,J) + INT_JACOB(K,J)*U3X(I,K,J)    &
            + TEMP1*(GMAT_11(K+1,J,2)*(U3X(I,K+1,J) - U3X(I,K,J))  &
                    -GMAT_11(K,J,2)*(U3X(I,K,J)   - U3X(I,K-1,J)))
          END DO
        END DO
      END DO
      ENDIF

      
      IF (LES.AND.((.NOT.CREATE_NEW_FLOW).OR.(TIME_STEP.GT.1))) THEN
       DO J=JSTART,JEND
        DO K=ZSTART,ZEND
          DO I=0,NXP
            R3X(I,K,J)=R3X(I,K,J)+TEMP2*(0.5d0*(NU_T(I,K+1,J)+NU_T(I,K,J)) &
                    *(GMAT_11_z(K+1,J,2)*(U3X(I,K+1,J) - U3X(I,K,J)))        &
                    -0.5d0*(NU_T(I,K,J)+NU_T(I,K-1,J))                  &
                    *(GMAT_11_z(K,J,2)*(U3X(I,K,J)   - U3X(I,K-1,J))))
         END DO
        END DO
       END DO
      ENDIF
 
! C Solve for U3
! C Note, here the matrix will be indexed from 1...NY+1 corresponding to U3(0:NY)

! Initialize the matrix used to store implicit coefficients
      DO J=0,NY+1
        DO I=0,NXP
          MATLY(I,J)=0.
          MATDY(I,J)=1.
          MATUY(I,J)=0.
          VECY(I,J)=0.
        END DO
      END DO 

! C Build the implicit system of equations for U3 
!$OMP PARALLEL DO PRIVATE(K,J,I)
      DO K=ZSTART,ZEND

        DO J=JSTART,JEND
          DO I=0,NXP
            MATLY(I,J)=-TEMP1*GMAT_22(K,J,1)
            MATDY(I,J)= INT_JACOB(K,J)-TEMP1*  &
              (-GMAT_22(K,J+1,1) - GMAT_22(K,J,1)) 
            MATUY(I,J)=-TEMP1*GMAT_22(K,J+1,1)
            VECY(I,J)=R3X(I,K,J)
          END DO
        END DO

! IF using a subgrid model (LES) then add the eddy viscosity part implicitly

       IF (LES.AND.((.NOT.CREATE_NEW_FLOW).OR.(TIME_STEP.GT.1))) THEN
         DO J=JSTART,JEND
           DO I=0,NXP
              MATLY(I,J) = MATLY(I,J) &
           - TEMP2 * 0.5d0*(NU_T(I,K,J)+NU_T(I,K,J-1))*GMAT_22_z(K,J,1)
              MATDY(I,J) = MATDY(I,J) &
           + TEMP2 * 0.5d0*(NU_T(I,K,J)+NU_T(I,K,J+1))*GMAT_22_z(K,J+1,1) &
           + TEMP2 * 0.5d0*(NU_T(I,K,J)+NU_T(I,K,J-1))*GMAT_22_z(K,J,1)
              MATUY(I,J) = MATUY(I,J) &
           - TEMP2 * 0.5d0*(NU_T(I,K,J)+NU_T(I,K,J+1))*GMAT_22_z(K,J+1,1)
          END DO
         END DO
        END IF        


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! C Set the boundary conditions for U1
          CALL APPLY_BC_3_LOWER(K)
          CALL APPLY_BC_3_UPPER(K)
! C Now, solve the tridiagonal system for U1(:,k,:)
          CALL THOMAS_REAL(MATLY,MATDY,MATUY,VECY,NY+1,NXP)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


        DO J=0,NY+1
          DO I=0,NXP
            U3X(I,K,J)=VECY(I,J)
          END DO
        END DO
! End do k
      END DO
!$OMP END PARALLEL DO

      
! C     requred for convective bc
!      C_avg = 0.0
!      DO J=JSTART,JEND
!	DO I=0,NXP
!	  C_avg = C_avg + U3X(I,NZ+1,J)
!        ENDDO
!      ENDDO	
        	        
!      C_avg = C_avg/real(NX*(JEND-JSTART+1))
      
!      DO J=JSTART,JEND
!        C_int(j) = 0.d0
!	C_int_le(j) = 0.d0
!	DO I=0,NXP
!	  C_int(J) = C_int(J) + U3X(I,NZ,J)
!	  C_int_le(J) = C_int_le(J) + U3X(I,1,J)
!        ENDDO
!	C_int(j) = C_int(j)/real(NX)
!	C_int_le(j) = C_int_le(j)/real(NX)
!      ENDDO
     	
      
! c      DO J=JSTART,JEND  	        
! c       C_int(j) = C_int(j)/real(NX)
! c      ENDDO
!            
! C Compute the horizontal viscous term in physical space and add to new RHS
! C already store in S1 
! C Important we have to multiply 1/2 with S1 before adding to RHS due
! C half time splitting during AD1

      DO J=JSTART,JEND
        DO K=ZSTART,ZEND
	   DO I=0,NXP
            R3X(I,K,J)=0.5*S1X(I,K,J) + INT_JACOB(K,J)*U3X(I,K,J)  &
           + TEMP1*(GMAT_22(K,J+1,1)*(U3X(I,K,J+1) - U3X(I,K,J)) &
                  - GMAT_22(K,J,1)*(U3X(I,K,J) - U3X(I,K,J-1)))

          END DO
        END DO
      END DO

      IF (LES.AND.((.NOT.CREATE_NEW_FLOW).OR.(TIME_STEP.GT.1))) THEN
       DO J=JSTART,JEND
        DO K=ZSTART,ZEND
          DO I=0,NXP
            R3X(I,K,J)=R3X(I,K,J)+TEMP2*(0.5d0*(NU_T(I,K,J+1)+NU_T(I,K,J)) &
                    *(GMAT_22_z(K,J+1,1)*(U3X(I,K,J+1) - U3X(I,K,J)))        &
                    -0.5d0*(NU_T(I,K,J-1)+NU_T(I,K,J))                   &
                    *(GMAT_22_z(K,J,1)*(U3X(I,K,J)   - U3X(I,K,J-1))))
         END DO
        END DO
       END DO
      ENDIF 


      IF ( IC_TYPE == 7) THEN
        DO K=ZSTART,ZEND
	   DO I=0,NXP
	     R3X(I,K,NY+1)=R3X(I,K,NY)
	   ENDDO
	ENDDO   
      ENDIF

! Initialize the matrix used to store implicit coefficients
      DO K=0,NZ+1
        DO I=0,NXP
          MATLX(I,K)=0.
          MATDX(I,K)=1.
          MATUX(I,K)=0.
          VECX(I,K)=0.
        END DO
      END DO 

! C Build the implicit system of equations for U3
!$OMP PARALLEL DO PRIVATE(K,J,I)
      DO J=JSTART,JEND+1
       IF (WAVE_ABS) THEN
        DO K=ZSTART,ZEND
          DO I=0,NXP
            MATLX(I,K)=-TEMP1*GMAT_11(K,J,2)*   &
                       (1.d0 + 5.d0*SPONGE_SIGMA_OUT(K,J))
            MATDX(I,K)=INT_JACOB(K,J)  - TEMP1*   &
                      (1.d0 + 5.d0*SPONGE_SIGMA_OUT(K,J))*  &
                       (-GMAT_11(K+1,J,2) - GMAT_11(K,J,2))
            MATUX(I,K)=-TEMP1*GMAT_11(K+1,J,2)*    &
                       (1.d0 + 5.d0*SPONGE_SIGMA_OUT(K,J))
            VECX(I,K)=R3X(I,K,J)
          END DO
        END DO
       ELSE
        DO K=ZSTART,ZEND
          DO I=0,NXP
            MATLX(I,K)=-TEMP1*GMAT_11(K,J,2) 
            MATDX(I,K)=INT_JACOB(K,J)  - TEMP1*  &
      	                (-GMAT_11(K+1,J,2) - GMAT_11(K,J,2)) 
            MATUX(I,K)=-TEMP1*GMAT_11(K+1,J,2)
            VECX(I,K)=R3X(I,K,J)
          END DO
        END DO
       ENDIF
! ! IF using a subgrid model (LES) then add the eddy viscosity part implicitly

        IF (LES.AND.((.NOT.CREATE_NEW_FLOW).OR.(TIME_STEP.GT.1))) THEN
          DO K=ZSTART,ZEND
            DO I=0,NXP
              MATLX(I,K) = MATLX(I,K) &
           - TEMP2 * 0.5d0*(NU_T(I,K,J)+NU_T(I,K-1,J))*GMAT_11_z(K,J,2)
              MATDX(I,K) = MATDX(I,K) &
           + TEMP2 * 0.5d0*(NU_T(I,K+1,J)+NU_T(I,K,J))*GMAT_11_z(K+1,J,2) &
           + TEMP2 * 0.5d0*(NU_T(I,K,J)+NU_T(I,K-1,J))*GMAT_11_z(K,J,2)
              MATUX(I,K) = MATUX(I,K)      &
           - TEMP2 * 0.5d0*(NU_T(I,K+1,J)+NU_T(I,K,J))*GMAT_11_z(K+1,J,2)
            END DO
          END DO
        END IF


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! C Set the boundary conditions for U1
          CALL APPLY_BC_3_LEFT(J)
          CALL APPLY_BC_3_RIGHT(J)
! C Now, solve the tridiagonal system for TH(:,k,:)
        IF ( W_BC_ZMAX  .EQ. 5 ) THEN
         D=1.
         CALL THOMAS_REAL_SP(MATLX,MATDX,MATUX,VECX,D,NZ,NZ+1,NXP)
        ELSE   
         CALL THOMAS_REAL(MATLX,MATDX,MATUX,VECX,NZ+1,NXP)
        ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        DO K=0,NZ+1
          DO I=0,NXP
            U3X(I,K,J) = VECX(I,K)
          END DO
        END DO
! End do J
      END DO
!$OMP END PARALLEL DO     

     
! C If Variable timestepping and done with one full R-K step, update
! C DELTA_T based on the specified CFL number
! C This is not parallelized and should be used only in the serial
! C version to ensure that each process uses the same timestep
      IF ((VARIABLE_DT).and.(RK_STEP.eq.3) &
             .and.(MOD(TIME_STEP,UPDATE_DT).EQ.0)) THEN
         IF (USE_MPI)THEN
          CALL COURANT_CURV_MPI
         ELSE
          CALL COURANT_CURV
         ENDIF 
      END IF 

      CALL REAL_FOURIER_TRANS_U1 (.true.)
      CALL REAL_FOURIER_TRANS_U2 (.true.) 
      CALL REAL_FOURIER_TRANS_U3 (.true.) 

      call allocation_R1 (.false.)
      call allocation_R2 (.false.)
      call allocation_R3 (.false.)  

!C      CALL FFT_X_TO_FOURIER(U1,CU1,0,NY+1,0,NZ+1)
!C      CALL FFT_X_TO_FOURIER(U2,CU2,0,NY+1,0,NZ+1)
!C      CALL FFT_X_TO_FOURIER(U3,CU3,0,NY+1,0,NZ+1)
      IF (N_TH .gt. 0 ) THEN
       CALL REAL_FOURIER_TRANS_TH (.true.)
       CALL allocation_Rth (.false.) 
!C        CALL FFT_X_TO_FOURIER(TH(0,0,0,N),CTH(0,0,0,N),0,NY+1,0,NZ+1)
      ENDIF

! C Begin second step of the Fractional Step algorithm, making u divergence free
! C The following subroutine projects Uhat onto divergence free space
      CALL REM_DIV_CURV
  
 
! C Now, phi is stored in CR1X, use this to update the pressure field
! C Note, here we divide by H_BAR since it was absorbed into PHI in REM_DIV
      DO J=JSTART-1,JEND+1
        DO K=ZSTART-1,ZEND+1
          DO I=0,NX2P
           CPX(I,K,J)=CPX(I,K,J)+CR1X(I,K,J)/TEMP4
          END DO
        END DO
      END DO

!      wtime = -(wtime - omp_get_wtime ( ))

      if (rank .eq. 0) then
      open(99,file='out_screen.txt',form='formatted', &
      status='unknown', position='append')
      write(99,*) 'time taken', wtime
      close(99)

      write(6,*) 'time taken', wtime
      endif


      RETURN
      END

! C--.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE REM_DIV_CURV
! C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
use ntypes
use Domain
use Grid
use Fft_var
use TIME_STEP_VAR
use run_variable
use mg_vari, only : INIT_FLAG, BC, CALL_CHEEK_DIV
use pr_rem
use mpi_var      
implicit none

!      REAL*8 P_TMP(0:NZ+1,0:NY+1), P_TMP2(0:NZ+1,0:NY+1), &
!       RHS_TMP(0:NZ+1,0:NY+1)
!      REAL*8 RTMP(0:NXP,0:NZV-1,0:NY+1)
!      COMPLEX*16 CRTMP(0:NX2P,0:NZV-1,0:NY+1)
!      EQUIVALENCE (RTMP,CRTMP)

      REAl*8 DIV  
      INTEGER I,J,K

! c     Set BCs for phi
! c     bc(1) : left, bc(2) : right, bc(3) : bottom, bc(4) : top
! c     For periodic and dirchlet set RHS to corresponding values
! c     For neumann, set RHS to the derivative times the grid spacing
! c     normal to the wall
      bc(:)=1
      CALL_CHEEK_DIV = .TRUE. 

      call allocation_ub(.FALSE.)
 
      CR1X(:,:,:) = (0.0d0,0.0d0)
!      CRTMP(:,:,:) = (0.0d0,0.0d0)
      CS1X(:,:,:) = (0.0d0,0.0d0)
      CU2bX(:,:,:) = (0.0d0,0.0d0)
      CU3bX(:,:,:) = (0.0d0,0.0d0)
      
      DO J=JSTART,JEND+1
        DO K=ZSTART,ZEND+1
          DO I=0,NX2P
! C      U2 required to define at bottom cell face 	  
            CU2bX(I,K,J)=0.5*CJOB_22(K,J,1)*(CU2X(I,K,J)+CU2X(I,K,J-1))  &
                     + 0.5*CJOB_21(K,J,1)*(CU3X(I,K,J)+CU3X(I,K,J-1))
! C      U3 required to define at side cell face     	    
	    CU3bX(I,K,J)=0.5*CJOB_11(K,J,2)*(CU3X(I,K,J)+CU3X(I,K-1,J)) &
     	              + 0.5*CJOB_12(K,J,2)*(CU2X(I,K,J)+CU2X(I,K-1,J))
          END DO
        END DO
      END DO
      
      
      
! C Now, create the RHS vector
      DO J=2,JEND
       DO K=ZSTART,ZEND
        DO I=0,NX2P
           CR1X(I,K,J)=CIKXP(I)*INT_JACOB(K,J)*CU1X(I,K,J) &
                + (CU2bX(I,K,J+1) - CU2bX(I,K,J))        &
                + (CU3bX(I,K+1,J) - CU3bX(I,K,J))
           CS1X(I,K,J)=CR1X(I,K,J)
        END DO
       END DO
      END DO

!      if (rank .eq. np-1)then
!       DO J=2,JEND
!        DO K=ZSTART,ZEND
!         write(666,*)dble(CR1X(0,K,J))
!        enddo
!       enddo
!       close(666)
!      endif

      IF (CALL_CHEEK_DIV) THEN 

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

!C      CALL FFT_X_TO_PHYSICAL(CRTMP,RTMP,0,NY+1,0,NZ+1)

      DIV = 0.0D0
      DO J = 1,NY
       DO K = 1,NZ
        DO I = 1,NXP
         DIV = DIV + DABS(S1X(I,K,J))
        ENDDO
       ENDDO
      ENDDO
      DIV = DIV/(DBLE(NXP*NP)*DBLE(NY)*DBLE(NZ))
      CALL MPI_COMBINE_STATS(DIV,1,1)

!C      CALL FFT_X_TO_FOURIER(RTMP,CRTMP,0,NY+1,0,NZ+1)
      IF (rank .eq. 0) THEN

      write(6,*) "THE DIVERGENCE IS ", DIV, &
                " BEFORE REMOVING THE DIVERGENCE"
     
      open(99,file='out_screen.txt',form='formatted', &
      status='unknown', position='append')

      write(99,*) "THE DIVERGENCE IS ", DIV, &
                " BEFORE REMOVING THE DIVERGENCE"
      ENDIF
      ENDIF


      IF ( INIT_FLAG .EQ. .TRUE. ) THEN
       DO I = 0,NX2P
        P_TMP(:,:) = 0.0d0; RHS_TMP(:,:) = 0.0d0
        DO K = 1,NZ
         DO J = 1,NY
          RHS_TMP(K,J) = dble(CR1X(I,K,J))
! C           RHS_TMP(K,J) =1.0
         ENDDO
        ENDDO
        CALL MULTIGRID(P_TMP,RHS_TMP,I)
       ENDDO
      ENDIF
      INIT_FLAG = .FALSE.

!c!$OMP PARALLEL PRIVATE(I,K,J)
!c!$OMP  DO 
      

      DO I = 0,NX2P
! c     Solve for complex part of phi
       RHS_TMP(:,:) = 0.0d0
       DO K = 1,NZ
        DO J = 1,NY
         RHS_TMP(K,J) =-dimag(CR1X(I,K,J))
! c          RHS_TMP(K,J) = 0.1*INT_JACOB(K,J) 
        ENDDO
       ENDDO
       P_TMP(:,:) = 0.0d0
      
       CALL MULTIGRID(P_TMP,RHS_TMP,I)


! C     Solve for real part of phi
       RHS_TMP(:,:) = 0.0d0
       DO K = 1,NZ
        DO J = 1,NY
         RHS_TMP(K,J) = -dble(CR1X(I,K,J))
! c          RHS_TMP(K,J) = 1.0
        ENDDO
       ENDDO

       P_TMP2(:,:) = 0.0d0

       CALL MULTIGRID(P_TMP2,RHS_TMP,I)

 
       DO J = 0,NY+1
        DO K = 0,NZ+1
         CR1X(I,K,J) = CMPLX(P_TMP2(K,J),P_TMP(K,J))
        ENDDO
       ENDDO 

      ENDDO !END OF I LOOP 
!c!$OMP END DO
 
!c!$OMP END PARALLEL

! C Now, Solve for CUi at cell face the divergenceless velocity field


      DO J=JSTART,JEND
        DO K=ZSTART,ZEND
          DO I=0,NX2P
! c	    CU1(I,K,J) = CU1(I,K,J)- CIKXP(I)*CR1X(I,K,J)
            CU2bX(I,K,J)= CU2bX(I,K,J) &
           	- GMAT_22(K,J,1)*(CR1X(I,K,J) - CR1X(I,K,J-1)) &
               - 0.25*GMAT_12(K,J,1)*(CR1X(I,K+1,J) + CR1X(I,K+1,J-1) &
                  - CR1X(I,K-1,J) - CR1X(I,K-1,J-1))
            CU3bX(I,K,J) = CU3bX(I,K,J) &
                  - GMAT_11(K,J,2)*(CR1X(I,K,J) - CR1X(I,K-1,J)) &
                  -0.25*GMAT_12(K,J,2)*(CR1X(I,K,J+1) + CR1X(I,K-1,J+1) &
                  - CR1X(I,K,J-1) - CR1X(I,K-1,J-1))
          END DO
        END DO
      END DO

! C  Calculating velocity at the cell center from the 
! C  velocity at the cell face. 

      DO J=JSTART,JEND
        DO K=ZSTART,ZEND
          DO I=0,NX2P
	    CU1X(I,K,J) = CU1X(I,K,J)- CIKXP(I)*CR1X(I,K,J)
            CU2X(I,K,J)= CU2X(I,K,J)  &
           	   - 0.5*(CJOB_22(K,J+1,1)*(CR1X(I,K,J) + CR1X(I,K,J+1)) &
                  - CJOB_22(K,J,1)*(CR1X(I,K,J) + CR1X(I,K,J-1)))       &
                   /INT_JACOB(K,J)                                    &
                  - 0.5*(CJOB_12(K+1,J,2)*(CR1X(I,K+1,J) + CR1X(I,K,J)) &
                  - CJOB_12(K,J,2)*(CR1X(I,K,J) + CR1X(I,K-1,J)))       &
                   /INT_JACOB(K,J)
     
            CU3X(I,K,J) = CU3X(I,K,J)   &
                  - 0.5*(CJOB_11(K+1,J,2)*(CR1X(I,K,J) + CR1X(I,K+1,J)) &
                  - CJOB_11(K,J,2)*(CR1X(I,K,J) + CR1X(I,K-1,J)))       &
                    /INT_JACOB(K,J)                                   &
                  - 0.5*(CJOB_21(K,J+1,1)*(CR1X(I,K,J+1) + CR1X(I,K,J)) &
                  - CJOB_21(K,J,1)*(CR1X(I,K,J) + CR1X(I,K,J-1)))       &
                   /INT_JACOB(K,J)
          END DO
        END DO
      END DO
      
      IF ((W_BC_ZMIN .EQ. 0) .OR. (W_BC_ZMIN .EQ. 9) .OR. (W_BC_ZMIN .EQ. 10)) THEN
      ELSE
       DO I=0,NX2P
        DO J=JSTART,JEND
         CU1X(I,0,J) = CU1X(I,1,J) !- CU1X(I,NZ-1,J)
! c      CU2X(I,0,J) = CU2X(I,1,J) !- CU2X(I,NZ-1,J)
         CU3X(I,0,J) = CU3X(I,1,J) !- CU3X(I,NZ-1,J)
        ENDDO
       ENDDO
      ENDIF 


      DO I=0,NX2P
       DO J=JSTART,JEND
        CU1X(I,NZ+1,J) = CU1X(I,NZ,J) !- CU1X(I,NZ-1,J)
! c	CU2X(I,NZ+1,J) = CU2X(I,NZ,J) !- CU2X(I,NZ-1,J)
        CU3X(I,NZ+1,J) = CU3X(I,NZ,J) !- CU3X(I,NZ-1,J)
       ENDDO
      ENDDO

      DO I=0,NX2P
       DO K=ZSTART-1,ZEND+1
        CU1X(I,K,NY+1) = CU1X(I,K,NY) !- CU1X(I,NZ-1,J)
! c        CU2X(I,K,NY+1) = CU2X(I,K,NY) !- CU2X(I,NZ-1,J)
        CU3X(I,K,NY+1) = CU3X(I,K,NY) !- CU3X(I,NZ-1,J)
       ENDDO
      ENDDO
      
      
     IF (CALL_CHEEK_DIV) THEN 
! C  Now, create the RHS vector
      DO J=1,NY
       DO K=1,NZ
        DO I=0,NX2P
           CS1X(I,K,J)=CIKXP(I)*INT_JACOB(K,J)*CU1X(I,K,J) &
               + (CU2bX(I,K,J+1)-CU2bX(I,K,J))   &
               + (CU3bX(I,K+1,J)-CU3bX(I,K,J))
        END DO
       END DO
      END DO

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

!      CALL FFT_X_TO_PHYSICAL(CRTMP,RTMP,0,NY+1,0,NZ+1)

      DIV = 0.0D0
      DO J = 2,NYM
       DO K = 2,NZM
        DO I = 1,NXP
         DIV = DIV + DABS(S1X(I,K,J))
        ENDDO
       ENDDO
      ENDDO
      DIV = DIV/(DBLE(NXP*NP)*DBLE(NY)*DBLE(NZ))

      CALL MPI_COMBINE_STATS(DIV,1,1)
!C       CALL FFT_X_TO_FOURIER(RTMP,CRTMP,0,NY+1,0,NZ+1)

      IF (rank .eq. 0) THEN
      write(6,*) "THE DIVERGENCE IS ", DIV,  &
                " AFTER REMOVING THE DIVERGENCE"

      write(99,*) "THE DIVERGENCE IS ", DIV,  &
                " AFTER REMOVING THE DIVERGENCE"


      close(99)
      ENDIF
      ENDIF
!ccccccccccccccccccccccccccccccccccccccccccccccccc
!     deallocating cubs
!ccccccccccccccccccccccccccccccccccccccccccccccccc

      call allocation_ub(.TRUE.)

                 
      RETURN
      END      

! C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE RK_DUCT_3

      RETURN
      END

! C--.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE REM_DIV_DUCT


      RETURN
      END


! C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE CREATE_TH_DUCT
! C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
! C Initialize the scalar fields
! C In this subroutine, you should initialize each scalar field for the
! C particular problem of interest

use ntypes
use Domain
use Grid
!use Fft_var, only : NKX
use TIME_STEP_VAR
use run_variable
use mg_vari, only : INIT_FLAG
use mpi_var, only : RANK      
implicit none

      INTEGER I,J,K,N,MM,NN,M
      REAL*8 RNUM1,RNUM2,RNUM3,DAMP_FACT,SUM1, &
            FACT1,FACT2,FACT3,FACT4,FACT5,FACT
      REAL*8 NGY(NY), NGZ(NZ)
!      REAL*8 xpoint(0:NZ+1,0:NY+1), ypoint(0:NZ+1,0:NY+1) 


      
!      open(21,file='GRID_XZ.dat',form='formatted',status='old') 
!      xpoint(:,:) =0.0d0 
!      ypoint(:,:) =0.0d0 
      open(31,file='density_data_input.dat',form='formatted', &
      status='old')
!      DO J=0,NY+1
!       DO K=0,NZ+1
!         READ(21,*)xpoint(K,J),ypoint(K,J)
!       ENDDO
!      ENDDO
!      close(21)


      DO N=1,N_TH
       IF (CREATE_NEW_TH(N)) THEN
        write(*,*) 'A new thfield has been created'
        IF ( IC_TYPE.eq.0 ) THEN
        ELSE IF ( (IC_TYPE.eq.5).OR.(IC_TYPE.eq.6) )then
         DO J=0,NY+1
          DO K=0,NZ+1
           DO I=0,NXP
             THX(I,K,J,N)=U3_BAR(1,J) 
           END DO
          END DO
         ENDDO
        ELSE IF ( IC_TYPE.eq.7 ) THEN
         
          DO J=0,NY+1
           DO K=0,NZ+1 
            DO I=0,NXP 
             THX(I,K,J,N)=0.d0 !TH_BC_YMAX_C1(N)*(ypoint(K,J)-ypoint(1,0))
    ! &               +TH_BC_YMIN_C1(N) 
            END DO
            IF (CONT_STRAT) THEN
             IF (N .eq. 1) THEN   
!            background condition for temperature 
              TH_BAR(K,J,N)= theta_0 !(-ypoint(k,j)+ypoint(0,NY+1))
             ELSEIF (N .eq. 2) THEN
              TH_BAR(K,J,N)= Sal_0 + dSaldz*(-ypoint(k,j)+ypoint(0,NY+1))
             ENDIF
            ELSE
             READ(31,*)TH_BAR(K,J,N)
            ENDIF
           
          END DO
         ENDDO
        ENDIF
        close(31)

!      CALL REAL_FOURIER_TRANS_TH (.true.)      
!      CALL REAL_FOURIER_TRANS_Rth (.true.)
!      call allocation_Fth(.false.)
 
!C       CALL FFT_X_TO_FOURIER(TH(0,0,0,n),CTH(0,0,0,n),0,NY+1,0,NZ+1)
      ENDIF
      ENDDO   

    
      IF (N_TH .gt. 1) THEN

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

      if(rank .eq. 0) then
      call plane_parav_sponge
      endif
   

      IF (CREATE_NEW_TH(1)) THEN
      CALL REAL_FOURIER_TRANS_TH (.true.)
      CALL REAL_FOURIER_TRANS_Rth (.true.)
      call allocation_Fth(.false.)
      ENDIF

!      stop  !!!! for checking 

      RETURN
      END 







! C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE CREATE_FLOW_DUCT
! C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
use ntypes
use Domain
use Grid
use Fft_var, only : pi
use TIME_STEP_VAR
use run_variable
use mg_vari, only : INIT_FLAG
      
implicit none

      INTEGER I,J,K,MM,NN,M,N
      REAL*8 RNUM1,RNUM2,RNUM3,DAMP_FACT,SUM1, &
            FACT1,FACT2,FACT3,FACT4,FACT5,FACT
      REAL*8 NGY(NY), NGZ(NZ)
      INTEGER, DIMENSION(:), ALLOCATABLE :: seed

!      real(r8),allocatable,dimension(:,:)   :: xpoint, ypoint

!      allocate (xpoint(0:NZ+1,0:NY+1), ypoint(0:NZ+1,0:NY+1))
! C Initialize the random number generator
      CALL RANDOM_SEED(SIZE = K)
      Allocate (seed(1:K))
      seed(1:K)=10
      CALL RANDOM_SEED(PUT = seed)

! C UBULK0 and KICK should be set in input.dat

!      open(21,file='GRID_XZ.dat',form='formatted',status='old')
!      xpoint(:,:) =0.0d0
!      ypoint(:,:) =0.0d0
!      DO J=0,NY+1
!       DO K=0,NZ+1
!         READ(21,*)xpoint(K,J),ypoint(K,J)
!       ENDDO
!      ENDDO
!      close(21)

! C Set the laminar velocity profile in physical space
       write(*,*) 'UBULK0: ',UBULK0

       IF (IC_TYPE.eq.0) then
! C For closed duct flow
      mm = 5 
      nn = 5 
      fact  = LZ/LY

      DO K =1,NZ
           NGZ(K) = 2.0*GZF(K)/LZ ;
      ENDDO 
      DO J =1,NY
           NGY(J) = 2.0*GYF(J)/LY ;
      ENDDO

       DO J=1,NY
         DO K=1,NZ
          SUM1 = 0.0 ;

         DO m=0,mm
           DO n=0,nn

           fact1 = (2.0*m+1) ;
           fact2 = (2.0*n+1) ;
           fact3 = fact1*fact2**3.0 + fact**2.0*fact2*fact1**3.0 ;
           fact5 = ((-1.0)**(m+n))*fact/fact3 ;
           sum1 = sum1 + fact5*cos(fact1*pi*NGZ(K)/2)*  &
                  cos(fact2*pi*NGY(J)/2) ;
          ENDDO
         ENDDO
        
           DO I=0,NXP
             U1X(I,K,J)=UBULK0*sum1 
             U2X(I,K,J)=0.
             U3X(I,K,J)=0.
           END DO
         END DO
      END DO
      else if ((IC_TYPE.eq.1) .OR. (IC_TYPE.eq.4) ) then 
! C For open channel flow :
       DO K=0,NZM
         DO I=0,NXP
           DO J=1,NY
!            U1(I,K,J)=-(3./2.)*UBULK0*GYF(J)**2.+3.*UBULK0*GYF(J)
!            U1(I,K,J)=(-GYF(J)**2.d0+(NU+LY**2.d0)*GYF(J)/LY)/NU
             U1X(I,K,J)= (-GYF(J)**2.d0+ 2.0*LY*GYF(J))/LY**2.0
!             U1(I,K,J)=0.d0
             U2X(I,K,J)=0.
             U3X(I,K,J)=0.
           END DO
           U1X(I,K,0)=0.
           U3X(I,K,0)=0.
           U1X(I,K,NY+1)=0.
           U3X(I,K,NY+1)=0.
         END DO
      END DO
      else if (IC_TYPE.eq.2) then
! C For Couette flow:
       DO J=0,NY
         DO K=0,NZM
           DO I=0,NXP
             U1X(I,K,J)=gyf(j)
             U2X(I,K,J)=0.
             U3X(I,K,J)=0.
           END DO
         END DO
      END DO
      else if (IC_TYPE.eq.3) then
! Shear layer
       DO J=0,NY+1
         DO K=0,NZ+1
           DO I=0,NXP
             U1X(I,K,J)=TANH(GYF(J)*20.d0)
             U2X(I,K,J)=0.d0
             U3X(I,K,J)=0.d0
            END DO
          END DO
        END DO
       else if ( (IC_TYPE.eq.5).OR.(IC_TYPE.eq.6) )then        
         DO J=0,NY+1
          DO K=0,NZ+1
           DO I=0,NXP
             U1X(I,K,J)=U1_BAR(J)
             U2X(I,K,J)=0.!U2_BAR(K,J)
             U3X(I,K,J)=1!U3_BAR(1,J) 
           END DO
         END DO	 
        ENDDO
      else if (IC_TYPE.eq.7) then
	DO J=0,NY+1
          DO K=0,NZ+1
           DO I=0,NXP
              U1X(I,K,J)=0. !U1_BAR(J)
              U3X(I,K,J)=0. !U3_BAR(k,J)
	      U2X(I,K,J)=0. !U2_BAR(k,J)
     	       
           END DO
         END DO	 
        ENDDO
	UBULK0=U3X(0,1,NY)
! c	DO J=0,NY+1
! c          DO K=0,NZ+1
! c           DO I=0,NXM
! c	    U3(I,K,J)=U3(I,K,J)/30.0
! c	   END DO
! c         END DO	 
! c        ENDDO   
      end if

! C Zero the ghost cells
       IF (.NOT.USE_MPI) THEN       
        IF ( IC_TYPE == 7 ) THEN
	 DO K=0,NZ+1
          DO I=0,NXP
	   U3X(I,K,NY+1)=U3X(I,K,NY)
	  ENDDO
	 ENDDO 
	 DO J=0,NY+1
          DO I=0,NXP
	   U3X(I,0,J)=U3X(I,1,J)
	   U3X(I,NZ+1,J)=U3X(I,NZ,J)
	  ENDDO
	 ENDDO
        ELSE
	 DO K=0,NZM
          DO I=0,NXP
           U1X(I,K,0)=0.
           U2X(I,K,0)=0.
           U3X(I,K,0)=0.
           U1X(I,K,NY+1)=0.
           U2X(I,K,NY+1)=0.
           U3X(I,K,NY+1)=0.
         END DO
        END DO
        ENDIF	 
      END IF

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

!C       CALL FFT_X_TO_FOURIER(U1,CU1,0,NY+1,0,NZ+1)
!C       CALL FFT_X_TO_FOURIER(U2,CU2,0,NY+1,0,NZ+1)
!C       CALL FFT_X_TO_FOURIER(U3,CU3,0,NY+1,0,NZ+1)

 
      DO I=1,min(NX2P,NX2P_L)
        DO J=1,NY
          DO K=1,NZ
! C Now, give the velocity field a random perturbation
            CALL RANDOM_NUMBER(RNUM1)
            CALL RANDOM_NUMBER(RNUM2)
            CALL RANDOM_NUMBER(RNUM3)
            IF (IC_TYPE.eq.3) THEN
! C If we are initializing with a shear layer 
              CU1X(I,K,J)=CU1X(I,K,J) &
                  +(RNUM1-0.5)*KICK*EXP(-(GYF(J)*20.d0)**2.d0)
              CU2X(I,K,J)=CU2X(I,K,J) &
                  +(RNUM1-0.5)*KICK*EXP(-(GY(J)*20.d0)**2.d0)
              CU3X(I,K,J)=CU3X(I,K,J) &
                  +(RNUM1-0.5)*KICK*EXP(-(GYF(J)*20.d0)**2.d0)

             ELSE
              CU1X(I,K,J)=CU1X(I,K,J)+(RNUM1-0.5)*KICK
              CU2X(I,K,J)=CU2X(I,K,J)+(RNUM2-0.5)*KICK
              CU3X(I,K,J)=CU3X(I,K,J)+(RNUM3-0.5)*KICK
            END IF
          END DO

        END DO
      END DO

      DO J=0,NY+1
       DO K=0,NZ+1
         DO I=1,min(NX2P,NX2P_L)
            IF (K .EQ. 0) THEN
               DAMP_FACT = 1.0
            ELSE
               IF (IC_TYPE.eq.7) THEN
!                DAMP_FACT = exp(-10.0d0*(ypoint(K,J)-0.0d0)**2.0)
                 DAMP_FACT = exp(-1000d0*(xpoint(K,J))**2.0)
               ELSE
                DAMP_FACT = exp(-20.d0*dble((J-1))/dble((NY-1)))
               ENDIF 
            ENDIF
               CU1X(I,K,J)=CU1X(I,K,J)*DAMP_FACT
               CU2X(I,K,J)=CU2X(I,K,J)*DAMP_FACT
               CU3X(I,K,J)=CU3X(I,K,J)*DAMP_FACT
         END DO    
        END DO
       END DO
      
      WRITE(*,*) 'KICK is : ',KICK, 'Bellow j', RNUM1

!      CALL MPI_TRANSPOSE_COMPLEX_X_TO_Z(CU1X,CU1Z)
!      CU1Z(NKX+1:NX2V-1,:,:)=(0.0,0.0)
!      CALL MPI_TRANSPOSE_COMPLEX_Z_TO_X(CU1Z,CU1X)

!      CALL MPI_TRANSPOSE_COMPLEX_X_TO_Z(CU2X,CU2Z)
!      CU2Z(NKX+1:NX2V-1,:,:)=(0.0,0.0)
!      CALL MPI_TRANSPOSE_COMPLEX_Z_TO_X(CU2Z,CU2X)
      
!      CALL MPI_TRANSPOSE_COMPLEX_X_TO_Z(CU3X,CU3Z)
!      CU3Z(NKX+1:NX2V-1,:,:)=(0.0,0.0)
!      CALL MPI_TRANSPOSE_COMPLEX_Z_TO_X(CU3Z,CU3X)
      
!      deallocate (xpoint, ypoint)
! C Remove the divergence of the velocity field
      IF ( IC_TYPE .EQ. 7) THEN
       CALL SAVE_STATS_CURVI(.FALSE.)
      ELSE
       CALL SAVE_STATS_DUCT(.FALSE.)
      ENDIF 



      INIT_FLAG = .TRUE.
      IF ( IC_TYPE .EQ. 7) THEN
        CALL REM_DIV_CURV
      ELSE
        CALL REM_DIV_DUCT
      ENDIF 	

      
      
! C Get the pressure from the poisson equation
! c      CALL POISSON_P_DUCT


      IF ( IC_TYPE .EQ. 7) THEN
      CALL SAVE_STATS_CURVI(.FALSE.)
      ELSE  
      CALL SAVE_STATS_DUCT(.FALSE.)
      ENDIF
      
      
      RETURN
      END



      subroutine sponge_th(N)
! This subroutine applies a sponge relaxation (Rayleigh damping) towards a
! specified background state for the temperature field
! The intention is to allow an open boundary
use ntypes
use Domain
use Grid
!use Fft_var, only : NKX
use TIME_STEP_VAR
use run_variable
use mg_vari, only : INIT_FLAG
use mpi_var
implicit none

! The following variables will store the background state
      real*8 TH_0(0:NY+1)
      real*8 TH_TOP,b,co_b
      integer i,j,k,N,NZ_C,NY_C

! Damp fluctuation flowfield

      do j=jstart,jend
       do k=zstart,zend
        do i=1,NX2P
        CFTHX(i,k,j,N)=CFTHX(i,k,j,N) - SPONGE_SIGMA_OUT(k,j) &
                           *(CTHX(i,k,j,N)-0.d0)*INT_JACOB(K,J)

        CFTHX(i,k,j,N)=CFTHX(i,k,j,N) - SPONGE_TEMP(k,j) &
                           *(CTHX(i,k,j,N)-0.d0)*INT_JACOB(K,J)
        end do
       end do
      end do


! Damp mean flow
      IF (RANK .EQ. 0) THEN

      do j=jstart,jend
       do k=zstart,zend
        CFTHX(0,k,j,N)=CFTHX(0,k,j,N) - SPONGE_SIGMA_OUT(k,j) &
                           *(CTHX(0,k,j,N)-0.d0)*INT_JACOB(K,J)


        CFTHX(0,k,j,N)=CFTHX(0,k,j,N) - SPONGE_SIGMA(k,j) &
                           *(CTHX(0,k,j,N)-0.d0)*INT_JACOB(K,J)
       end do
      end do


      ELSE
        i = 0

      do j=jstart,jend
       do k=zstart,zend
        CFTHX(i,k,j,N)=CFTHX(i,k,j,N) - SPONGE_SIGMA_OUT(k,j) &
                           *(CTHX(i,k,j,N)-0.d0)*INT_JACOB(K,J)

        CFTHX(i,k,j,N)=CFTHX(i,k,j,N) - SPONGE_TEMP(k,j) &
                           *(CTHX(i,k,j,N)-0.d0)*INT_JACOB(K,J)
       end do
      end do

      ENDIF

      return
      end


       subroutine sponge_vel
! This subroutine applies a sponge relaxation (Rayleigh damping) towards a
! specified background state for the velocity field
! The intention is to allow an open boundary
use ntypes
use Domain
use Grid
!use Fft_var, only : NKX
use TIME_STEP_VAR
use run_variable
use mg_vari, only : INIT_FLAG
use mpi_var
implicit none

! The following variables will store the background state
      real*8 TH_0(0:NY+1)
      real*8 TH_TOP,b,co_b
      integer i,j,k,N,NZ_C,NY_C

! Damp fluctuation flowfield

      do j=jstart,jend
       do k=zstart,zend
        do i=1,NX2P
        CF1X(i,k,j)=CF1X(i,k,j)-SPONGE_SIGMA_OUT(k,j)*(CU1X(i,k,j)-0.d0) &
                           *INT_JACOB(K,J)
        CF2X(i,k,j)=CF2X(i,k,j)-SPONGE_SIGMA_OUT(k,j)*(CU2X(i,k,j)-0.d0) &
                           *INT_JACOB(K,J)
        CF3X(i,k,j)=CF3X(i,k,j)-SPONGE_SIGMA_OUT(k,j)*(CU3X(i,k,j)-0.d0) &
                           *INT_JACOB(K,J)

        CF1X(i,k,j)=CF1X(i,k,j)-SPONGE_TEMP(k,j)*(CU1X(i,k,j)-0.d0) &
                           *INT_JACOB(K,J)
        CF2X(i,k,j)=CF2X(i,k,j)-SPONGE_TEMP(k,j)*(CU2X(i,k,j)-0.d0) &
                           *INT_JACOB(K,J)
        CF3X(i,k,j)=CF3X(i,k,j)-SPONGE_TEMP(k,j)*(CU3X(i,k,j)-0.d0) &
                           *INT_JACOB(K,J)
        end do
       end do
      end do


! Damp mean flow
      IF (RANK .EQ. 0) THEN

      do j=jstart,jend
       do k=zstart,zend
        CF1X(0,k,j)=CF1X(0,k,j)-SPONGE_SIGMA_OUT(k,j)*(CU1X(0,k,j)-0.d0)  &
                           *INT_JACOB(K,J)
        CF2X(0,k,j)=CF2X(0,k,j)-SPONGE_SIGMA_OUT(k,j)*(CU2X(0,k,j)-0.d0)  &
                           *INT_JACOB(K,J)
!        CF3X(0,k,j)=CF3X(0,k,j)-SPONGE_SIGMA_OUT(k,j)*(CU3X(0,k,j)-0.d0)  &
!                           *INT_JACOB(K,J)


        CF1X(0,k,j)=CF1X(0,k,j)-SPONGE_SIGMA(k,j)*(CU1X(0,k,j)-0.d0)  &
                           *INT_JACOB(K,J)
        CF2X(0,k,j)=CF2X(0,k,j)-SPONGE_SIGMA(k,j)*(CU2X(0,k,j)-0.d0)  &
                           *INT_JACOB(K,J)
!        CF3X(0,k,j)=CF3X(0,k,j)-SPONGE_SIGMA(k,j)*(CU3X(0,k,j)-0.d0)  &
!                           *INT_JACOB(K,J)
       end do
      end do

      ELSE
        i = 0

      do j=jstart,jend
       do k=zstart,zend
        CF1X(i,k,j)=CF1X(i,k,j)-SPONGE_SIGMA_OUT(k,j)*(CU1X(i,k,j)-0.d0) &
                           *INT_JACOB(K,J)
        CF2X(i,k,j)=CF2X(i,k,j)-SPONGE_SIGMA_OUT(k,j)*(CU2X(i,k,j)-0.d0) &
                           *INT_JACOB(K,J)
        CF3X(i,k,j)=CF3X(i,k,j)-SPONGE_SIGMA_OUT(k,j)*(CU3X(i,k,j)-0.d0) &
                           *INT_JACOB(K,J)

        CF1X(i,k,j)=CF1X(i,k,j)-SPONGE_TEMP(k,j)*(CU1X(i,k,j)-0.d0) &
                           *INT_JACOB(K,J)
        CF2X(i,k,j)=CF2X(i,k,j)-SPONGE_TEMP(k,j)*(CU2X(i,k,j)-0.d0) &
                           *INT_JACOB(K,J)
        CF3X(i,k,j)=CF3X(i,k,j)-SPONGE_TEMP(k,j)*(CU3X(i,k,j)-0.d0) &
                           *INT_JACOB(K,J)
       end do
      end do

      ENDIF

      return
      end
!------------------------------------------------------------------.----|--------|


      
! C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE VIS_FLOW_DUCT
! C----*|--.---------.---------.---------.---------.---------.---------.-|-------|

      RETURN
      END



           
      subroutine courant
! This subroutine sets the timestep based on the specified CFL number
! The subroutine should be called with the velocity in physical space


      return
      end

      
      subroutine courant_curv
! This subroutine sets the timestep based on the specified CFL number
! The subroutine should be called with the velocity in physical space

use Domain
use Grid
!use Fft_var, only : NKX
use TIME_STEP_VAR
use run_variable
use mg_vari, only : INIT_FLAG
      
implicit none

      real*8 vel
      real*8 dt
      real*8 dt_x,dt_y,dt_z,min_x,dt_dif
      integer i,j,k
      integer imin,jmin,kmin

! Set the initial dt to some arbitrary large number
      dt=999.d0

!    dt based on Diffusion no 
      min_x = LX/dble(NX)       
      dt_dif    = DFN*min_x**2/NU
     
      do j=1,NY
        do k=1,NZ
          do i=0,NXP
! C      U2 required to define at bottom cell face 	  
            U2bX(I,K,J)=0.5*CJOB_22(K,J,1)*(U2X(I,K,J)+U2X(I,K,J-1)) &
                     + 0.5*CJOB_21(K,J,1)*(U3X(I,K,J)+U3X(I,K,J-1))
! C      U3 required to define at side cell face     	    
	    U3bX(I,K,J)=0.5*CJOB_11(K,J,2)*(U3X(I,K,J)+U3X(I,K-1,J)) &
     	              + 0.5*CJOB_12(K,J,2)*(U2X(I,K,J)+U2X(I,K-1,J))
     
            dt_x=cfl*dx(i)/abs(U1X(i,k,j))
            dt_y=cfl*INT_JACOB(K,J)/abs(U2bX(i,k,j))
            dt_z=cfl*INT_JACOB(K,J)/abs(U3bX(i,k,j))
! c            if( dt_z .le. 0 ) then
! c             write(6,*) 'dt=', dt_x, 'i,j,k=',i,j,k,dz(k),abs(U3(i,k,j))
! c            endif       
            dt=min(dt,dt_x,dt_y,dt_z)
          end do
        end do
      end do

! c      stop
      
      if (dt.le.0) then
        write(*,*) 'Error: dt<=0 in courant'
! Set DELTA_T to some small default value
        write(*,*) dt
! c        DELTA_T=0.0001d0
      else if (dt.ge.0.002) then
         write(*,*) 'WARNING: DELTA_T > 0.002, value capped at 0.002, dt=', dt
! c        DELTA_T=0.0001d0
      else
        DELTA_T=dt
        H_BAR(1)=DELTA_T*(8.0/15.0)
        H_BAR(2)=DELTA_T*(2.0/15.0)
        H_BAR(3)=DELTA_T*(5.0/15.0)
      end if

      return
      end   

      subroutine courant_curv_mpi
! This subroutine sets the timestep based on the specified CFL number
! The subroutine should be called with the velocity in physical space

use Domain
use Grid
!use Fft_var, only : NKX
use TIME_STEP_VAR
use run_variable
use mg_vari, only : INIT_FLAG
use mpi_var, only : rank      
implicit none

      real*8 vel
      real*8 dt
      real*8 dt_x,dt_y,dt_z,min_x,dt_dif
      integer i,j,k
      integer imin,jmin,kmin

! Set the initial dt to some arbitrary large number
      dt=999.d0

!    dt based on Diffusion no 
      min_x = LX/dble(NX)       
      dt_dif    = DFN*min_x**2/NU
     
      do j=1,NY
        do k=1,NZ
          do i=0,NXP
! C      U2 required to define at bottom cell face 	  
            U2bX(I,K,J)=0.5*CJOB_22(K,J,1)*(U2X(I,K,J)+U2X(I,K,J-1)) &
                     + 0.5*CJOB_21(K,J,1)*(U3X(I,K,J)+U3X(I,K,J-1))
! C      U3 required to define at side cell face     	    
	    U3bX(I,K,J)=0.5*CJOB_11(K,J,2)*(U3X(I,K,J)+U3X(I,K-1,J)) &
     	              + 0.5*CJOB_12(K,J,2)*(U2X(I,K,J)+U2X(I,K-1,J))
     
            dt_x=cfl*dx(i)/abs(U1X(i,k,j))
            dt_y=cfl*INT_JACOB(K,J)/abs(U2bX(i,k,j))
            dt_z=cfl*INT_JACOB(K,J)/abs(U3bX(i,k,j))
! c            if( dt_z .le. 0 ) then
! c             write(6,*) 'dt=', dt_x, 'i,j,k=',i,j,k,dz(k),abs(U3(i,k,j))
! c            endif       
            dt=min(dt,dt_x,dt_y,dt_z)
          end do
        end do
      end do

! c      stop
      CALL MPI_COURANT(dt)

      if (dt.le.0) then
        write(*,*) 'Error: dt<=0 in courant'
! Set DELTA_T to some small default value
        write(*,*) dt
! c        DELTA_T=0.0001d0
      else if (dt.ge.DELTA_T_in) then
       if (rank .eq. 0) then  
        write(*,*) 'WARNING: DELTA_T > ',DELTA_T_in, ', value capped at',DELTA_T_in,',dt=', dt
       endif
! c        DELTA_T=0.0001d0
      else
        DELTA_T=dt
        H_BAR(1)=DELTA_T*(8.0/15.0)
        H_BAR(2)=DELTA_T*(2.0/15.0)
        H_BAR(3)=DELTA_T*(5.0/15.0)
      end if

      return
      end   



! C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE VEL_PROFILE
! C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      use ntypes
use Domain
use Grid
use Fft_var
use TIME_STEP_VAR
use run_variable
use mg_vari, only : INIT_FLAG
      
implicit none

      INTEGER I,J,K 
!      real*8  xpoint(0:NZ+1,0:NY+1), ypoint(0:NZ+1,0:NY+1), &
      real*8  rhs_w(0:nz+1,0:ny+1),phi_int(0:NZ+1,0:NY+1)       
      real*8  div,Q,H1, H2,sigma,omega_0,din1,z1,y1 

!      open(21,file='GRID_XZ.dat',form='formatted',status='old')
!      xpoint(:,:) =0.0d0
!      ypoint(:,:) =0.0d0
!      DO J=0,NY+1
!       DO K=0,NZ+1
!         READ(21,*)xpoint(K,J),ypoint(K,J)
!       ENDDO
!      ENDDO
!      close(21)

111   format(2f16.12)
      
      DO J=0,NY+1
        DO K=0,NZ+1
           U3_BAR(K,J) = 0.0
           U2_BAR(K,J) = 0.0
	   phi_int(k,j)= 0.0
           u2bx(:,k,j)  = 0.d0
           u3bx(:,k,j)  = 0.d0
        ENDDO
      ENDDO
     
      H1 = (ypoint(1,NY) - ypoint(1,1))
      Q =  0.d0 !-(PX0/(12.d0*NU))*H1**3.0
      H2 = (ypoint(NZ,NY) - ypoint(NZ,1))
            

      DO J=1,NY
       DO K=0,1
        U3_BAR(K,J)=6.d0*(Q/H1**3.0)*( (ypoint(1,NY)+ypoint(1,1) &
            - ypoint(1,J) )*ypoint(1,J) &
            - ypoint(1,NY)*ypoint(1,1) )
       ENDDO
       
       DO K=NZ-1,NZ+1
        U3_BAR(K,J)=6.d0*(Q/H2**3.0)*( (ypoint(NZ,NY)+ypoint(NZ,1) &
            - ypoint(NZ,J) )*ypoint(NZ,J)                         &
            - ypoint(NZ,NY)*ypoint(NZ,1) )
       ENDDO
      ENDDO      
       
      DO J=1,NY
       DO K=2,NZ-2
       H1 = (ypoint(k,NY) - ypoint(k,1))
        U3_BAR(K,J)=6.d0*(Q/H1**3.0)*( (ypoint(k,NY)+ypoint(k,1) &
            - ypoint(k,J) )*ypoint(k,J)                         &
            - ypoint(k,NY)*ypoint(k,1) )  
       ENDDO
       ENDDO     
 
       z1=3.0 
       y1=1.5 
       omega_0 = 0.5 
       sigma = 0.3	 

! c       DO J=2,NY-1
! c       DO K=2,NZ+1
! c        din1 = (ypoint(k,j)-y1)**2.0 + (xpoint(k,j)-z1)**2.0 ;
! c        U2_BAR(K,J)=U2_BAR(K,J)-(sigma**2.0*omega_0*(xpoint(k,j)-z1))
! c     &   *( 1.0 - exp(-(din1/(2.0*sigma**2.0))) )/din1
! c        U3_BAR(K,J)=U3_BAR(K,J)+(sigma**2.0*omega_0*(ypoint(k,j)-y1))
! c     &   *( 1.0 - exp(-(din1/(2.0*sigma**2.0))) )/din1
! 
! c       ENDDO
! c       ENDDO

      do j=1,NY
        do k=1,NZ
          do i=0,NXP
! C      U2 required to define at bottom cell face        
            U2bx(I,K,J)=0.5*CJOB_22(K,J,1)*(U2_bar(K,J)+U2_bar(K,J-1)) &
                     + 0.5*CJOB_21(K,J,1)*(U3_bar(K,J)+U3_bar(K,J-1))
! C      U3 required to define at side cell face          
            U3bx(I,K,J)=0.5*CJOB_11(K,J,2)*(U3_bar(K,J)+U3_bar(K-1,J)) &
                     + 0.5*CJOB_12(K,J,2)*(U2_bar(K,J)+U2_bar(K-1,J))

           enddo
         enddo
       enddo
       
       
        DO J=1,NY  
	 do i=0,NXP
	  U3bx(I,0,J) = U3bx(I,1,J)
	 enddo 
        enddo
      
      
!       
! c      call SOR (phi_int,rhs_w)
!       
! c      DO J=2,NY-1
! c        DO K=2,NZ
! c         U2_bar(K,J)=0.5*(CJOB_22(K,J+1,1)*(phi_int(K,J)+phi_int(K,J+1))
! c     &             - CJOB_22(K,J,1)*(phi_int(K,J)+phi_int(K,J-1)))
! c     &              /INT_JACOB(K,J)
! c     &             - 0.5*(CJOB_12(K+1,J,2)*(phi_int(K+1,J)+phi_int(K,J))
! c     &             - CJOB_12(K,J,2)*(phi_int(K,J)+phi_int(K-1,J)))
! c     &              /INT_JACOB(K,J)
! 
! c         U3_bar(K,J)=0.5*(CJOB_11(K+1,J,2)*(phi_int(K,J)+phi_int(K+1,J))
! c     &             - CJOB_11(K,J,2)*(phi_int(K,J)+phi_int(K-1,J)))
! c     &               /INT_JACOB(K,J)
! c     &             - 0.5*(CJOB_21(K,J+1,1)*(phi_int(K,J+1)+phi_int(K,J))
! c     &             - CJOB_21(K,J,1)*(phi_int(K,J)+phi_int(K,J-1)))
! c     &              /INT_JACOB(K,J)
! c        END DO
! c      END DO

      do j=1,NY
        do k=1,NZ
          do i=0,NXP
! C      U2 required to define at bottom cell face        
            U2bx(I,K,J)=0.5*CJOB_22(K,J,1)*(U2_bar(K,J)+U2_bar(K,J-1)) &
                     + 0.5*CJOB_21(K,J,1)*(U3_bar(K,J)+U3_bar(K,J-1))
! C      U3 required to define at side cell face          
            U3bx(I,K,J)=0.5*CJOB_11(K,J,2)*(U3_bar(K,J)+U3_bar(K-1,J)) &
                     + 0.5*CJOB_12(K,J,2)*(U2_bar(K,J)+U2_bar(K-1,J))

           enddo
         enddo
       enddo

! c       DO J=1,NY
! c        DO K=1,NZ
! c          DO I=0,NXM
! c            U2b(I,K,J)= GMAT_22(K,J,1)*(phi_int(K,J)-phi_int(K,J-1))
! c     &          + 0.25*GMAT_12(K,J,1)*(phi_int(K+1,J)+phi_int(K+1,J-1)
! c     &             - phi_int(K-1,J)-phi_int(K-1,J-1))
! c            U3b(I,K,J) = GMAT_11(K,J,2)*(phi_int(K,J)-phi_int(K-1,J))
! c     &          + 0.25*GMAT_12(K,J,2)*(phi_int(K,J+1)+phi_int(K-1,J+1)
! c     &             - phi_int(K,J-1)-phi_int(K-1,J-1))
! c          END DO
! c        END DO
! c      END DO


      S1X(:,:,:)=0.d0

      DO J=2,NY-1
       DO K=2,NZ-1
          do i=0,NXP
           S1X(I,K,J)=  (U2bx(I,K,J+1) - U2bx(I,K,J)) &
                + (U3bx(I,K+1,J) - U3bx(I,K,J))
        END DO
       END DO
      END DO

      DIV = 0.0D0
      DO J = 1,NY
       DO K = 1,NZ
        do i=0,NXP
           S1X(I,K,J)=  (U2bx(I,K,J+1) - U2bx(I,K,J)) &
                + (U3bx(I,K+1,J) - U3bx(I,K,J))
        END DO
       END DO
      END DO

      DIV = 0.0D0
      DO J = 1,NY
       DO K = 1,NZ
        DO I = 1,NX
         DIV = DIV + DABS(S1X(I,K,J))
        ENDDO
       ENDDO
      ENDDO
      DIV = DIV/(DBLE(NX)*DBLE(NY)*DBLE(NZ))

      
      write(6,*) "THE DIVERGENCE IS ", DIV, &
                "At INITIAL FLOW field"
 

      open(10,file='sample.dat',form='formatted',status='unknown')
      write(10,*) 'title = "sample mesh" '
      write(10,*) 'variables = "x", "y" "A" '
      write(10,*) 'zone i=',NZ+2 ,', j=', NY+2,'DATAPACKING=POINT'

      do j=0,NY+1
       do i=0,NZ+1
         write(10,113) xpoint(i,j),ypoint(i,j), phi_int(i,j)
       enddo
      enddo
113   format(3f16.8)
      close(10)
    

      
! c      stop 

      RETURN
      END  




       
! CCC BLANK SUBROUTINE FOR COMPILATION
! C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE RK_CHAN_1(FINAL)
! C----*|--.---------.---------.---------.---------.---------.---------.-|-------|

      return
      end

! C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE RK_CHAN_2(FINAL)
! C----*|--.---------.---------.---------.---------.---------.---------.-|-------|

      return
      end

! C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE RK_CHAN_3(FINAL)
! C----*|--.---------.---------.---------.---------.---------.---------.-|-------|

      return
      end

! C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE REM_DIV_CHAN
! C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
       

      RETURN
      END

! C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE POISSON_P_CHAN
! C----*|--.---------.---------.---------.---------.---------.---------.-|-------|

      RETURN
      END

! C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE create_flow_chan
! C----*|--.---------.---------.---------.---------.---------.---------.-|-------|

      RETURN
      END

! C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE INIT_CHAN

      RETURN
      END

subroutine filter_all_real
      
use ntypes
use Domain
use Grid
use Fft_var
use TIME_STEP_VAR
use run_variable

implicit none

      INTEGER I,J,K,N
      
              ! C convert to physical space.


      CS1X=CU1X
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
      call  filter_chan


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
      CU1X=CS1X


      CS1X=CU2X
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
      call  filter_chan


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
      CU2X=CS1X
            
      CS1X=CU3X
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
      call  filter_chan


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
      CU3X=CS1X
      
      return
      end


! C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
        SUBROUTINE filter_chan
! C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
use ntypes
use Domain
use Grid
use run_variable
use mpi_var, only : RANK      
implicit none     

      INTEGER I,J,K
      integer km2(0:NZ+1),km1(0:NZ+1),kp1(0:NZ+1),kp2(0:NZ+1)

! These are the weights for the filtering operation used
      real*8 W0,W1,W2,Wm1,Wm2,Wm1_j,W0_j,W1_j,alf_y(0:NY+1), alf_z(0:NZ+1)

! filter type = 1: grid filter operation
! filter type = 2: test filter operation
! filter type = 3: grid filter operation designed for ABL problem
! The following is for the 3-point trapezoidal rule, alpha*beta=sqrt(6)
      
      alf_y=0.0d0 
      alf_z=0.0d0
!       open(321,file='alf2.dat',status='unknown',form='formatted')
      If ( filter_type .eq. 1 ) then
       Wm2=0.d0
       Wm1=1.d0/8.d0
       W0=6.d0/8.d0
       W1=1.d0/8.d0
       W2=0.d0
      elseif ( filter_type .eq. 2 ) then
       Wm2=0.d0
       Wm1=1.d0/4.d0
       W0=1.d0/2.d0
       W1=1.d0/4.d0
       W2=0.d0 
       elseif ( filter_type .eq. 3 ) then
       Wm2=0.d0
       Wm1=1.d0/4.d0
       W0=1.d0/2.d0
       W1=1.d0/4.d0
       W2=0.d0

       do j=0,NY+1
       alf_Y(j)=1.0d0-(0.5*(1+tanh(0.075*(ypoint(1,J)+50))))
       if (ypoint(1,J) > 0 ) alf_Y(j) = 0.0d0
       If (RANK == 0) THEN
!        write(113,*)ypoint(1,J),alf_Y(j)
       ENDIF
       enddo

       do k=0,NZ+1
       alf_Z(k)=1.0d0-(0.5*(1+tanh(0.03*(xpoint(K,1)-60))))
       if (xpoint(k,1) > 200.0d0 ) alf_Z(k) = 0.d0
       If (RANK == 0) THEN
!        write(114,*)xpoint(K,1),alf_Z(K)
       ENDIF  
       enddo
!        write(*,*) 'stopped'
!        close(321)
!        stop
       else
       pause 'Error, unsupported LES_FILTER_TYPE chosen'
      endif


! Filter in the periodic directions using cshift
! Note, cshift if not used since it appears to be slower
! Apply filter in the x-direction
!      B=Wm2*CSHIFT(B,-2,1)+Wm1*CSHIFT(B,-1,1)+W0*B+W1*CSHIFT(B,1,1)
!     &       +W2*CSHIFT(B,2,1)

! Filter using more efficient F77 syntax:
! Set up array to loop around periodic directions
 
      do k=2,NZ+1
        km2(k)=k-2
      end do
      km2(1)=NZ+1
      km2(0)=NZ
      
      do k=1,NZ+1
        km1(k)=k-1
      end do
      km1(0)=NZ+1
      
      do k=0,NZ
        kp1(k)=k+1
      end do
      kp1(NZ+1)=0
      
      do k=0,NZ-1
        kp2(k)=k+2    
      end do
      kp2(NZ)=0
      kp2(NZ+1)=1

      DO J=JSTART,JEND-1
        DO K=ZSTART,ZEND-1
          DO I=0,NXP
            S2X(i,k,j)=Wm2*S1X(i,km2(k),j)+Wm1*S1X(i,km1(k),j)+W0*S1X(i,k,j) &
              +W1*S1X(i,kp1(k),j)+W2*S1X(i,kp2(k),j)
          end do
        end do  
      end do

      DO J=JSTART,JEND-1
        DO K=ZSTART,ZEND-1
          DO I=0,NXP
            S1X(i,k,j) = (1.0d0-alf_y(J)*alf_z(K))*S1X(i,k,j)+alf_z(K)* alf_y(J)*S2X(i,k,j)
          end do
        end do
      end do
      
      RETURN
      END

! ! C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE create_TH_chan
! ! C----*|--.---------.---------.---------.---------.---------.---------.-|-------|

      RETURN
      END

! ! C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE SAVE_STATS_CHAN(FINAL)
! ! C----*|--.---------.---------.---------.---------.---------.---------.-|-------|

      return
      end 
