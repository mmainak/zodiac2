
      
! C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE THOMAS_REAL_SP(A,B,C,G,D,K,NY,NX)
! C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
! C Uses the Thomas algorithm to solve Ax=b for tridiagonal A
! C The RHS vector and solution are real
! C Input lower, main, and upper diagonals, ld, md, ud, and rhs x
! C Returns solution in x
! C The indexing should be done by ROW, ie.
! C [ b1  c1   0   0   0 ...
! C [ a2  b2  c2   0   0 ...
! C [  0  a3  b3   c3  0 ...

      INTEGER I, J, NX, NY,k
      REAL*8 A(0:NX,0:NY), B(0:NX,0:NY), C(0:NX,0:NY), G(0:NX,0:NY)
      REAL*8 D

       J= K
! c       write(6,*) 'I am here', J
       DO I=0,NX
! c          D=-D/A(I,J+1)
          A(I,J+1)=A(I,J+1)-B(I,J)*D/A(I,J)
          B(I,J+1)=B(I,J+1)-C(I,J)*D/A(I,J)
          G(I,J+1)=G(I,J+1)-G(I,J)*D/A(I,J)
       END DO

! c       write(6,*) 'I am here after', J

      DO J=0,NY-1
        DO I=0,NX
          A(I,J+1)=-A(I,J+1)/B(I,J)
          B(I,J+1)=B(I,J+1)+A(I,J+1)*C(I,J)
          G(I,J+1)=G(I,J+1)+A(I,J+1)*G(I,J)
        END DO
      END DO
      DO I=0,NX
        G(I,NY)=G(I,NY)/B(I,NY)
      END DO
      DO J=NY-1,0,-1
        DO I=0,NX
          G(I,J)=(G(I,J)-C(I,J)*G(I,J+1))/B(I,J)
        END DO
      END DO

      RETURN
      END




! C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE THOMAS_REAL(A,B,C,G,NY,NX)
! C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
! C Uses the Thomas algorithm to solve Ax=b for tridiagonal A
! C The RHS vector and solution are real
! C Input lower, main, and upper diagonals, ld, md, ud, and rhs x
! C Returns solution in x
! C The indexing should be done by ROW, ie.
! C [ b1  c1   0   0   0 ...
! C [ a2  b2  c2   0   0 ...
! C [  0  a3  b3   c3  0 ...

      INTEGER I, J, NX, NY
      REAL*8 A(0:NX,0:NY), B(0:NX,0:NY), C(0:NX,0:NY), G(0:NX,0:NY)

      DO J=0,NY-1
        DO I=0,NX
          A(I,J+1)=-A(I,J+1)/B(I,J)
          B(I,J+1)=B(I,J+1)+A(I,J+1)*C(I,J)
          G(I,J+1)=G(I,J+1)+A(I,J+1)*G(I,J)
        END DO
      END DO
      DO I=0,NX
        G(I,NY)=G(I,NY)/B(I,NY)
      END DO
      DO J=NY-1,0,-1
        DO I=0,NX
          G(I,J)=(G(I,J)-C(I,J)*G(I,J+1))/B(I,J)
        END DO
      END DO

      RETURN
      END

! C----*|--.---------.---------.---------.---------.---------.---------.-|-------|    
      SUBROUTINE THOMAS_COMPLEX(A,B,C,G,NY,NX)
! C----*|--.---------.---------.---------.---------.---------.---------.-|-------|

! C Uses the Thomas algorithm to solve Ax=b for tridiagonal A
! C The RHS vector and solution is complex
! C Input lower, main, and upper diagonals, ld, md, ud, and rhs x
! C Returns solution in x
! C The indexing should be done by ROW, ie.
! C [ b1  c1   0   0   0 ...
! C [ a2  b2  c2   0   0 ...
! C [  0  a3  b3   c3  0 ...

      INTEGER I, J, NY, NX
      REAL*8 A(0:NX,0:NY), B(0:NX,0:NY), C(0:NX,0:NY)
      COMPLEX*16 G(0:NX,0:NY)

      DO J=0,NY-1
        DO I=0,NX
          A(I,J+1)=-A(I,J+1)/B(I,J)
          B(I,J+1)=B(I,J+1)+A(I,J+1)*C(I,J)
          G(I,J+1)=G(I,J+1)+A(I,J+1)*G(I,J)
        END DO
      END DO
      DO I=0,NX
        G(I,NY)=G(I,NY)/B(I,NY)
      END DO
      DO I=0,NX
        DO J=NY-1,0,-1
          G(I,J)=(G(I,J)-C(I,J)*G(I,J+1))/B(I,J)
        END DO
      END DO

      RETURN
      END


! C----*|--.---------.---------.---------.---------.---------.---------.-|------
      SUBROUTINE APPLY_BC_TH_LOWER(N,K)
! C----*|--.---------.---------.---------.---------.---------.---------.-|-----
use ntypes
use Domain
use Grid
! use Fft_var
use TIME_STEP_VAR
use run_variable
use ADI_var
use mg_vari, only : INIT_FLAG
      
implicit none
      INTEGER I,N,K

! C Bottom Wall:
      IF ( TH_BC_YMIN(N) .EQ.0) THEN
! C Dirichlet
        DO I=0,NXP
          MATLY(I,0)=0. 
          MATDY(I,0)=1.
          MATUY(I,0)=0.                   
          VECY(I,0)=0.

          MATLY(I,1)=0. 
          MATDY(I,1)=1.
          MATUY(I,1)=0.                   
          VECY(I,1)=TH_BC_YMIN_C1(N) 
        END DO
      ELSE IF (TH_BC_YMIN(N) .eq. 1) THEN
! C Neumann
        DO I=0,NXP
          MATLY(I,1)=0.
          MATDY(I,1)=-1.
          MATUY(I,1)=1.
          VECY(I,1)= -(TH_BAR(K,2,N)-TH_BAR(K,1,N))
!0.5d0*TH_BC_YMIN_C1(N)*  &
!               (INT_JACOB(K,1)+INT_JACOB(K,2))/CJOB_22(K,2,1)
        END DO

        DO I=0,NXP
          MATLY(I,0)=0.
          MATDY(I,0)=-1.
          MATUY(I,0)=1.
          VECY(I,0)=-(TH_BAR(K,1,N)-TH_BAR(K,0,N))
!0.5d0*TH_BC_YMIN_C1(N)* &
!               (INT_JACOB(K,0)+INT_JACOB(K,1))/CJOB_22(K,1,1)
        END DO 
      ElSE IF (TH_BC_YMIN(N) .EQ. 6) THEN
! C Periodic
        DO I=0,NXP
          MATLY(I,0)=0. 
          MATDY(I,0)=1.
          MATUY(I,0)=0.                   
          VECY(I,0)=THX(I,K,NY,N)
	ENDDO        	 
      ElSE IF (TH_BC_YMIN(N) .EQ. 7) THEN
! C Mixed boundary condition for w at the bottom wall 
       IF ( K .lE. 80 ) THEN
        DO I=0,NXP
          MATLY(I,1)=0.
          MATDY(I,1)=-1.
          MATUY(I,1)=1.
          VECY(I,1)=DY(2)*TH_BC_YMIN_C1(N)          

          MATLY(I,0)=0.
          MATDY(I,0)=-1.
          MATUY(I,0)=1.
          VECY(I,0)=DY(1)*TH_BC_YMIN_C1(N)
        END DO
       ELSE
        DO I=0,NXP
          MATLY(I,0)=0.
          MATDY(I,0)=1.
          MATUY(I,0)=0.
          VECY(I,0)=0.

          MATLY(I,1)=0.
          MATDY(I,1)=1.
          MATUY(I,1)=0.
          VECY(I,1)=TH_BC_YMIN_C1(N)
        END DO   
       END IF   
      ELSE IF (TH_BC_YMIN(N) .eq. 8) THEN
! C Neumann when sponge is used at left and right boundary
        DO I=0,NXP
          MATLY(I,1)=0.
          MATDY(I,1)=-1.
          MATUY(I,1)=1.
! c          VECY(I,1)=(0.5d0*TH_BC_YMIN_C1(N)*
! c     &          (INT_JACOB(K,1)+INT_JACOB(K,2))/CJOB_22(K,2,1))
! c     &          /(1.d0 + SPONGE_SIGMA_OUT(K))
          VECY(I,1)=-(TH_BAR(K,2,N)-TH_BAR(K,1,N)) &
               /(1.d0 + SPONGE_SIGMA_OUT(K,0))
        END DO

        DO I=0,NXP
          MATLY(I,0)=0.
          MATDY(I,0)=-1.
          MATUY(I,0)=1.
! c          VECY(I,0)=(0.5d0*TH_BC_YMIN_C1(N)*
! c     &          (INT_JACOB(K,0)+INT_JACOB(K,1))/CJOB_22(K,1,1))
! c     &          /(1.d0 + SPONGE_SIGMA_OUT(K)) 
          VECY(I,0)=-(TH_BAR(K,1,N)-TH_BAR(K,0,N)) &
               /(1.d0 + SPONGE_SIGMA_OUT(K,0)) 
        END DO
      ELSE
        write(*,*) 'WARNING: TH_BC_LOWER is of unknown type'
      END IF

      RETURN 
      END



! C----*|--.---------.---------.---------.---------.---------.---------.-|----
      SUBROUTINE APPLY_BC_TH_UPPER(N,K)
! C----*|--.---------.---------.---------.---------.---------.---------.-|--
use ntypes
use Domain
use Grid
! use Fft_var
use TIME_STEP_VAR
use run_variable
use ADI_var
use mg_vari, only : INIT_FLAG
      
implicit none

      INTEGER I, N,K

! C Top wall
      IF (TH_BC_YMAX(N) .EQ. 0) THEN
! C Dirichlet
        DO I=0,NXP
          MATLY(I,NY+1)=0.
          MATDY(I,NY+1)=1.
          MATUY(I,NY+1)=0.
          VECY(I,NY+1)=TH_BC_YMAX_C1(N)!TH_BC_ZMIN_C1(N)*SPONGE_temp(K,0)

          MATLY(I,NY)=0.
          MATDY(I,NY)=1.
          MATUY(I,NY)=0.
          VECY(I,NY)=TH_BC_YMAX_C1(N)!TH_BC_ZMIN_C1(N)*SPONGE_temp(K,0)
        END DO
      ELSE IF (TH_BC_YMAX(N) .eq. 1) THEN
! C Neumann
        DO I=0,NXP
          MATLY(I,NY)=-1.
          MATDY(I,NY)=1.
          MATUY(I,NY)=0.
          VECY(I,NY)=0.5d0*TH_BC_YMAX_C1(N)* &
               (INT_JACOB(K,NY)+INT_JACOB(K,NY-1))/CJOB_22(K,NY-1,1)
        END DO
        DO I=0,NXP
          MATLY(I,NY+1)=-1
          MATDY(I,NY+1)=1.
          MATUY(I,NY+1)=0.
          VECY(I,NY+1)=0.5d0*TH_BC_YMAX_C1(N)* &
               (INT_JACOB(K,NY+1)+INT_JACOB(K,NY))/CJOB_22(K,NY,1)
        END DO 
      ELSE IF (TH_BC_YMAX(N) .eq. 6) THEN	
! c periodic 
       DO I=0,NXP
          MATLY(I,NY+1)=0.
          MATDY(I,NY+1)=1.
          MATUY(I,NY+1)=0.
          VECY(I,NY+1)=THX(I,K,1,N)
       ENDDO	  
      
      END IF

      RETURN
      END   

! C----*|--.---------.---------.---------.---------.---------.---------.-|------
      SUBROUTINE APPLY_BC_TH_LEFT(N,J)
! C----*|--.---------.---------.---------.---------.---------.---------.-|-----
use ntypes
use Domain
use Grid
! use Fft_var
use TIME_STEP_VAR
use run_variable
use ADI_var
use forcing
use mg_vari, only : INIT_FLAG
      
implicit none
      INTEGER I,N,J

! C Left Wall:
      IF (TH_BC_ZMIN(N) .EQ.0) THEN
! C Dirichlet
       IF (IC_TYPE .EQ. 5) THEN
        DO I=0,NXP
          MATLX(I,0)=0.
          MATDX(I,0)=1.
          MATUX(I,0)=0.
          VECX(I,0)=0.

          MATLX(I,1)=0.
          MATDX(I,1)=1.
          MATUX(I,1)=0.
          VECX(I,1)=U3_BAR(1,J)
        END DO
       ElSE  
        DO I=0,NXP          
          MATLX(I,0)=0. 
          MATDX(I,0)=1.
          MATUX(I,0)=0.                   
          VECX(I,0)=TH_BC_ZMIN_C1(N)*((THX(I,1,J,2)+Sal_0)*m_FP -theta_0) !/(1.0d0 + 10.0d0*SPONGE_SIGMA(0,J))

          MATLX(I,1)=0. 
          MATDX(I,1)=1.
          MATUX(I,1)=0.                   
          VECX(I,1)=TH_BC_ZMIN_C1(N)*((THX(I,1,J,2)+Sal_0)*m_FP-theta_0)  !/(1.0d0 + 10.0d0*SPONGE_SIGMA(0,J)) 
        END DO
       ENDIF
      ELSE IF (TH_BC_ZMIN(N).eq.1) THEN
! C Neumann
! c        DO I=0,NXP
! c          MATLX(I,1)=0.
! c          MATDX(I,1)=-1.
! c          MATUX(I,1)=1.
! c          VECX(I,1)=TH_BC_ZMIN_C1(N)
! c        END DO
       IF(MELTING_MODEL)THEN
        DO I=0,NXP
          MATLX(I,0)=0.
          MATDX(I,0)=-1.
          MATUX(I,0)=1.
          VECX(I,0)=TH_BC_ZMIN_C1(N)*(25.d0*PR(2)/PR(1))*(C_sp_heat/L_heat)* &
                    (THX(I,1,J,1)+theta_0)*dthdz_mean(J,1)/m_FP !* &
!                   0.50d0*(INT_JACOB(0,J)+INT_JACOB(1,J))/CJOB_11(1,J,1)
        END DO
       ELSE
        DO I=0,NXP
          MATLX(I,0)=0.
          MATDX(I,0)=-1.
          MATUX(I,0)=1.
          VECX(I,0)=TH_BC_ZMIN_C1(N) 
        END DO
       ENDIF
      ELSE IF (TH_BC_ZMIN(N).eq.6) THEN
! C Periodic
       DO I=0,NXP          
          MATLX(I,0)=0. 
          MATDX(I,0)=1.
          MATUX(I,0)=0.                   
          VECX(I,0)=THX(I,NZ,J,N)
! 
! c          MATLX(I,1)=0. 
! c          MATDX(I,1)=1.
! c          MATUX(I,1)=0.                   
! c          VECX(I,1)=TH_BC_ZMIN_C1(N) 
        END DO
      !varying dirichlet: wave beam
       ELSE IF (TH_BC_ZMIN(N).EQ.9) THEN
        DO I=0,NXP
          MATLX(I,0)=0.
          MATDX(I,0)=1.
          MATUX(I,0)=0.
          VECX(I,0)=FACT_AMP*real(beam_th(j)*exp(-im*OMEGA0*TIME))

        END DO
      ELSE
        write(*,*) 'WARNING: TH_BC_LEFT is of unknown type'
      END IF

      RETURN 
      END

! C----*|--.---------.---------.---------.---------.---------.---------.-|----
      SUBROUTINE APPLY_BC_TH_RIGHT(N,J)
! C----*|--.---------.---------.---------.---------.---------.---------.-|--
use ntypes
use Domain
use Grid
! use Fft_var
use TIME_STEP_VAR
use run_variable
use ADI_var
use mg_vari, only : INIT_FLAG
      
implicit none
      INTEGER I,N,J

! C Right wall
      IF (TH_BC_ZMAX(N).EQ.0) THEN
! C Dirichlet
        DO I=0,NXP
          MATLX(I,NZ+1)=0.
          MATDX(I,NZ+1)=1.
          MATUX(I,NZ+1)=0.
          VECX(I,NZ+1)=0.

          MATLX(I,NZ)=0.
          MATDX(I,NZ)=1.
          MATUX(I,NZ)=0.
          VECX(I,NZ)=TH_BC_ZMAX_C1(N)
        END DO
      ELSE IF (TH_BC_ZMAX(N) .eq. 1) THEN
! C Neumann
! c        DO I=0,NXP
! c          MATLX(I,NZ)=-1.
! c          MATDX(I,NZ)=1.
! c          MATUX(I,NZ)=0.
! c          VECX(I,NZ)=TH_BC_ZMAX_C1(N)
! c        END DO
        DO I=0,NXP
          MATLX(I,NZ+1)=-1.
          MATDX(I,NZ+1)=1.
          MATUX(I,NZ+1)=0.
          VECX(I,NZ+1)=TH_BC_ZMAX_C1(N)
        END DO
      ELSE IF (TH_BC_ZMAX(N) .eq. 5) THEN
        DO I=0,NXP
          MATLX(I,NZ+1)=-2.0
          MATDX(I,NZ+1)=1.d0
          MATUX(I,NZ+1)=0.
          VECX(I,NZ+1) =0.
        END DO
      ELSE IF (TH_BC_ZMAX(N) .eq. 6) THEN
! C periodic condition      	
	DO I=0,NXP
          MATLX(I,NZ+1)=0.
          MATDX(I,NZ+1)=1.
          MATUX(I,NZ+1)=0.
          VECX(I,NZ+1)=THX(I,1,J,N)
	ENDDO  
      END IF

      RETURN
      END

! C----*|--.---------.---------.---------.---------.---------.---------.-|------
      SUBROUTINE APPLY_BC_1_LOWER(K)
! C----*|--.---------.---------.---------.---------.---------.---------.-|-----
use ntypes
use Domain
use Grid
! use Fft_var
use TIME_STEP_VAR
use run_variable
use ADI_var
use mg_vari, only : INIT_FLAG
      
implicit none
      INTEGER I,K

! C Bottom Wall:
      IF (U_BC_YMIN.EQ.0) THEN
! C Dirichlet
       If (F_TYPE.EQ.5)THEN
        DO I=0,NXP
          MATLY(I,0)=0.
          MATDY(I,0)=1.
          MATUY(I,0)=0.
          VECY(I,0)=0.

          MATLY(I,1)=0.
          MATDY(I,1)=1.
          MATUY(I,1)=0.
          VECY(I,1)=-Q_H0*(1.d0/cos(ANG_BETA)+GX(NX/2)*tan(ANG_BETA) &
                  *In_H0)*sin(OMEGA0*TIME)           
! c           VECY(I,1)=0.0
        END DO
       ELSE
        DO I=0,NXP
          MATLY(I,0)=0. 
          MATDY(I,0)=1.
          MATUY(I,0)=0.                   
          VECY(I,0)=0.

          MATLY(I,1)=0. 
          MATDY(I,1)=1.
          MATUY(I,1)=0.                   
          VECY(I,1)=U_BC_YMIN_C1 
        END DO
       ENDIF
      ELSE IF (U_BC_YMIN.eq.1) THEN
! C Neumann
        DO I=0,NXP
          MATLY(I,0)=0.
          MATDY(I,0)=-1.
          MATUY(I,0)=1.
          VECY(I,0)=DY(1)*U_BC_YMIN_C1
        END DO
      ELSE IF (U_BC_YMIN.eq.3) THEN
! C Wall model artificial slip boundary condition
! C Neumann
        DO I=0,NXP
          MATLY(I,1)=0.
          MATDY(I,1)=-1.
          MATUY(I,1)=1.
          VECY(I,1)=DY(2)*U_BC_LOWER_NWM(I,K)
        END DO
! C This gridpoint is not used here
        DO I=0,NXP
          MATLY(I,0)=0.
          MATDY(I,0)=1.
          MATUY(I,0)=0.
          VECY(I,0)=0.
        END DO
      ElSE IF (U_BC_YMIN .EQ. 6) THEN
! C Periodic
        DO I=0,NXP
          MATLY(I,0)=0. 
          MATDY(I,0)=1.
          MATUY(I,0)=0.                   
          VECY(I,0)=U1X(I,K,NY)
	ENDDO	
      ELSE
        write(*,*) 'WARNING: U_BC_LOWER is of unknown type'
      END IF

      RETURN 
      END



! C----*|--.---------.---------.---------.---------.---------.---------.-|----
      SUBROUTINE APPLY_BC_1_UPPER(K)
! C----*|--.---------.---------.---------.---------.---------.---------.-|--
use ntypes
use Domain
use Grid
! use Fft_var
use TIME_STEP_VAR
use run_variable
use ADI_var
use mg_vari, only : INIT_FLAG
      
implicit none
      INTEGER I,K

! C Top wall
      IF (U_BC_YMAX.EQ.0) THEN
! C Dirichlet
        DO I=0,NXP
          MATLY(I,NY+1)=0.
          MATDY(I,NY+1)=1.
          MATUY(I,NY+1)=0.
          VECY(I,NY+1)=0.

          MATLY(I,NY)=0.
          MATDY(I,NY)=1.
          MATUY(I,NY)=0.
          VECY(I,NY)=U_BC_YMAX_C1
        END DO
      ELSE IF (U_BC_YMAX.eq.1) THEN
! C Neumann
        DO I=0,NXP
          MATLY(I,NY+1)=-1.
          MATDY(I,NY+1)=1.
          MATUY(I,NY+1)=0.
          VECY(I,NY+1)=U_BC_YMAX_C1
        END DO
        
      ELSE IF (U_BC_YMAX.EQ.3) THEN
! C Use wall model artificial boundary condition
! C Neumann
        DO I=0,NXP
          MATLY(I,NY)=-1.
          MATDY(I,NY)=1.
          MATUY(I,NY)=0.
          VECY(I,NY)=DY(NY)*U_BC_UPPER_NWM(I,K)
        END DO
! C This gridpoint is not used here
        DO I=0,NXP
          MATLY(I,NY+1)=0.
          MATDY(I,NY+1)=1.
          MATUY(I,NY+1)=0.
          VECY(I,NY+1)=U1X(I,K,NY)
        END DO
      ELSE IF (U_BC_YMAX .eq. 6) THEN	
! c periodic 
       DO I=0,NXP
          MATLY(I,NY+1)=0.
          MATDY(I,NY+1)=1.
          MATUY(I,NY+1)=0.
          VECY(I,NY+1)=U1X(I,K,1)
       ENDDO	
       
      END IF

      RETURN
      END   

! C----*|--.---------.---------.---------.---------.---------.---------.-|------
      SUBROUTINE APPLY_BC_1_LEFT(J)
! C----*|--.---------.---------.---------.---------.---------.---------.-|-----
use ntypes
use Domain
use Grid
! use Fft_var
use TIME_STEP_VAR
use run_variable
use ADI_var
use mg_vari, only : INIT_FLAG
      
implicit none

      INTEGER I,J

! C Left Wall:
      IF (U_BC_ZMIN.EQ.0) THEN
! C Dirichlet
        DO I=0,NXP
          MATLX(I,0)=0. 
          MATDX(I,0)=1.
          MATUX(I,0)=0.                   
          VECX(I,0)=U_BC_ZMIN_C1
        IF (IC_TYPE .EQ. 5) THEN
          MATLX(I,1)=0.
          MATDX(I,1)=1.
          MATUX(I,1)=0.
          VECX(I,1)= U1_BAR(J)
        ELSE  
          MATLX(I,1)=0. 
          MATDX(I,1)=1.
          MATUX(I,1)=0.                   
          VECX(I,1)=U_BC_ZMIN_C1
         ENDIF 
        END DO
      ELSE IF (U_BC_ZMIN.eq.1) THEN
! C Neumann
       IF (IC_TYPE .EQ. 7) THEN
        DO I=0,NXP
          MATLX(I,0)=0.
          MATDX(I,0)=-1.
          MATUX(I,0)=1.
          VECX(I,0)=U_BC_ZMIN_C1
        END DO
       ELSE
        DO I=0,NXP
          MATLX(I,0)=0.
          MATDX(I,0)=-1.
          MATUX(I,0)=1.
          VECX(I,0)=DZ(1)*U_BC_ZMIN_C1
        END DO	  
       ENDIF	
      ELSE IF (U_BC_ZMIN.eq.3) THEN
! C Wall model artificial slip boundary condition
! C Neumann
        DO I=0,NXP
          MATLX(I,1)=0.
          MATDX(I,1)=-1.
          MATUX(I,1)=1.
          VECX(I,1)=DZ(2)*U_BC_LOWER_NWM(I,J)
        END DO
! C This gridpoint is not used here
        DO I=0,NXP
          MATLX(I,0)=0.
          MATDX(I,0)=1.
          MATUX(I,0)=0.
          VECX(I,0)=0.
        END DO
      ELSE IF (U_BC_ZMIN .eq. 6) THEN
! C Periodic
       DO I=0,NXP          
          MATLX(I,0)=0. 
          MATDX(I,0)=1.
          MATUX(I,0)=0.                   
          VECX(I,0)=U1X(I,NZ,J)
        END DO	
      ELSE
        write(*,*) 'WARNING: U_BC_LEFT is of unknown type'
      END IF

      RETURN 
      END

! C----*|--.---------.---------.---------.---------.---------.---------.-|----
      SUBROUTINE APPLY_BC_1_RIGHT(J)
! C----*|--.---------.---------.---------.---------.---------.---------.-|--
use ntypes
use Domain
use Grid
! use Fft_var
use TIME_STEP_VAR
use run_variable
use ADI_var
use mg_vari, only : INIT_FLAG
      
implicit none

      INTEGER I,J

! C Right wall
      IF (U_BC_ZMAX.EQ.0) THEN
! C Dirichlet
        DO I=0,NXP
          MATLX(I,NZ+1)=0.
          MATDX(I,NZ+1)=1.
          MATUX(I,NZ+1)=0.
          VECX(I,NZ+1)=0.

          MATLX(I,NZ)=0.
          MATDX(I,NZ)=1.
          MATUX(I,NZ)=0.
          VECX(I,NZ)=U_BC_ZMAX_C1
        END DO
      ELSE IF (U_BC_ZMAX.eq.1) THEN
! C Neumann
        IF (IC_TYPE == 7) THEN
	 DO I=0,NXP
           MATLX(I,NZ+1)=-1.
           MATDX(I,NZ+1)=1.
           MATUX(I,NZ+1)=0.
           VECX(I,NZ+1)=U_BC_ZMAX_C1
         END DO
	ELSE
         DO I=0,NXP
           MATLX(I,NZ)=-1.
           MATDX(I,NZ)=1.
           MATUX(I,NZ)=0.
           VECX(I,NZ)=DZ(NZ)*U_BC_ZMAX_C1
         END DO
	ENDIF
        DO I=0,NXP
          MATLX(I,NZ+1)=0.
          MATDX(I,NZ+1)=1.
          MATUX(I,NZ+1)=0.
          VECX(I,NZ+1)=0.
        END DO
      ELSE IF (U_BC_ZMAX.EQ.3) THEN
! C Use wall model artificial boundary condition
! C Neumann
        DO I=0,NXP
          MATLX(I,NZ)=-1.
          MATDX(I,NZ)=1.
          MATUX(I,NZ)=0.
          VECX(I,NZ)=DZ(NZ)*U_BC_UPPER_NWM(I,J)
        END DO
! C This gridpoint is not used here
        DO I=0,NXP
          MATLX(I,NZ+1)=0.
          MATDX(I,NZ+1)=1.
          MATUX(I,NZ+1)=0.
          VECX(I,NZ+1)=U1X(I,J,NZ)
        END DO
	
      ELSE IF (U_BC_ZMAX .eq. 6) THEN
! C periodic condition      	
	DO I=0,NXP
          MATLX(I,NZ+1)=0.
          MATDX(I,NZ+1)=1.
          MATUX(I,NZ+1)=0.
          VECX(I,NZ+1)=U1X(I,1,J)
	ENDDO	
	
      END IF

      RETURN
      END   

! C----*|--.---------.---------.---------.---------.---------.---------.-|------
      SUBROUTINE APPLY_BC_2_LOWER(K)
! C----*|--.---------.---------.---------.---------.---------.---------.-|-----
use ntypes
use Domain
use Grid
! use Fft_var
use TIME_STEP_VAR
use run_variable
use ADI_var
use mg_vari, only : INIT_FLAG
      
implicit none

      INTEGER I,K

! C Bottom Wall:
      IF (V_BC_YMIN.EQ.0) THEN
! C Dirichlet
       IF ( IC_TYPE .EQ. 7) THEN
        DO I=0,NXP
          MATLY(I,0)=0. 
          MATDY(I,0)=1.
          MATUY(I,0)=0.                   
          VECY(I,0)=0.

          MATLY(I,1)=0. 
          MATDY(I,1)=1.
          MATUY(I,1)=0.                   
          VECY(I,1)=V_BC_YMIN_C1 
        END DO 
       ELSE
        DO I=0,NXP
          MATLY(I,1)=0. 
          MATDY(I,1)=1.
          MATUY(I,1)=0.                   
          VECY(I,1)=V_BC_YMIN_C1

          MATLY(I,2)=0. 
          MATDY(I,2)=1.
          MATUY(I,2)=0.                   
          VECY(I,2)=V_BC_YMIN_C1 
        END DO
       ENDIF	
      ELSE IF (V_BC_YMIN.eq.1) THEN
! C Neumann
        DO I=0,NXP
          MATLY(I,1)=0.
          MATDY(I,1)=-1.
          MATUY(I,1)=1.
          VECY(I,1)=DYF(1)*V_BC_YMIN_C1
        END DO
      ELSE IF (V_BC_YMIN.eq.3) THEN
! C Wall model artificial slip boundary condition
! C Neumann
        DO I=0,NXP
          MATLY(I,1)=0.
          MATDY(I,1)=1.
          MATUY(I,1)=0.
          VECY(I,1)=0.
  
          MATLY(I,2)=0.
          MATDY(I,2)=1.
          MATUY(I,2)=0.
          VECY(I,2)=0.
        END DO
! C This gridpoint is not used here
        DO I=0,NXP
          MATLY(I,0)=0.
          MATDY(I,0)=1.
          MATUY(I,0)=0.
          VECY(I,0)=0.
        END DO
      ElSE IF (V_BC_YMIN .EQ. 6) THEN
! C Periodic
        DO I=0,NXP
          MATLY(I,0)=0. 
          MATDY(I,0)=1.
          MATUY(I,0)=0.                   
          VECY(I,0)=U2X(I,K,NY)
	ENDDO	
      ELSE
        write(*,*) 'WARNING: V_BC_LOWER is of unknown type'
      END IF
! 
! C The following is only a placeholder, this row is used for U1 and U3
! c      DO I=0,NXP
! c        MATLY(I,0) = 0.
! c        MATDY(I,0) = 1.
! c        MATUY(I,0) = 0.
! c        VECY(I,0) = 0.
! c      END DO



      RETURN 
      END



! C----*|--.---------.---------.---------.---------.---------.---------.-|----
      SUBROUTINE APPLY_BC_2_UPPER(K)
! C----*|--.---------.---------.---------.---------.---------.---------.-|--
use ntypes
use Domain
use Grid
! use Fft_var
use TIME_STEP_VAR
use run_variable
use ADI_var
use mg_vari, only : INIT_FLAG
      
implicit none

      INTEGER I,K

! C Top wall
      IF (V_BC_YMAX.EQ.0) THEN
! C Dirichlet
        DO I=0,NXP
          MATLY(I,NY+1)=0.
          MATDY(I,NY+1)=1.
          MATUY(I,NY+1)=0.
          VECY(I,NY+1)=V_BC_YMAX_C1

          MATLY(I,NY)=0.
          MATDY(I,NY)=1.
          MATUY(I,NY)=0.
          VECY(I,NY)=V_BC_YMAX_C1
        END DO
      ELSE IF (V_BC_YMAX.eq.1) THEN
! C Neumann
      IF (IC_TYPE ==7 ) THEN
       DO I=0,NXP
          MATLY(I,NY+1)=-1.
          MATDY(I,NY+1)=1.
          MATUY(I,NY+1)=0.
          VECY(I,NY+1)=V_BC_YMAX_C1
        END DO 
      ELSE	 
        DO I=0,NXP
          MATLY(I,NY+1)=-1.
          MATDY(I,NY+1)=1.
          MATUY(I,NY+1)=0.
          VECY(I,NY+1)=DYF(NY)*V_BC_YMAX_C1
        END DO
       ENDIF	
      ELSE IF (V_BC_YMAX.EQ.3) THEN
! C Use wall model artificial boundary condition
! C Neumann
        DO I=0,NXP
          MATLY(I,NY+1)=0.
          MATDY(I,NY+1)=1.
          MATUY(I,NY+1)=0.
          VECY(I,NY+1)=0.
        END DO
! C This gridpoint is not used here
        DO I=0,NXP
          MATLY(I,NY)=0.
          MATDY(I,NY)=1.
          MATUY(I,NY)=0.
          VECY(I,NY)=0.
        END DO
      ELSE IF (V_BC_YMAX .eq. 6) THEN	
! c periodic 
       DO I=0,NXP
          MATLY(I,NY+1)=0.
          MATDY(I,NY+1)=1.
          MATUY(I,NY+1)=0.
          VECY(I,NY+1)=U2X(I,K,1)
       ENDDO	
      END IF

      RETURN
      END   

! C----*|--.---------.---------.---------.---------.---------.---------.-|------
      SUBROUTINE APPLY_BC_2_LEFT(J)
! C----*|--.---------.---------.---------.---------.---------.---------.-|-----
use ntypes
use Domain
use Grid
! use Fft_var
use TIME_STEP_VAR
use run_variable
use ADI_var
use forcing
use mg_vari, only : INIT_FLAG
      
      implicit none

      INTEGER I,J

! C Left Wall:
      IF (V_BC_ZMIN.EQ.0) THEN
! C Dirichlet
        DO I=0,NXP
	IF (IC_TYPE .EQ. 7) THEN
          MATLX(I,0)=0. 
          MATDX(I,0)=1.
          MATUX(I,0)=0.                   
          VECX(I,0)=U2_BAR(1,j)
         
          MATLX(I,1)=0.
          MATDX(I,1)=1.
          MATUX(I,1)=0.
          VECX(I,1)=U2_BAR(1,j) 
         ELSE
	  MATLX(I,0)=0. 
          MATDX(I,0)=1.
          MATUX(I,0)=0.                   
          VECX(I,0)=0.
	 
          MATLX(I,1)=0. 
          MATDX(I,1)=1.
          MATUX(I,1)=0.                   
          VECX(I,1)=V_BC_ZMIN_C1
         ENDIF 
        END DO
      ELSE IF (V_BC_ZMIN.eq.1) THEN
! C Neumann
      IF (IC_TYPE .EQ. 7) THEN
        DO I=0,NXP
          MATLX(I,0)=0.
          MATDX(I,0)=-1.
          MATUX(I,0)=1.
          VECX(I,0)=V_BC_ZMIN_C1
        END DO
      ELSE
        DO I=0,NXP
          MATLX(I,0)=0.
          MATDX(I,0)=-1.
          MATUX(I,0)=1.
          VECX(I,0)=DZ(1)*V_BC_ZMIN_C1
        END DO 
      ENDIF 
      ELSE IF (V_BC_ZMIN.eq.3) THEN
! C Wall model artificial slip boundary condition
! C Neumann
! c        DO I=0,NXP
! c          MATLX(I,1)=0.
! c          MATDX(I,1)=-1.
! c          MATUX(I,1)=1.
! c          VECX(I,1)=0.
! c        END DO
! C This gridpoint is not used here
        DO I=0,NXP
          MATLX(I,0)=0.
          MATDX(I,0)=1.
          MATUX(I,0)=0.
          VECX(I,0)=DZ(1)*V_BC_ZMIN_C1
        END DO
      ELSE IF (V_BC_ZMIN .eq. 6) THEN
! C Periodic
       DO I=0,NXP          
          MATLX(I,0)=0. 
          MATDX(I,0)=1.
          MATUX(I,0)=0.                   
          VECX(I,0)=U2X(I,NZ,J)
        END DO
      ELSE IF (V_BC_ZMIN.EQ.7) THEN
         DO I=0,NXP
          MATLX(I,0)=0.
          MATDX(I,0)=(1. - dtc*C_int_le(J)*0.25*  &
          ( CJOB_11(0,J,2) + CJOB_11(1,J,2)        &
          + CJOB_11(0,J,1) + CJOB_11(0,J+1,1))    &
           /INT_JACOB(0,J))
          MATUX(I,0)=(dtc*C_int_le(J)*0.25*       &
         ( CJOB_11(0,J,2) + CJOB_11(1,J,2)         &
          + CJOB_11(0,J,1) + CJOB_11(0,J+1,1))    &
          /INT_JACOB(0,J))
          VECX(I,0) =U2X(I,0,J)
         END DO
      ELSE IF (V_BC_ZMIN .eq. 9) THEN     
!  Varying dirchlet: wave beam
        do I=0,NXP
          MATLX(I,0)=0.
          MATDX(I,0)=1.
          MATUX(I,0)=0.
          VECX(I,0)=FACT_AMP*real(v_prime(j)*exp(-im*OMEGA0*TIME))
        enddo  
      ELSE
        write(*,*) 'WARNING: V_BC_LEFT is of unknown type'
      END IF

      RETURN 
      END

! C----*|--.---------.---------.---------.---------.---------.---------.-|----
      SUBROUTINE APPLY_BC_2_RIGHT(J)
! C----*|--.---------.---------.---------.---------.---------.---------.-|--
use ntypes
use Domain
use Grid
! use Fft_var
use TIME_STEP_VAR
use run_variable
use ADI_var
use mg_vari, only : INIT_FLAG
      
implicit none

      INTEGER I,J

! C Right wall
      IF (V_BC_ZMAX.EQ.0) THEN
! C Dirichlet
        DO I=0,NXP
          MATLX(I,NZ+1)=0.
          MATDX(I,NZ+1)=1.
          MATUX(I,NZ+1)=0.
          VECX(I,NZ+1)=0.

          MATLX(I,NZ)=0.
          MATDX(I,NZ)=1.
          MATUX(I,NZ)=0.
          VECX(I,NZ)=V_BC_ZMAX_C1
        END DO
      ELSE IF (V_BC_ZMAX.eq.1) THEN
! C Neumann
       IF (IC_TYPE .EQ. 7) THEN
        DO I=0,NXP
          MATLX(I,NZ+1)=-1. 
          MATDX(I,NZ+1)=1.
          MATUX(I,NZ+1)=0.
          VECX(I,NZ+1)=V_BC_ZMAX_C1
        END DO
       ELSE
        DO I=0,NXP
          MATLX(I,NZ)=-1.
          MATDX(I,NZ)=1.
          MATUX(I,NZ)=0.
          VECX(I,NZ)=DZ(NZ)*V_BC_ZMAX_C1
        END DO
        DO I=0,NXP
          MATLX(I,NZ+1)=0.
          MATDX(I,NZ+1)=1.
          MATUX(I,NZ+1)=0.
          VECX(I,NZ+1)=0.
        END DO
       ENDIF
      ELSE IF (V_BC_ZMAX.EQ.4) THEN
! C  outlet Boundary condition
        DO I=0,NXP
          MATLX(I,NZ+1)=0.
          MATDX(I,NZ+1)=1.
          MATUX(I,NZ+1)=0.
          VECX(I,NZ+1)=U2X(I,NZ-1,J)+(U2X(I,NZ,J)-U2X(I,NZ-1,J)) &
                   *DZ(NZ+1)/DZ(NZ)

          MATLX(I,NZ)=-1.d0/DZ(NZ)
          MATDX(I,NZ)=1.d0/DZ(NZ)+1.d0/DZ(NZ+1)
          MATUX(I,NZ)=-1.d0/DZ(NZ+1)
          VECX(I,NZ) = 0.
        END DO

! c        DO I=0,NXP
! c          MATLX(I,NZ+1)=0.
! c          MATDX(I,NZ+1)=1.
! c          MATUX(I,NZ+1)=0.
! c          VECX(I,NZ+1)=U2(I,NZ-1,J)+(U2(I,NZ,J)-U2(I,NZ-1,J))
! c     +              *DZ(NZ+1)/DZ(NZ)
! 
! c          MATLX(I,NZ)=-(1.d0/DZ(NZ-1)+1.d0/DZ(NZ))
! c          MATDX(I,NZ)=1.d0/DZ(NZ)
! c          MATUX(I,NZ)=0.
! c          VECX(I,NZ) =0.
! c        END DO
      ELSE IF ( V_BC_ZMAX.EQ.5 ) THEN
       IF ( IC_TYPE == 7) THEN
        DO I=0,NXP
          MATLX(I,NZ+1)=-2.
          MATDX(I,NZ+1)=1.
          MATUX(I,NZ+1)=0.
          VECX(I,NZ+1) =0.
        END DO
       ELSE
        DO I=0,NXP
! c         MATLX(I,NZ)=-(1.d0/DZ(NZ-1)+1.d0/DZ(NZ))
! c          MATDX(I,NZ)=1.d0/DZ(NZ)
! c          MATUX(I,NZ)=0.
! c          VECX(I,NZ) =0.

          MATLX(I,NZ+1)=-(1.d0/DZ(NZ)+1.d0/DZ(NZ+1))
          MATDX(I,NZ+1)=1.d0/DZ(NZ+1)
          MATUX(I,NZ+1)=0.
          VECX(I,NZ+1) =0.
        END DO
       ENDIF
      ELSE IF (V_BC_ZMAX.EQ.3) THEN
! C Use wall model artificial boundary condition
! C Neumann
        DO I=0,NXP
          MATLX(I,NZ)=-1.
          MATDX(I,NZ)=1.
          MATUX(I,NZ)=0.
          VECX(I,NZ)=0.0
        END DO
! C This gridpoint is not used here
        DO I=0,NXP
          MATLX(I,NZ+1)=0.
          MATDX(I,NZ+1)=1.
          MATUX(I,NZ+1)=0.
          VECX(I,NZ+1)=U2X(I,NZ,J)
        END DO
      ELSE IF (V_BC_ZMAX .eq. 6) THEN
! C periodic condition      	
	DO I=0,NXP
          MATLX(I,NZ+1)=0.
          MATDX(I,NZ+1)=1.
          MATUX(I,NZ+1)=0.
          VECX(I,NZ+1)=U2X(I,1,J)
	ENDDO
      ELSE IF (V_BC_ZMAX.EQ.7) THEN
	 DO I=0,NXP
	  MATLX(I,NZ+1)=-(dtc*C_int(J)*0.25*             &
         ( CJOB_11(NZ+1,J,2) + CJOB_11(NZ+2,J,2)          &
          + CJOB_11(NZ+1,J,1) + CJOB_11(NZ+1,J+1,1))	  & 
     	  /INT_JACOB(NZ+1,J))
          MATDX(I,NZ+1)=(1. + dtc*C_int(J)*0.25*       &
          ( CJOB_11(NZ+1,J,2) + CJOB_11(NZ+2,J,2)       &
          + CJOB_11(NZ+1,J,1) + CJOB_11(NZ+1,J+1,1))	&  
     	  /INT_JACOB(NZ+1,J))
          MATUX(I,NZ+1)=0.
          VECX(I,NZ+1) =U2X(I,NZ+1,J)
         END DO			
      END IF

      RETURN
      END

! C----*|--.---------.---------.---------.---------.---------.---------.-|------
      SUBROUTINE APPLY_BC_3_LOWER(K)
! C----*|--.---------.---------.---------.---------.---------.---------.-|-----
use ntypes
use Domain
use Grid
! use Fft_var
use TIME_STEP_VAR
use run_variable
use ADI_var
use mg_vari, only : INIT_FLAG
      
implicit none

      INTEGER I,K

! C Bottom Wall:
      IF (W_BC_YMIN.EQ.0) THEN
! C Dirichlet
        DO I=0,NXP
          MATLY(I,0)=0. 
          MATDY(I,0)=1.
          MATUY(I,0)=0.                   
          VECY(I,0)=0.

          MATLY(I,1)=0. 
          MATDY(I,1)=1.
          MATUY(I,1)=0.                   
          VECY(I,1)=W_BC_YMIN_C1 
        END DO
      ELSE IF (W_BC_YMIN.eq.1) THEN
! C Neumann
        DO I=0,NXP
          MATLY(I,0)=0.
          MATDY(I,0)=-1.
          MATUY(I,0)=1.
          VECY(I,0)=DY(1)*W_BC_YMIN_C1
        END DO
      ELSE IF (W_BC_YMIN.eq.3) THEN
! C Wall model artificial slip boundary condition
! C Neumann
        DO I=0,NXP
          MATLY(I,1)=0.
          MATDY(I,1)=-1.
          MATUY(I,1)=1.
          VECY(I,1)=DY(2)*W_BC_LOWER_NWM(I,K)
        END DO
! C This gridpoint is not used here
        DO I=0,NXP
          MATLY(I,0)=0.
          MATDY(I,0)=1.
          MATUY(I,0)=0.
          VECY(I,0)=0.
        END DO
      ElSE IF (W_BC_YMIN .EQ. 6) THEN
! C Periodic
        DO I=0,NXP
          MATLY(I,0)=0. 
          MATDY(I,0)=1.
          MATUY(I,0)=0.                   
          VECY(I,0)=U3X(I,K,NY)
	ENDDO	
      ElSE IF (W_BC_YMIN .EQ. 7) THEN
! C Mixed boundary condition for w at the bottom wall 
       IF ( K .lE. 80 ) THEN
        DO I=0,NXP
          MATLY(I,0)=0.
          MATDY(I,0)=-1.
          MATUY(I,0)=1.
          VECY(I,0)=DY(1)*W_BC_YMIN_C1
        END DO
       ELSE
        DO I=0,NXP
          MATLY(I,0)=0.
          MATDY(I,0)=1.
          MATUY(I,0)=0.
          VECY(I,0)=0.

          MATLY(I,1)=0.
          MATDY(I,1)=1.
          MATUY(I,1)=0.
          VECY(I,1)=W_BC_YMIN_C1
        END DO
       ENDIF
 
      ELSE
        write(*,*) 'WARNING: W_BC_LOWER is of unknown type'
      END IF

      RETURN 
      END



! C----*|--.---------.---------.---------.---------.---------.---------.-|----
      SUBROUTINE APPLY_BC_3_UPPER(K)
! C----*|--.---------.---------.---------.---------.---------.---------.-|--
use ntypes
use Domain
use Grid
! use Fft_var
use TIME_STEP_VAR
use run_variable
use ADI_var
use mg_vari, only : INIT_FLAG
      
implicit none

      INTEGER I,K

! C Top wall
      IF (W_BC_YMAX.EQ.0) THEN
! C Dirichlet
        DO I=0,NXP
          MATLY(I,NY+1)=0.
          MATDY(I,NY+1)=1.
          MATUY(I,NY+1)=0.
          VECY(I,NY+1)=W_BC_YMAX_C1

          MATLY(I,NY)=0.
          MATDY(I,NY)=1.
          MATUY(I,NY)=0.
          VECY(I,NY)=W_BC_YMAX_C1
        END DO
      ELSE IF (W_BC_YMAX.eq.1) THEN
! C Neumann
       If (IC_TYPE == 7) THEN
        DO I=0,NXP
          MATLY(I,NY+1)=-1.0
          MATDY(I,NY+1)=1.0
          MATUY(I,NY+1)=0.
          VECY(I,NY+1)=W_BC_YMAX_C1
        END DO
       ELSE
        DO I=0,NXP
          MATLY(I,NY)=-1.
          MATDY(I,NY)=1.
          MATUY(I,NY)=0.
          VECY(I,NY)=DY(NY)*W_BC_YMAX_C1
        END DO
        DO I=0,NXP
          MATLY(I,NY+1)=0.
          MATDY(I,NY+1)=1.
          MATUY(I,NY+1)=0.
          VECY(I,NY+1)=0.
        END DO
       ENDIF	
      ELSE IF (W_BC_YMAX.EQ.3) THEN
! C Use wall model artificial boundary condition
! C Neumann
        DO I=0,NXP
          MATLY(I,NY)=-1.
          MATDY(I,NY)=1.
          MATUY(I,NY)=0.
          VECY(I,NY)=DY(NY)*W_BC_UPPER_NWM(I,K)
        END DO
! C This gridpoint is not used here
        DO I=0,NXP
          MATLY(I,NY+1)=0.
          MATDY(I,NY+1)=1.
          MATUY(I,NY+1)=0.
          VECY(I,NY+1)=U3X(I,K,NY)
        END DO
      ELSE IF (W_BC_YMAX .eq. 6) THEN	
! c periodic 
       DO I=0,NXP
          MATLY(I,NY+1)=0.
          MATDY(I,NY+1)=1.
          MATUY(I,NY+1)=0.
          VECY(I,NY+1)=U3X(I,K,1)
       ENDDO	
      END IF

      RETURN
      END   

! C----*|--.---------.---------.---------.---------.---------.---------.-|------
      SUBROUTINE APPLY_BC_3_LEFT(J)
! C----*|--.---------.---------.---------.---------.---------.---------.-|-----
use ntypes
use Domain
use Grid
! use Fft_var
use TIME_STEP_VAR
use run_variable
use ADI_var
use forcing
use mg_vari, only : INIT_FLAG
      
      implicit none

      INTEGER I,J,K

! C Bottom Wall:
      IF (W_BC_ZMIN.EQ.0) THEN
! C Dirichlet
       
        DO I=0,NXP
         IF (IC_TYPE .EQ. 5) THEN
          MATLX(I,1)=0.
          MATDX(I,1)=1.
          MATUX(I,1)=0.
          VECX(I,1)= U3_BAR(1,J)

          MATLX(I,2)=0.
          MATDX(I,2)=1.
          MATUX(I,2)=0.
          VECX(I,2)=U3_BAR(1,J)
        ELSEIF (IC_TYPE .EQ. 7) THEN

          IF(MELTING_MODEL)THEN
          MATLX(I,0)=0.
          MATDX(I,0)=1.
          MATUX(I,0)=0.
          VECX(I,0)=W_BC_ZMIN_C1*(C_sp_heat/L_heat)*(0.70*NU/Pr(1))*(1.0d0/1.0d0)* &
                    dthdz_mean(J,1)/abs(xpoint(2,1)-xpoint(1,1))
        
          MATLX(I,1)=0.
          MATDX(I,1)=1.
          MATUX(I,1)=0.
          VECX(I,1)=W_BC_ZMIN_C1*(C_sp_heat/L_heat)*(0.70*NU/Pr(1))*(1.0d0/1.0d0)* & 
                    dthdz_mean(J,1)/abs(xpoint(2,1)-xpoint(1,1))
          ELSE
	  MATLX(I,0)=0. 
          MATDX(I,0)=1.
          MATUX(I,0)=0.                   
          VECX(I,0)=U3_BAR(1,J)
	  
	  MATLX(I,1)=0. 
          MATDX(I,1)=1.
          MATUX(I,1)=0.                   
          VECX(I,1)= U3_BAR(1,J) !W_BC_ZMIN_C1
          ENDIF
	ELSE
          MATLX(I,1)=0. 
          MATDX(I,1)=1.
          MATUX(I,1)=0.                   
          VECX(I,1)=W_BC_ZMIN_C1

          MATLX(I,2)=0. 
          MATDX(I,2)=1.
          MATUX(I,2)=0.                   
          VECX(I,2)=W_BC_ZMIN_C1
         ENDIF 
        END DO
      ELSE IF (W_BC_ZMIN.eq.1) THEN
! C Neumann
      IF (IC_TYPE .EQ. 7) THEN 
	DO I=0,NXP
          MATLX(I,0)=0.
          MATDX(I,0)=-1.
          MATUX(I,0)=1.
          VECX(I,0)=W_BC_ZMIN_C1
        END DO
       ELSE
        DO I=0,NXP
          MATLX(I,1)=0.
          MATDX(I,1)=-1.
          MATUX(I,1)=1.
          VECX(I,1)=DZF(1)*W_BC_ZMIN_C1
        END DO
       ENDIF	
	
      ELSE IF (W_BC_ZMIN.eq.3) THEN
! C Wall model artificial slip boundary condition
! C Neumann
        DO I=0,NXP
          MATLX(I,1)=0.
          MATDX(I,1)=1.
          MATUX(I,1)=0.
          VECX(I,1)=0.
  
          MATLX(I,2)=0.
          MATDX(I,2)=1.
          MATUX(I,2)=0.
          VECX(I,2)=0.
        END DO
! C This gridpoint is not used here
        DO I=0,NXP
          MATLX(I,0)=0.
          MATDX(I,0)=1.
          MATUX(I,0)=0.
          VECX(I,0)=0.
        END DO
      ELSE IF (W_BC_ZMIN .eq. 6) THEN
! C Periodic
       DO I=0,NXP          
          MATLX(I,0)=0. 
          MATDX(I,0)=1.
          MATUX(I,0)=0.                   
          VECX(I,0)=U3X(I,NZ,J)
        END DO
      ELSE IF (W_BC_ZMIN.EQ.7) THEN
	 DO I=0,NXP
	  MATLX(I,0)=0.
          MATDX(I,0)=(1. - dtc*C_int_le(J)*0.25* &
          ( CJOB_11(0,J,2) + CJOB_11(1,J,2)      &
          + CJOB_11(0,J,1) + CJOB_11(0,J+1,1))	   &
     	  /INT_JACOB(0,J))
          MATUX(I,0)=(dtc*C_int_le(J)*0.25*       &
         ( CJOB_11(0,J,2) + CJOB_11(1,J,2)        &
          + CJOB_11(0,J,1) + CJOB_11(0,J+1,1))	  &
     	  /INT_JACOB(0,J))
          VECX(I,0) =U3X(I,0,J)
         END DO
      ELSE IF (W_BC_ZMIN .eq. 9) THEN
       do I=0,NXP
          MATLX(I,0)=0.
          MATDX(I,0)=1.
          MATUX(I,0)=0.
          VECX(I,0)=FACT_AMP*real(u_prime(j)*exp(-im*OMEGA0*TIME))
       end do
      ELSE
        write(*,*) 'WARNING: W_BC_LEFT is of unknown type'
      END IF

! C The following is only a placeholder, this row is used for U1 and U2
! c      DO I=0,NXP
! c        MATLX(I,0) = 0.
! c        MATDX(I,0) = 1.
! c        MATUX(I,0) = 0.
! c        VECX(I,0) = 0.
! c      END DO

      RETURN 
      END



! C----*|--.---------.---------.---------.---------.---------.---------.-|----
      SUBROUTINE APPLY_BC_3_RIGHT(J)
! C----*|--.---------.---------.---------.---------.---------.---------.-|--
use ntypes
use Domain
use Grid
! use Fft_var, only : NX2P
use TIME_STEP_VAR
use run_variable
use ADI_var
use mg_vari, only : INIT_FLAG
      
implicit none
      INTEGER I,K,J

! C Top wall
      IF (W_BC_ZMAX.EQ.0) THEN
! C Dirichlet
        DO I=0,NXP
          MATLX(I,NZ+1)=0.
          MATDX(I,NZ+1)=1.
          MATUX(I,NZ+1)=0.
          VECX(I,NZ+1)=W_BC_ZMAX_C1

          MATLX(I,NZ)=0.
          MATDX(I,NZ)=1.
          MATUX(I,NZ)=0.
          VECX(I,NZ)=W_BC_ZMAX_C1
        END DO
      ELSE IF (W_BC_ZMAX.eq.1) THEN
! C Neumann
      IF (IC_TYPE == 7) THEN
        DO I=0,NXP
          MATLX(I,NZ+1)=-1.
          MATDX(I,NZ+1)=1.
          MATUX(I,NZ+1)=0.
          VECX(I,NZ+1)=W_BC_ZMAX_C1
        END DO
      ELSE  
	DO I=0,NXP
          MATLX(I,NZ+1)=-1.
          MATDX(I,NZ+1)=1.
          MATUX(I,NZ+1)=0.
          VECX(I,NZ+1)=DZF(NZ)*W_BC_ZMAX_C1
        END DO
       ENDIF	
      ELSE IF (W_BC_ZMAX.EQ.4) THEN
! C  outlet Boundary condition
         DO I=0,NXP
           MATLX(I,NZ+1)=0.
           MATDX(I,NZ+1)=1.
           MATUX(I,NZ+1)=0.
           VECX(I,NZ+1)=U3X(I,NZ-1,J)+(U3X(I,NZ,J)-U3X(I,NZ-1,J)) &
                   *DZF(NZ)/DZF(NZ-1)

           MATLX(I,NZ)=-1.d0/DZF(NZ-1)
           MATDX(I,NZ)=1.d0/DZF(NZ)+1.d0/DZF(NZ-1)
           MATUX(I,NZ)=-1.d0/DZF(NZ)
           VECX(I,NZ) = 0.
          END DO 
	  
       ELSE IF (W_BC_ZMAX.EQ.5) THEN
        IF (IC_TYPE ==7 )THEN
	 DO I=0,NXP
	  MATLX(I,NZ+1)=-2.
          MATDX(I,NZ+1)=1.
          MATUX(I,NZ+1)=0.
          VECX(I,NZ+1) =0.
         END DO
	ELSE
         DO I=0,NXP
          MATLX(I,NZ+1)=-(1.d0/DZF(NZ-1)+1.d0/DZF(NZ))
          MATDX(I,NZ+1)=1.d0/DZF(NZ)
          MATUX(I,NZ+1)=0.
          VECX(I,NZ+1) =0.

! c          MATLX(I,NZ)=-(1.d0/DZF(NZ-1)+1.d0/DZF(NZ-2))
! c          MATDX(I,NZ)=1.d0/DZF(NZ-1)
! c          MATUX(I,NZ)=0.
! c          VECX(I,NZ) =0.
        END DO
	ENDIF
      
      ELSE IF (W_BC_ZMAX.EQ.3) THEN
! C Use wall model artificial boundary condition
! C Neumann
        DO I=0,NXP
          MATLX(I,NZ+1)=-1.
          MATDX(I,NZ+1)=1.
          MATUX(I,NZ+1)=0.
          VECX(I,NZ+1)=DZF(NZ)*W_BC_UPPER_NWM(I,K)
        END DO
! C This gridpoint is not used here
        DO I=0,NXP
          MATLX(I,NZ)=0.
          MATDX(I,NZ)=1.
          MATUX(I,NZ)=0.
          VECX(I,NZ)=0.
        END DO
      ELSE IF (W_BC_ZMAX .eq. 6) THEN
! C periodic condition      	
	DO I=0,NXP
          MATLX(I,NZ+1)=0.
          MATDX(I,NZ+1)=1.
          MATUX(I,NZ+1)=0.
          VECX(I,NZ+1)=U3X(I,1,J)
	ENDDO
      ELSE IF (W_BC_ZMAX.EQ.7) THEN
	 DO I=0,NXP
	  MATLX(I,NZ+1)=-(dtc*C_int(J)*0.25*  &
         ( CJOB_11(NZ+1,J,2) + CJOB_11(NZ+2,J,2) &
          + CJOB_11(NZ+1,J,1) + CJOB_11(NZ+1,J+1,1))&
     	  /INT_JACOB(NZ+1,J))
          MATDX(I,NZ+1)=(1. + dtc*C_int(J)*0.25* &
         ( CJOB_11(NZ+1,J,2) + CJOB_11(NZ+2,J,2) &
          + CJOB_11(NZ+1,J,1) + CJOB_11(NZ+1,J+1,1))	&  
     	  /INT_JACOB(NZ+1,J))
          MATUX(I,NZ+1)=0.
          VECX(I,NZ+1) =U3X(I,NZ+1,J)
         END DO		
      END IF

      RETURN
      END
! C----*|--.---------.---------.---------.---------.---------.---------.-|--
      SUBROUTINE APPLY_BC_VEL_LOWER
! C----*|--.---------.---------.---------.---------.---------.---------.-|--
! C This subroutine is called after initializing the flow
! C It sets the appropriate boundary conditions including ghost cell values
! C  on the velocity field in Fourier space
use ntypes
use Domain
use Grid
!use Fft_var, only : NX2P
use TIME_STEP_VAR
use run_variable
use ADI_var
use mg_vari, only : INIT_FLAG
      
implicit none
      INTEGER I,K      


      RETURN
      END

     
! C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE  APPLY_BC_VEL_UPPER
! C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
! C This subroutine is called after initializing the flow
! C It sets the appropriate boundary conditions including ghost cell values
! C  on the velocity field in Fourier space
use ntypes
use Domain
use Grid
!use Fft_var, only : NX2P
use TIME_STEP_VAR
use run_variable
use mg_vari, only : INIT_FLAG
      
implicit none

   
      RETURN
      END

! C----*|--.---------.---------.---------.---------.---------.---------.-|--
      SUBROUTINE APPLY_BC_VEL_LEFT
! C----*|--.---------.---------.---------.---------.---------.---------.-|--
! C This subroutine is called after initializing the flow
! C It sets the appropriate boundary conditions including ghost cell values
! C  on the velocity field in Fourier space
use ntypes
use Domain
use Grid
!use Fft_var, only : NX2P
use TIME_STEP_VAR
use run_variable
use mg_vari, only : INIT_FLAG
      
implicit none


      RETURN
      END

    
! C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE  APPLY_BC_VEL_RIGHT
! C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
! C This subroutine is called after initializing the flow
! C It sets the appropriate boundary conditions including ghost cell values
! C  on the velocity field in Fourier space
use ntypes
use Domain
use Grid
!use Fft_var, only : NX2P
use TIME_STEP_VAR
use run_variable
use mg_vari, only : INIT_FLAG
      
implicit none

      INTEGER I,J      
   
      RETURN
      END
