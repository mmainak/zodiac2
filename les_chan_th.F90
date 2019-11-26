      subroutine les_chan_th(n)
!C This subroutine models the subgridscale terms 
!c in the scalar advection equation for scalar number n
!C if the computation is to be treated as an LES not a DNS
!C This subroutine should be called when the velocity is in fourier space 
!C   in the periodic directions
!C S1 should contain |S| which was calculated in les_chan

use ntypes
use Domain
use Grid
use Fft_var
use TIME_STEP_VAR
use run_variable
use les_chan_var
use mpi_var, only : rank
implicit none

      CHARACTER*35 FNAME_TH
      integer i,j,k,l,m,ij,n

      real*8 S1_mean(1:NY)
 

! Variables for Dynamic Smagorinsky model:
      real*8 C_SMAG
      parameter (C_SMAG=0.13d0)

      real*8 alpha,beta 
      

! Here, alpha is the test/LES filter width ratio
      parameter (alpha=2.44950)
! beta is the LES/grid filter width ratio
      parameter (beta=1.d0)
      CHARACTER*31   file_name

      I = 1

! First, for all models, apply boundary conditions to the velocity field
! (fill ghost cells) to ensure accurate calculation of gradients
      IF (USE_MPI) THEN
        J1=JSTART
        J2=JEND 
      ELSE
        J1= 0
        J2= NY+1
      END IF

      if (LES_MODEL_TYPE_TH .EQ. 1) then
!     Constant Smagorinsky model

! First, compute the rate of strain tensor S_ij

      call compute_scalar_grad(n)

! Now, compute |S|*dTH/dx_i, storing in Sij
! First compute at GYF points 
      DO J=0,NY+1
        DO K=0,NZ+1
          DO I=0,NXP
            Sij(I,K,J,1)=Sij(I,K,J,1)
            Sij(I,K,J,3)=Sij(I,K,J,3)
          END DO
        END DO
      END DO
! Convert |S|*S_ij to Fourier space
      ij=1 

      S1X=Sij(:,:,:,ij)
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
      CSij(:,:,:,ij)=CS1X

!      CALL FFT_X_TO_FOURIER(Sij(0,0,0,ij),CSij(0,0,0,ij),0,NY+1,0,NZ+1)

! Sij(:,:,:,2) is added through an implicit eddy viscosity
      DO J=1,NY+1
        DO K=0,NZ+1
          DO I=0,NX2P
             CSij(I,K,J,2)=0.d0
             CSij(I,K,J,3)=0.d0
          END DO
        END DO
      END DO

! We now have |S|*dTH/dx_i stored in Sij(:,:,:,1..3) in Physical space

! Compute the filter lengthscale
! Absorb -2.d0*C_SMAG**2.d0 here for effienciency

      DO J=0,NY+1
       DO K=0,NZ+1
! At GY points:
        DELTA_Y(K,J)=-C_SMAG**2.d0*(beta*DX(1)*INT_JACOB(K,J))**(2.d0/3.d0)
!       DELTA_Y(J)=sqrt((beta*DX(1))**2.d0+(DY(J)*2.d0)**2.d0
!     &          +(beta*DZ(1))**2.d0)
       END DO
      END DO


! Get the eddy diffusivity at GY points
! NU_T = (C_S^2 * DELTA^2)*|S| With |S| interpolated to GY points
      DO J=JSTART,JEND
        DO K=ZSTART,ZEND
          DO I=0,NXP
!            KAPPA_T(I,K,J,N)=-1.d0*DELTA_Y(K,J)*S1(I,K,J)
             KAPPA_T(I,K,J,N) = NU_T(I,K,J)/0.85
          END DO
        END DO
      END DO

! Now, compute TAU, store in the corresponging Sij
      DO K=0,NZ
        DO I=0,NX2P
          DO J=0,NY
!            CSij(I,K,J,1)=DELTA_Y(K,J)*CSij(I,K,J,1)*S1(I,K,J)
             CSij(I,K,J,1) = -2.0d0*S2X(I,K,J)*KAPPA_T(I,K,J,N)*CSij(I,K,J,1)
          END DO
        END DO
      END DO

! tau_ij is now contained in CSij in Fourier space

! Convert the scalar to physical space
      CALL REAL_FOURIER_TRANS_TH (.false.)

      else if ((LES_MODEL_TYPE_TH.EQ.2).or.(LES_MODEL_TYPE_TH.eq.3)) then
! Here, use a dynamic smagorinsky model
! Note, there is no scale similar model for the scalar,
! so model type choice 2 and 3 are identical for the scalar equation

! Compute the filter width


      DO J=1,NY+1
       DO K=0,NZ+1
! At GY points:
       DELTA_Y(K,J)=beta*(beta*DX(1)*INT_JACOB(K,J))**(1.d0/3.d0)
!       DELTA_Y(J)=sqrt((beta*DX(1))**2.d0+(DY(J)*2.d0)**2.d0
!     &          +(beta*DZ(1))**2.d0)
       END DO
      END DO

! We need to calculate the components of C, the dynamic coefficient
! C_DYN_TH will be defined at GYF points

! Compute the scalar gradient, store in Sij(:,:,:,1..3)
      call compute_scalar_grad(n)

! Convert the scalar to physical space
      CALL REAL_FOURIER_TRANS_TH (.false.)


      
! Compute C_DYN_TH only every x # of timesteps
      if ( ((MOD(TIME_STEP,10).eq.0).AND.(RK_STEP.eq.1)).OR.FIRST_TIME) THEN

      if (rank .eq. 0) then
      write(6,*) '*************************'
      write(6,*) 'C_dyn_th is recalculating'
      write(6,*) '*************************'
      endif

      call allocate_les_tmp

! Store TH in Sij(:,:,:,4) and apply the test filter
      do j=0,NY+1
        do k=0,NZ+1
          do i=0,NXP
            Sij(i,k,j,4)=THX(i,k,j,n)
          end do
        end do
      end do
      S1X = Sij(:,:,:,4)
      call FILTER_VAR(2)
      Sij(:,:,:,4)=S1X

!      call les_filter_chan(Sij(0,0,0,4),0,NY+1,2)

! Zero C_DYN_TH
      do j=0,NY+1
       do k=0,NZ+1
        C_DYN_TH(k,j,N)=0.d0
        denominator_sum(k,j)=0.d0 
        numerator_sum(k,j)  =0.d0
       end do
      end do

! Do over all non-repeating components of the scalar gradient
      do ij=1,3

! Zero the numerator and denominator:
        do j=0,NY+1
          do k=0,NZ+1
            do i=0,NXP
              numerator(i,k,j)=0.d0
              denominator(i,k,j)=0.d0
            end do
          end do
        end do
        
! First, compute Mij

       do k=0,NZ+1
         do i=0,NXP
           do j=0,NY+1
              temp(i,k,j)=Sij(i,k,j,ij)
           end do
         end do
       end do
! Filter temp
      S1X = temp
      call FILTER_VAR(2)
      temp=S1X

!       call les_filter_chan(temp,0,NY+1,2)
! Multiply by |S| filtered        
       do j=0,NY+1
         do k=0,NZ+1
           do i=0,NXP
             temp(i,k,j)=temp(i,k,j)*(alpha*DELTA_Y(k,j))**2.d0 &
                                         *S_2BAR(i,k,j) 
           end do
         end do
       end do
! Get second term of Mij
         do i=0,NXP
           do k=0,NZ+1 
             do j=0,NY+1
               Mij(i,k,j)=DELTA_Y(k,j)**2.d0*S2X(i,k,j)*Sij(i,k,j,ij)
             end do
           end do
         end do
! Filter Mij
      S1X = Mij
      call FILTER_VAR(2)
      Mij=S1X

!       call les_filter_chan(Mij,0,NY+1,2)
 
! Add the second term of Mij stored in temp
       Mij=temp-Mij    
! Now, compute Lij and add Lij*Mij to the numerator
! temp=Ui*Uj:
       SELECT CASE (ij)
       CASE(1)
         do j=0,NY+1
           do k=0,NZ+1
             do i=0,NXP
               temp(i,k,j)=U1X(i,k,j)*THX(i,k,j,n)
             end do
           end do 
         end do
       CASE(2)
         do j=0,NY+1
           do k=0,NZ+1
             do i=0,NXP
               temp(i,k,j)=THX(i,k,j,n)*U2X(i,k,j)
             end do
           end do
         end do
       CASE(3)
         do j=0,NY+1
           do k=0,NZ+1
             do i=0,NXP
               temp(i,k,j)=U3X(i,k,j)*THX(i,k,j,n)
             end do
           end do 
         end do
       END SELECT
! Filter temp
      S1X = temp
      call FILTER_VAR(2)
      temp=S1X

!       call les_filter_chan(temp,0,NY+1,2)
! Add Lij*Mij to numerator
! Recall that Sij(:,:,:,4) holds TH_2BAR
       do j=0,NY+1
         do k=0,NZ+1
           do i=0,NXP 
             numerator(i,k,j)=Mij(i,k,j) &
                   *(temp(i,k,j)-Sij(i,k,j,4)*U_BAR_TIL(i,k,j,ij))
           end do
         end do
       end do

! Now, the denominator for this ij  
       do j=0,NY+1
         do k=0,NZ+1
           do i=0,NXP
             denominator(i,k,j)=Mij(i,k,j)*Mij(i,k,j)
           end do
         end do
       end do

! Get plane average of numerator and denominator, add to C
! Note, since both the numerator and denominator are averaged, there
! is not need to divide the sum by NX*NZ

      do j=0,NY+1
       do k=0,NZ+1
        denominator_sum(k,j)=denominator_sum(k,j)+   &
         SUM(denominator(0:min(NXP,NXP_L),k,j))
        numerator_sum(k,j)=numerator_sum(k,j)+SUM(numerator(0:min(NXP,NXP_L),k,j))
       end do
      end do

! End to ij
      end do
      
      CALL MPI_COMBINE_STATS(denominator_sum,NZ+2,NY+2)
      CALL MPI_COMBINE_STATS(numerator_sum,NZ+2,NY+2)

       do j=jstart,jend
        do k=zstart,zend
        if (denominator_sum(k,j) .ne. 0.) then
          C_DYN_TH(j,k,n)=-0.5d0*numerator_sum(k,j)/denominator_sum(k,j)
        else
          C_DYN_TH(j,k,n)=0.d0
        endif
       enddo
       enddo

! We are now done with the dynamic procedure to calculate C

! If C_DYN_TH < 0 at any level, set C_DYN_TH=0 for numerical stability
      do j=0,NY+1
       do k=0,NZ+1
        if (C_DYN_TH(k,j,n).lt.0) C_DYN_TH(k,j,n)=0.d0
       end do
      end do

      CALL MPI_BCAST_COMPLEX(C_DYN_TH(0,0,n), NZ+2, NY+2)

! End if compute C_DYN_TH
       call deallocate_les_tmp

       END IF

       

! Get the eddy diffusivity at GY points
! KAPPA_T = C_DYN_TH * (DELTA^2)*|S|
! Use exact second order interpolation to get C_DYN_TH and S1 interpolated to
! GY points
!c      DO J=J1,J2
      DO J=JSTART_TH(N),JEND_TH(N)
        DO K=ZSTART_TH(N),ZEND_TH(N)
          DO I=0,NXP
            KAPPA_T(I,K,J,N)=C_DYN_TH(k,j,n)*DELTA_Y(k,j)**2.d0*S2X(i,k,j)
          END DO
        END DO
      END DO

       

! At this point we have C_DYN_TH and dTH/dx_i (stored in Sij(:,:,:,1...3)
! Calculate lambda_i in physical space, stored in Sij(:,:,:,1..3)

      do ij=1,3
! Dynamic Smagorinsky model, no scale similar part
      do j=0,NY+1
       do k=0,NZ+1
        do i=0,NXP
            Sij(i,k,j,ij)=   &
             -C_DYN_TH(k,j,n)*DELTA_Y(k,j)**2.d0*S2X(i,k,j)*Sij(i,k,j,ij)
          end do 
        end do
       end do
       
       do k=0,NZ+1
          do i=0,NXP
           Sij(i,k,NY+1,ij)=Sij(i,k,NY,ij)
           Sij(i,k,0,ij)=Sij(i,k,1,ij)
          end do
         end do

        do j=0,NY+1
          do i=0,NXP
           Sij(i,NZ+1,j,ij)=Sij(i,NZ,j,ij)
           Sij(i,0,j,ij)=Sij(i,1,j,ij)
          end do
        end do

! Convert TAUij, now stored in Sij to Fourier space
       call FFT_X_TO_FOURIER(Sij(0,0,0,ij),CSij(0,0,0,ij),0,NY+1,0,NZ+1)

! End do ij
        end do  

! End if LES_MODEL_TYPE dynamic Smagorinsky or Dynamic model
      else
        pause 'Error, unsupported LES_MODEL_TYPE chosen'
      end if

! Now, add the subgrid scale forcing to CFi
! (This includes the subgrid scale stress as an explicit R-K term
! Include only CSij terms 1 and 3 since term 2 is accounted for
! as an implicit eddy diffusivity through KAPPA_T

!      DO J=J2+1,NY+1
!        DO K=0,TNKZ
!          DO I=0,NKX
!            DO ij =1,6
!             CSij(I,K,J,ij) = 0.d0
!            ENDDO
!           ENDDO
!         ENDDO
!       ENDDO

      DO J=JSTART_TH(N),JEND_TH(N) 
       DO K=ZSTART_TH(N),ZEND_TH(N)
        DO I=0,NX2P
          CFTHX(I,K,J,n)=CFTHX(I,K,J,n)       &
                     - INT_JACOB(K,J)*CIKXP(I)*CSij(I,K,J,1) 
         END DO
       END DO
      END DO

!      write(981,*) 'kappa_t_corner', SUM(KAPPA_T(0:NXM,NZ+1,NY+1,1))     
 
! Periodically, output mean quantities
      IF ((MOD(TIME_STEP,SAVE_STATS_INT).EQ.0).AND.(RK_STEP.EQ.1)) THEN

      do j=0,NY+1
        do k=0,NZ+1
         KAPPA_T_MEAN(k,j,n)=SUM(KAPPA_T(0:NXP,k,j,n))/dble(NX)
        end do
       end do

       CALL MPI_COMBINE_STATS(KAPPA_T_MEAN(0,0,n),NZ+2,NY+2)

      if (rank .eq. 0) then
      k = time_step/SAVE_STATS_INT

     
      

      file_name = 'plane_les_th/data_tec_'  &
             //CHAR(MOD(k,100000)/10000+48) &
             //CHAR(MOD(k,10000)/1000+48)   &
             //CHAR(MOD(k,1000)/100+48)     &
             //CHAR(MOD(k,100)/10+48)       &
             //CHAR(MOD(k,10)+48) //        &
             '.plt'

      call plot_les_th_tecplot(file_name)     
      endif

      ENDIF

      RETURN
      END


      subroutine compute_scalar_grad(n)
!C This subroutine computes dTH/dx_i for the filtered scalar field
!C The input velocity field should be in fourier space in the periodic
!C directions.
!C For use in the LES model in channel flow (2 periodic directions)
!C Store in Sij(:,:,:,1..3)
use ntypes
use Domain
use Grid
use Fft_var
use TIME_STEP_VAR
use run_variable
use les_chan_var

implicit none

      integer I,J,K,ij,n

      DO J=JSTART,JEND
        DO K=ZSTART,ZEND
          DO I=0,NX2P
            CSij(I,K,J,1)=CIKXP(I)*CTHX(I,K,J,n)
            CSij(I,K,J,2)=( 0.5*CJOB_12(K+1,J,2)*(CTHX(I,K,J,n) + CTHX(I,K+1,J,n))   &
                        - 0.5*CJOB_12(K,J,2)*(CTHX(I,K,J,n)   + CTHX(I,K-1,J,n))         &
                        + 0.5*CJOB_22(K,J+1,1)*(CTHX(I,K,J,n) + CTHX(I,K,J+1,n))         &
                        - 0.5*CJOB_22(K,J,1)*(CTHX(I,K,J,n)   + CTHX(I,K,J-1,n)) ) /     &
                             INT_JACOB(K,J)
            CSij(I,K,J,3)= ( 0.5*CJOB_11(K+1,J,2)*(CTHX(I,K,J,n) + CTHX(I,K+1,J,n)) &
                        - 0.5*CJOB_11(K,J,2)*(CTHX(I,K,J,n)   + CTHX(I,K-1,J,n))    &
                        + 0.5*CJOB_21(K,J+1,1)*(CTHX(I,K,J,n) + CTHX(I,K,J+1,n))    &
                        - 0.5*CJOB_21(K,J,1)*(CTHX(I,K,J,n)   + CTHX(I,K,J-1,n))  )/&
                             INT_JACOB(K,J)
          END DO
        END DO
      END DO

! Convert the scalar gradients to physical space
      do ij=1,3
      CS1X=CSij(:,:,:,ij)
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
      Sij(:,:,:,ij)=S1X

!        call FFT_X_TO_PHYSICAL(CSij(0,0,0,ij),Sij(0,0,0,ij),0,NY+1,0,NZ+1)
      end do

! We now have dTH/dx_i in Physical space
      RETURN
      END


! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! C *** Subroutine to write a .plt file (uses TecPlot tecio.a library)
! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
subroutine plot_les_th_tecplot(file_name)

use ntypes
use Domain
use run_variable, only : CU1X, CU2X, CU3X
use les_chan_var

      implicit none

      integer mm,nk,k

      integer  i,j,imax,jmax,kmax
      integer  debug,ier,itot
      integer  tecini,tecdat,teczne,tecnod,tecfil,tecend
      integer  visdouble,disdouble
      character*1 nulchar
      
      CHARACTER*31 file_name


      return
      end


