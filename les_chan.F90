      subroutine les_chan
! This subroutine models the terms owing to the subgrid scale stress
! if the computation is to be treated as an LES not a DNS
!C This subroutine should be called when the velocity is in fourier space 
!C in the periodic directions, on output, the velocity will be 
!C in physical space.
!C It is assumed that the test filter and the LES filter are performed
!C by the same operation
!C On output S1 should contain |S| which may be used again in les_chan_th
!C if for the subgrid scalar dissipation

use ntypes
use Domain
use Grid
use Fft_var
use TIME_STEP_VAR
use run_variable
use les_chan_var
use mpi_var, only : rank
implicit none

      integer i,j,k,l,m,ij
      CHARACTER*35 FNAME, FNAME_tke
      CHARACTER*28 file_name
      real*8 S1_mean(1:NY)

! Variables for Dynamic Smagorinsky model:
      real*8 C_SMAG
      parameter (C_SMAG=0.13d0)
       
      real*8 alpha,beta 
! Array to store the velocity index for each component of the strain rate tensor
      integer U_index1(6)
      integer U_index2(6)
      logical ADD_LES,LES_IMPLICIT

! Here, alpha is the test/LES filter width ratio
      parameter (alpha=2.44950)
! beta is the LES/grid filter width ratio
      parameter (beta=1.d0)


      ADD_LES = .TRUE.
      LES_IMPLICIT =.TRUE. 


! Set the velocity index for each component of the stress tensor
      U_index1(1)=1
      U_index2(1)=1

      U_index1(2)=2
      U_index2(2)=2 

      U_index1(3)=3
      U_index2(3)=3

      U_index1(4)=1
      U_index2(4)=2

      U_index1(5)=1
      U_index2(5)=3

      U_index1(6)=2
      U_index2(6)=3



! When using a Near-wall model, don't use LES at the wall
!      IF ((U_BC_YMIN.EQ.3).or.(W_BC_YMIN.EQ.3)) then
!C         J1=2
!C      ELSE
!C         J1=JSTART
!C      END IF

!C      IF ((U_BC_YMAX.EQ.3).or.(W_BC_YMAX.EQ.3)) then
!C         J2=NY-1
!C      ELSE
!C         J2=JEND
!C      END IF


! First, for all models, apply boundary conditions to the velocity field
! (fill ghost cells) to ensure accurate calculation of gradients
!C Apply Boundary conditions to velocity field
      IF (USE_MPI) THEN
        J1  = JSTART
        J2  = JEND
        J1i= 1
        J2e=NY+1
      ELSE
        J1  = JSTART
        J2  = JEND
        J1i= 1
        J2e=NY+1
      END IF


      if (LES_MODEL_TYPE.EQ.1) then
!     Constant Smagorinsky model

! First, compute the rate of strain tensor S_ij

      call compute_strain_chan

! Compute |S| at GYF points, store in S2X
! Interpolation to GYF points is easy since by definition
! GYF points are exactly midway between neighboring GY points
      DO J=JSTART,JEND
        DO K=ZSTART,ZEND
          DO I=0,NXP
            S2X(I,K,J)=SQRT(                                            &
                     2.d0*Sij(I,K,J,1)**2.d0                            & 
                    +4.d0*Sij(I,K,J,4)**2.d0                            &
                    +4.d0*Sij(I,K,J,5)**2.d0                            &
                    +2.d0*Sij(I,K,J,2)**2.d0                            & 
                    +4.d0*Sij(I,K,J,6)**2.d0                            &
                    +2.d0*Sij(I,K,J,3)**2.d0 )
          END DO
        END DO
      END DO
! Extend |S| to ghost cells
      DO K=0,NZ+1
        DO I=0,NXP
          S2X(I,K,0)=S2X(I,K,1)
          S2X(I,K,NY+1)=S2X(I,K,NY)
        END DO
      END DO


! Now, compute |S|*S_ij, storing in Sij
! First compute at GYF points 
      DO J=JSTART,JEND
        DO K=ZSTART,ZEND
          DO I=0,NXP
            Sij(I,K,J,1)=S2X(I,K,J)*Sij(I,K,J,1)
            Sij(I,K,J,5)=S2X(I,K,J)*Sij(I,K,J,5)
! CSij(:,:,:,2) is added through an implicit eddy viscosity
            Sij(I,K,J,2)=0.d0
            Sij(I,K,J,3)=S2X(I,K,J)*Sij(I,K,J,3)
          END DO
        END DO
      END DO
! Now, compute at GY points, interpolating |S|
      DO J=1,NY+1
        DO K=0,NZ+1
          DO I=0,NXP
! |S| interpolated to GY point 
            TEMP(I,K,J)=(S2X(I,K,J)*DYF(j-1)+S2X(I,K,J-1)*DYF(j)) &
                       /(2.d0*DY(j))
! The terms dU1/dy and dU3/dy in CSij(:,:,:,4) and CSij(:,:,:,6) respectively
! are subtracted out from Sij here since they are treated implicitly
! in eddy viscosity terms
            Sij(I,K,J,4)=TEMP(I,K,J)             &
             *(Sij(I,K,J,4)-0.5*(CU1X(I,K,J)-CU1X(I,K,J-1))/DY(j))
            Sij(I,K,J,6)=TEMP(I,K,J)             &
             *(Sij(I,K,J,6)-0.5*(CU3X(I,K,J)-CU3X(I,K,J-1))/DY(j))
          END DO
        END DO
      END DO



! We now have |S|*S_ij stored in Sij in Physical space

! Convert |S|*S_ij to Fourier space
      do ij=1,6
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
!        CALL FFT_X_TO_FOURIER(Sij(0,0,0,ij),CSij(0,0,0,ij),0,NY+1,0,NZ+1)
      end do
! Compute the filter lengthscale
! Absorb -2.d0*C_SMAG**2.d0 here for effienciency
      DO J=1,NY
! At GYF points:
! Constant Smagorinsky
!        DELTA_YF(J)=-2.d0*C_SMAG**2.d0
!     &     *(DX(1)*beta*DYF(J)*2.d0*DZ(1)*beta)**(2.d0/3.d0)
! Wall Damping
        DELTA_YF(J)=                                                 &
         -2.d0*(0.1d0*(1.d0-exp((-GYF(J)/(NU*25.d0))**3.d0)))**2.d0  &
                 *(DX(1)*beta*DYF(J)*2.d0*DZ(1)*beta)**(2.d0/3.d0)

      END DO
! Extend to ghost cells 
      DELTA_YF(0)=DELTA_YF(1)
      DELTA_YF(NY+1)=DELTA_YF(NY)      

      DO J=1,NY+1
! At GY points:
! Constant Smagorinsky
!        DELTA_Y(J)=-2.d0*C_SMAG**2.d0
!     &        *(DX(1)*beta*DY(J)*2.d0*DZ(1)*beta)**(2.d0/3.d0)
! Wall Damping
        DELTA_Y(k,j)=                                                 &
         -2.d0*(0.1d0*(1.d0-exp((-GY(J)/(NU*25.d0))**3.d0)))**2.d0  &
                 *(DX(1)*beta*DY(J)*2.d0*DZ(1)*beta)**(2.d0/3.d0)
      END DO

! Get the eddy viscosity at GY points
! NU_T = (C_S^2 * DELTA^2)*|S|

      DO J=1,NY+1
        DO K=0,NZ+1
          DO I=0,NXP
            NU_T(I,K,J)=-0.5d0*DELTA_Y(K,J)*TEMP(I,K,J)
          END DO
        END DO
      END DO

! Now, compute TAU, store in the corresponging Sij
      DO J=1,NY
        DO K=1,NZ
          DO I=0,NX2P
            CSij(I,K,J,1)=DELTA_Y(k,J)*CSij(I,K,J,1)
            CSij(I,K,J,5)=DELTA_Y(k,J)*CSij(I,K,J,5)
! CSij(:,:,:,2) is added through an implicit eddy viscosity
!            CSij(I,K,J,2)=DELTA_Y(K,J)*CSij(I,K,J,2)
            CSij(I,K,J,3)=DELTA_Y(K,J)*CSij(I,K,J,3)
          END DO
        END DO
      END DO
      DO J=1,NY+1 
        DO K=1,NZ+1
          DO I=0,NX2P
            CSij(I,K,J,4)=DELTA_Y(K,J)*CSij(I,K,J,4)
            CSij(I,K,J,6)=DELTA_Y(K,J)*CSij(I,K,J,6)
          END DO
        END DO
      END DO

! tau_ij is now contained in CSij in Fourier space

! Convert the velocity to physical space

      CALL REAL_FOURIER_TRANS_U1 (.false.)
      CALL REAL_FOURIER_TRANS_U2 (.false.)
      CALL REAL_FOURIER_TRANS_U3 (.false.)


!     CALL FFT_X_TO_PHYSICAL(CU1,U1,0,NY+1,0,NZ+1)
!     CALL FFT_X_TO_PHYSICAL(CU2,U2,0,NY+1,0,NZ+1)
!     CALL FFT_X_TO_PHYSICAL(CU3,U3,0,NY+1,0,NZ+1)


      else if ((LES_MODEL_TYPE.EQ.2).or.(LES_MODEL_TYPE.eq.3)) then
! Here, use a dynamic smagorinsky model with or without scale similar part

! Compute the filter width
      DO J=0,NY+1
! At GYF points:
        DELTA_YF(J)=(beta*DX(1)*DYF(J)*beta*DZ(1))**(1.d0/3.d0)
!        DELTA_YF(J)=sqrt((beta*DX(1))**2.d0+(DYF(J)*2.d0)**2.d0
!     &      +(beta*DZ(1))**2.d0)
      END DO


      DO J=0,NY+1
       DO K=0,NZ+1
! At cell centered  points:
        DELTA_Y(K,J)=(beta*DX(1)*INT_JACOB(K,J))**(1.d0/3.d0)
       END DO
      END DO

! We need to calculate the components of C, the dynamic coefficient
      
! Compute the rate of strain tensor, store in Sij
      call compute_strain_chan
      
 
! Compute |S| , store in S2X, bcz S2X will be used for ffts and filter variables
     DO J=JSTART,JEND
        DO K=ZSTART,ZEND
          DO I=0,NXP
            S2X(I,K,J)=SQRT(                                            &
                     2.d0*Sij(I,K,J,1)**2.d0                            &
                    +4.d0*Sij(I,K,J,4)**2.d0                            &
                    +4.d0*Sij(I,K,J,5)**2.d0                            &
                    +2.d0*Sij(I,K,J,2)**2.d0                            &
                    +4.d0*Sij(I,K,J,6)**2.d0                            &
                    +2.d0*Sij(I,K,J,3)**2.d0 )
          END DO
        END DO
      END DO
! Extend |S| to ghost cells
      DO K=0,NZ+1
        DO I=0,NXP
          S2X(I,K,0)   =S2X(I,K,1)
          S2X(I,K,NY+1)=S2X(I,K,NY)
        END DO
      END DO   


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!      Storing for tke les
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      DO ij=1,6
       DO J=0,NY+1
         DO K=0,NZ+1
          DO I=0,NXM
            St_rij(I,K,J,ij)=Sij(I,K,J,ij)
          ENDDO
         ENDDO
       ENDDO
      ENDDO      


! Need to calculate Mean Strain Rate

      do ij=1,6 
       do j=0,NY+1
        do k=0,NZ+1
         Sij_mean(k,j,ij) = SUM(Sij(0:min(NXP,NXP_L),k,j,ij))/dble(NX)
        enddo
       enddo
       CALL MPI_COMBINE_STATS(Sij_mean(0,0,ij),NZ+2,NY+2)  
      enddo



      DO K=0,NZ+1
      DO J=0,NY+1
       U1_BAR_les(k,j) = CU1X(0,k,j)
       U2_BAR_les(k,j) = CU2X(0,k,j)
       U3_BAR_les(k,j) = CU3X(0,k,j)
      ENDDO
      ENDDO

      CALL MPI_BCAST_REAL(U1_BAR_les,NZ+2,NY+2)
      CALL MPI_BCAST_REAL(U2_BAR_les,NZ+2,NY+2)
      CALL MPI_BCAST_REAL(U3_BAR_les,NZ+2,NY+2)     

! Convert Ui to physical space


      CALL REAL_FOURIER_TRANS_U1 (.false.)
      CALL REAL_FOURIER_TRANS_U2 (.false.)
      CALL REAL_FOURIER_TRANS_U3 (.false.)

!      CALL FFT_X_TO_PHYSICAL(CU1,U1,0,NY+1,0,NZ+1)
!      CALL FFT_X_TO_PHYSICAL(CU2,U2,0,NY+1,0,NZ+1)
!      CALL FFT_X_TO_PHYSICAL(CU3,U3,0,NY+1,0,NZ+1)
      
! Apply the filter to the LES velocity and save
      do j=0,NY+1
        do k=0,NZ+1
          do i=0,NXP
            U_BAR_TIL(i,k,j,1)=U1X(i,k,j)
            U_BAR_TIL(i,k,j,2)=U2X(i,k,j)
            U_BAR_TIL(i,k,j,3)=U3X(i,k,j)
          end do
        end do
      end do

      if (les_model_type.eq.3) then
! storing in the U_2BAR
       do j=0,NY+1
         do k=0,NZ+1
           do i=0,NXP
             do ij =1,3
               U_2BAR(i,k,j,ij) = U_BAR_TIL(i,k,j,ij)
             end do
           end do
         end do
       end do
    
       do i=1,3
! application of grid filter i.e. filter_type = 1
        S1X = U_2BAR(:,:,:,i)         
        call FILTER_VAR(1)
        U_2BAR(:,:,:,i)=S1X
!        call les_filter_chan(U_2BAR(0,0,0,i),0,NY+1,1)
       enddo
      endif
   
      
! Now, filter the velocity
      do i=1,3
! application of test filter i.e. filter_type = 2
        S1X = U_BAR_TIL(:,:,:,i)         
        call FILTER_VAR(2)
        U_BAR_TIL(:,:,:,i)=S1X

!        call les_filter_chan(U_BAR_TIL(0,0,0,i),0,NY+1,2)
      end do
      

! Compute C_DYN only every x # of timesteps
      if (((MOD(TIME_STEP,10).eq.0).AND.(RK_STEP.eq.1)) &
            .or.FIRST_TIME) THEN

      if (rank .eq. 0)then
      write(6,*)'##############################'
      write(6,*)'C_dyn is recalculating' 
      write(6,*)'##############################'
      endif

      call allocate_les_tmp
      
! Filter |S| and store in S_2BAR
      DO J=0,NY+1 
        DO K=0,NZ+1
          DO I=0,NXP
            S_2BAR(I,K,J)=S2X(I,K,J)
          END DO
        END DO
      END DO
      
! Test filtering operation filter type = 2
        S1X = S_2BAR
        call FILTER_VAR(2)
        S_2BAR=S1X

!      call les_filter_chan(S_2BAR,0,NY+1,2)

! Save a copy of the velocity which will be filtered twice more 

      if (les_model_type.eq.3) then 
! Do only if a scale similar part is needed
        do ij=1,3  
        do j=0,NY+1
          do k=0,NZ+1
            do i=0,NXP
              U_4BAR(i,k,j,ij)=U_BAR_TIL(i,k,j,ij)
            end do
          end do
        end do
        end do

        do i=1,3
! application of test and bar filter togather i.e. two filtering
! operation first filter type = 1, second filter type = 2
          S1X = U_4BAR(:,:,:,i)
          call FILTER_VAR(1)
          call FILTER_VAR(2)
          U_4BAR(:,:,:,i)=S1X
!          call les_filter_chan(U_4BAR(0,0,0,i),0,NY+1,1)
!          call les_filter_chan(U_4BAR(0,0,0,i),0,NY+1,2)
        end do
      end if


! Zero C_DYN
      do j=0,NY+1
       do k=0,NZ+1
        C_DYN(k,j)=0.d0
       end do
      end do

      

! The prep. work is now done, we are ready to start the algorithm to compute
! the dynamic model coefficient

      DO j =0,NY+1
       DO k =0,NZ+1
        denominator_sum(k,j) = 0.d0
        numerator_sum(k,j)   = 0.d0
       ENDDO
      ENDDO
 
! Do over all non-repeating components of the stress tensor
      do ij=1,6

! Here        ij=1 -> l=1,m=1
!             ij=2 -> l=2,m=2
!             ij=3 -> l=3,m=3
!             ij=4 -> l=1,m=2
!             ij=5 -> l=1,m=3
!             ij=6 -> l=2,m=3       
  
! Zero the numerator and denominator:
        do j=0,NY+1
          do k=0,NZ+1
            do i=0,NXP
              numerator(i,k,j)=0.d0
              denominator(i,k,j)=0.d0
            end do
          end do
        end do

! cross is used to multiply by two to include the contribution from the
! other symmetric term if we are dealing with a cross-term
        if (ij.le.3) then 
! We are computing a diagonal term
          cross=1.d0
        else
! We are computing a cross term
          cross=2.d0 
        end if 

! First, compute Mij

       do j=0,NY+1 
         do k=0,NZ+1
           do i=0,NXP
! Sij is defined at GYF points, no interpolation needed
             temp(i,k,j)=Sij(i,k,j,ij)
           end do
         end do
        end do
       
! Filter temp test filter operation filter type = 2
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

      
 
! Sij is used for Mij 
         do j=0,NY+1
           do i=0,NXP
             do k=0,NZ+1 
               Mij(i,k,j)=DELTA_Y(k,j)**2.d0*S2X(i,k,j)*Sij(i,k,j,ij)
             end do
           end do
         end do

! Filter Mij test filter operation filter type = 2
        S1X = Mij
        call FILTER_VAR(2)
        Mij=S1X 
!       call les_filter_chan(Mij,0,NY+1,2)
 
! Add the second term of Mij stored in temp
       do j=0,NY+1
         do k=0,NZ+1
           do i=0,NXP
             Mij(i,k,j)=temp(i,k,j)-Mij(i,k,j)    
           end do
         end do
       end do
! Now, compute Lij and add Lij*Mij to the numerator
! temp=Ui*Uj:
       SELECT CASE (ij)
       CASE(1)
         do j=0,NY+1
           do k=0,NZ+1
             do i=0,NXP
               temp(i,k,j)=U1X(i,k,j)*U1X(i,k,j)
             end do
           end do 
         end do
       CASE(2)
         do j=0,NY+1
           do k=0,NZ+1
             do i=0,NXP
               temp(i,k,j)=U2X(i,k,j)*U2X(i,k,j)
             end do
           end do
         end do
       CASE(3)
         do j=0,NY+1
           do k=0,NZ+1
             do i=0,NXP
               temp(i,k,j)=U3X(i,k,j)*U3X(i,k,j)
             end do
           end do 
         end do
       CASE(4) 
         do j=0,NY+1
           do k=0,NZ+1
             do i=0,NXP
               temp(i,k,j)=U1X(i,k,j)*U2X(i,k,j)
             end do
           end do
         end do
       CASE(5)
         do j=0,NY+1
           do k=0,NZ+1
             do i=0,NXP
               temp(i,k,j)=U1X(i,k,j)*U3X(i,k,j)
             end do
           end do 
         end do
       CASE(6) 
         do j=0,NY+1
           do k=0,NZ+1
             do i=0,NXP
               temp(i,k,j)=U2X(i,k,j)*U3X(i,k,j)
             end do
           end do
         end do
!end select
       END SELECT

       do k=0,NZ+1
           do i=0,NXP
             do j=0,NY+1
              temp_1(i,k,j) = temp(i,k,j)
             end do
           end do
       end do
! Filter temp: test filter operation i.e. filter type = 2
       S1X = temp
       call FILTER_VAR(2)
       temp=S1X
!       call les_filter_chan(temp,0,NY+1,2)
! Add Lij*Mij to numerator
       do j=0,NY+1
         do k=0,NZ+1
           do i=0,NXP 
             numerator(i,k,j)=Mij(i,k,j) &
             *(temp(i,k,j)               &
          -U_BAR_TIL(i,k,j,U_index1(ij))*U_BAR_TIL(i,k,j,U_index2(ij)))
           end do
         end do
       end do

       if (LES_MODEL_TYPE.eq.3) then
! If mixed model, include the scale similar part  
! Add Nij*Mij to the numerator piece-by-piece
! Term3 

! similarly grid filter operation is done on temp i.e. filter type = 1
         S1X = temp_1
         call FILTER_VAR(1)
!         call les_filter_chan(temp_1,0,NY+1,1)
! Filter temp_1 test filter operation filter type = 2
         call FILTER_VAR(2)
         temp_1=S1X
!         call les_filter_chan(temp_1,0,NY+1,2)
         do j=0,NY+1
           do k=0,NZ+1
             do i=0,NXP
              numerator(i,k,j)=numerator(i,k,j)+temp_1(i,k,j)*Mij(i,k,j)
             end do
           end do
         end do
! Term 4
         do j=0,NY+1
           do k=0,NZ+1
             do i=0,NXP
               temp(i,k,j)=U_2BAR(i,k,j,U_index1(ij)) &
                         *U_2BAR(i,k,j,U_index2(ij))
               temp_1(i,k,j) = U_BAR_TIL(i,k,j,U_index1(ij)) &
                         *U_BAR_TIL(i,k,j,U_index2(ij))
             end do
           end do
         end do
! Filter temp test filter operation filter type = 2
         S1X = temp
         call FILTER_VAR(2)
         temp=S1X 
!         call les_filter_chan(temp,0,NY+1,2)
         do j=0,NY+1
           do k=0,NZ+1
             do i=0,NXP
               numerator(i,k,j)=numerator(i,k,j)-temp(i,k,j)*Mij(i,k,j)
             end do
           end do
         end do 
! Terms 1 and 2
! Filter temp_1 grid filter operation filter type = 1
         S1X = temp_1
         call FILTER_VAR(1)
         call FILTER_VAR(2)
         temp_1=S1X
!         call les_filter_chan(temp_1,0,NY+1,1)
! Filter temp test filter operation filter type = 2
!         call les_filter_chan(temp_1,0,NY+1,2)
         do j=0,NY+1
           do k=0,NZ+1
             do i=0,NXP
              numerator(i,k,j)=numerator(i,k,j)                       &
             -(temp_1(i,k,j)                                          &
              -U_4BAR(i,k,j,U_index1(ij))*U_4BAR(i,k,j,U_index2(ij))) &
                *Mij(i,k,j)
             end do
           end do
         end do
       end if
! Now, the denominator for this ij
       do j=0,NY+1
         do k=0,NZ+1
           do i=0,NXP
             denominator(i,k,j)=Mij(i,k,j)*Mij(i,k,j)
           end do
         end do
       end do

        DO j=0,NY+1
         do k=0,NZ+1
          denominator_sum(k,j) = denominator_sum(k,j)   +          &
                     cross*SUM(denominator(0:min(NXP,NXP_L),k,j))
          numerator_sum(k,j)   = numerator_sum(k,j)     +          &
                     cross*SUM(numerator(0:min(NXP,NXP_L),k,j))
         ENDDO
        ENDDO 


  
! End to ij
      end do


      CALL MPI_COMBINE_STATS(denominator_sum,NZ+2,NY+2)
      CALL MPI_COMBINE_STATS(numerator_sum,NZ+2,NY+2)

! Get plane average of numerator and denominator, add to C
! Note, since both the numerator and denominator are averaged, there
! is not need to divide the sum by NX*NZ

       do j=jstart,jend
         do k=zstart,zend
          if (denominator_sum(k,j).ne.0.) then
            C_DYN(k,j) =-0.5d0* numerator_sum(k,j)/denominator_sum(k,j)
          else 
            C_DYN(k,j)=0.d0  
          end if
         end do
       end do




! We are now done with the dynamic procedure to calculate C
      

! If C_DYN < 0 at any level, set C_DYN=0 for numerical stability
      do j=0,NY+1
      do k=0,NZ+1
        if (C_DYN(k,j).lt.0) C_DYN(k,j)=0.d0
      end do
      end do

! At this point we have C_DYN and Sij
      
      CALL MPI_BCAST_COMPLEX(C_DYN, NZ+2, NY+2)
 
      
! End if compute C_DYN
       
       call deallocate_les_tmp

       END IF

      

! Get the eddy viscosity at cell centered points
! NU_T = C_DYN * (DELTA^2)*|S|

!      DO J=J1,J2
       DO J=JSTART,JEND
        DO K=ZSTART,ZEND
          DO I=0,NXP
            NU_T(I,K,J)=C_DYN(k,j)*DELTA_Y(k,j)**2.d0*S2X(i,k,j)
          END DO
        END DO
      END DO
      
       

! Calculate TAUij in physical space, stored in Sij

      do ij=1,6
      if (LES_MODEL_TYPE.eq.2) then
! Dynamic Smagorinsky model, no scale similar part
      if ( ij.eq.1 ) then
! Here, Sij is defined at GYF points 
      do j=1,NY+1
        do k=1,NZ+1
          do i=0,NXP
            Sij(i,k,j,ij)= &
             -2.d0*C_DYN(k,j)*DELTA_Y(k,j)**2.d0*S2X(i,k,j)*Sij(i,k,j,ij)
          end do 
        end do
      end do
      else if (ij.eq.2) then
! Sij(:,:,:,2) = du2/dy, but this term will be added implicitly through
! an eddy viscosity, so set it equal to zero here
        do j=1,NY+1
          do k=1,NZ+1
            do i=0,NXP
              Sij(i,k,j,ij)=0.d0
            end do
          end do
        end do
       else if (ij.eq.3) then
! Sij(:,:,:,3) = du3/dz, but this term will be added implicitly through
! an eddy viscosity, so set it equal to zero here
        do j=0,NY+1
          do k=1,NZ+1
            do i=0,NXP
              Sij(i,k,j,ij)=0.d0
            end do
          end do
        end do
       else if (ij.eq.4) then 
! Here, Sij is defined at GY points, interpolate C_DYN, etc 
! Use exact second order interpolation
! Sij(:,:,:,4)=0.5*(dU1/dy + dU2/dx)
! But dU1/dy will be accounted for as an implicit eddy viscosity term,
! So, subtract if off from Sij here
       
      do j=JSTART,JEND
        do k=ZSTART,ZEND
          do i=0,NXP
            temp_1(i,k,j)=-2.d0                                                     &
             *C_DYN(k,j)*DELTA_Y(k,j)**2.d0*S2X(i,k,j)*Sij(i,k,j,ij)

            Sij(i,k,j,ij)=-2.d0                                                     &
             *C_DYN(k,j)*DELTA_Y(k,j)**2.d0*S2X(i,k,j)*(Sij(i,k,j,ij)-              &
                        0.5* ( 0.5*CJOB_12(K+1,J,2)*(U1X(I,K,J)   + U1X(I,K+1,J))     &
                        - 0.5*CJOB_12(K,J,2)*(U1X(I,K,J)   + U1X(I,K-1,J))            &
                        + 0.5*CJOB_22(K,J+1,1)*(U1X(I,K,J) + U1X(I,K,J+1))            &
                        - 0.5*CJOB_22(K,J,1)*(U1X(I,K,J)   + U1X(I,K,J-1)))/          &
                        INT_JACOB(K,J)  )
          end do
        end do
      end do
      
 
      else if (ij.eq.5) then
! Here, Sij is defined at GY points, interpolate C_DYN, etc 
! Use exact second order interpolation
! Sij(:,:,:,5)=0.5*(dU1/dz + dU3/dx)
! But dU1/dz will be accounted for as an implicit eddy viscosity term,
! So, subtract if off from Sij here
      do j=JSTART,JEND
        do k=ZSTART,ZEND
          do i=0,NXP
            temp_2(i,k,j)=-2.d0                                    &
             *C_DYN(k,j)*DELTA_Y(k,j)**2.d0*S2X(i,k,j)*Sij(i,k,j,ij)

            Sij(i,k,j,ij)=-2.d0                                    &
             *C_DYN(k,j)*DELTA_Y(k,j)**2.d0*S2X(i,k,j)*(Sij(i,k,j,ij)-0.5d0*      &
                        ( 0.5*CJOB_11(K+1,J,2)*(U1X(I,K,J)+ U1X(I,K+1,J))          &
                        - 0.5*CJOB_11(K,J,2)*(U1X(I,K,J)   + U1X(I,K-1,J))         &
                        + 0.5*CJOB_21(K,J+1,1)*(U1X(I,K,J) + U1X(I,K,J+1))         &
                        - 0.5*CJOB_21(K,J,1)*(U1X(I,K,J)   + U1X(I,K,J-1)) )/      & 
                       INT_JACOB(K,J)   )                         
          end do
        end do
      end do
      else if (ij.eq.6) then     
! Here, Sij is defined at GY points, interpolate C_DYN, etc 
! Use exact second order interpolation
! Sij(:,:,:,6)=0.5*(dU3/dy + dU2/dz)
! But dU3/dy as well as  dU2/dz will be accounted for as an implicit eddy viscosity term,
! So, subtract if off from Sij here
      do j=JSTART,JEND
        do k=ZSTART,ZEND
          do i=0,NXP
             temp_3(i,k,j)=-2.d0                                                  &
             *C_DYN(k,j) *DELTA_Y(k,j)**2.d0*S2X(i,k,j)*(Sij(i,k,j,ij)-           &
                       0.5* ( 0.5*CJOB_12(K+1,J,2)*(U3X(I,K,J)   + U3X(I,K+1,J))  &
                        - 0.5*CJOB_12(K,J,2)*(U3X(I,K,J)   + U3X(I,K-1,J))        &
                        + 0.5*CJOB_22(K,J+1,1)*(U3X(I,K,J) + U3X(I,K,J+1))        &
                        - 0.5*CJOB_22(K,J,1)*(U3X(I,K,J)   + U3X(I,K,J-1)))/ &
                        INT_JACOB(K,J) )

            Sij(i,k,j,ij)=-2.d0                                    &
             *C_DYN(k,j)*DELTA_Y(k,j)**2.d0*S2X(i,k,j)*(Sij(i,k,j,ij)-           &
                       0.5*( 0.5*CJOB_11(K+1,J,2)*(U2X(I,K,J)+ U2X(I,K+1,J))       &
                        - 0.5*CJOB_11(K,J,2)*(U2X(I,K,J)   + U2X(I,K-1,J))         &
                        + 0.5*CJOB_21(K,J+1,1)*(U2X(I,K,J) + U2X(I,K,J+1))         &
                        - 0.5*CJOB_21(K,J,1)*(U2X(I,K,J)   + U2X(I,K,J-1)))/       &
                       INT_JACOB(K,J)  )
          end do
        end do
      end do 
       
       do k=0,NZ+1
        do i=0,NXP
           temp_3(i,k,NY+1) = temp_3(i,k,NY)
           temp_3(i,k,0)    = temp_3(i,k,1)
         end do
       end do

        do j=0,NY+1
          do i=0,NXP
           temp_3(i,NZ+1,j) = temp_3(i,NZ,j)
           temp_3(i,0,j)    = temp_3(i,1,j)
          end do
         end do

! End if ij
      end if

!    EXTRAPOLATING VALUES 

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

       

      else if (LES_MODEL_TYPE.eq.3) then
! Model type = 3, dynamic mixed model with scale similar part
! Always define temp at GYF points to match U_2BAR
! temp=Ui*Uj:
       SELECT CASE (ij)
       CASE(1)
         do j=0,NY+1
           do k=0,NZ+1
             do i=0,NXP 
               temp(i,k,j)=U1X(i,k,j)*U1X(i,k,j)
             end do
           end do
         end do
       CASE(2)
         do j=0,NY+1
           do k=0,NZ+1
             do i=0,NXP 
               temp(i,k,j)=U2X(i,k,j)*U2X(i,k,j)
             end do
           end do
         end do
       CASE(3)
         do j=0,NY+1
           do k=0,NZ+1
             do i=0,NXP 
               temp(i,k,j)=U3X(i,k,j)*U3X(i,k,j)
             end do
           end do
         end do
       CASE(4) 
         do j=0,NY+1
           do k=0,NZ+1
             do i=0,NXP 
               temp(i,k,j)=U2X(i,k,j)*U1X(i,k,j)
             end do
           end do
         end do
       CASE(5)
         do j=0,NY+1
           do k=0,NZ+1
             do i=0,NXP 
               temp(i,k,j)=U1X(i,k,j)*U3X(i,k,j)
             end do
           end do
         end do
       CASE(6) 
         do j=0,NY+1
           do k=0,NZ+1
             do i=0,NXP 
               temp(i,k,j)=U2X(i,k,j)*U3X(i,k,j)
             end do
           end do
         end do
       END SELECT

! Filter temp grid filter operation filter type = 1
       S1X = temp
       call FILTER_VAR(1)
       temp=S1X

!       call les_filter_chan(temp,0,NY+1,1)


      IF(LES_IMPLICIT) THEN


      if ( ij.eq.1 ) then
! Here, Sij is defined at GYF points 

      do j=0,NY+1
        do k=1,NZ+1
          do i=0,NXP
            Sij(i,k,j,ij)=temp(i,k,j)   &
              -U_2BAR(i,k,j,U_index1(ij))*U_2BAR(i,k,j,U_index2(ij)) &
              -2.d0*C_DYN(k,j)*DELTA_Y(k,j)**2.d0*S2X(i,k,j)*Sij(i,k,j,ij)
          end do
        end do 
      end do
       else if ( (ij.eq.2) .or. (ij.eq.3) ) then
! Sij(:,:,:,2) = du2/dy, 
! Sij(:,:,:,3) = du3/dz, but this term will be added implicitly through
! an eddy viscosity, so set the Smagorinsky pert equal to zero here
        do j=0,NY+1
          do k=1,NZ+1
            do i=0,NXP
              Sij(i,k,j,ij)=temp(i,k,j)                              &
             -U_2BAR(i,k,j,U_index1(ij))*U_2BAR(i,k,j,U_index2(ij))
            end do
          end do
        end do
        write(6,*) 'I am here 6666' 
      else if (ij.eq.4) then
! Sij(:,:,:,4)=0.5*(dU1/dy + dU2/dx)
! But the dU1/dy term will be accounted for as an implicit eddy viscosity
! so subtract this term from Sij in the Smagorinsky part
      do j=1,NY+1
        do k=1,NZ+1
          do i=0,NXP
            temp_1(i,k,j)=                                                        &
               temp(i,k,j) -U_2BAR(i,k,j,U_index1(ij))*U_2BAR(i,k,j,U_index2(ij)) &
                    -2.d0*C_DYN(k,j)*DELTA_Y(k,j)**2.d0*S2X(i,k,j)*Sij(i,k,j,ij)

            Sij(i,k,j,ij)=                                                        &
              temp(i,k,j) -U_2BAR(i,k,j,U_index1(ij))*U_2BAR(i,k,j,U_index2(ij))  &
           -2.d0*C_DYN(k,j)*DELTA_Y(k,j)**2.d0*S2X(i,k,j)*(Sij(i,k,j,ij)-0.5d0*   &
                         ( 0.5*CJOB_12(K+1,J,2)*(U1X(I,K,J)   + U1X(I,K+1,J))     &
                        - 0.5*CJOB_12(K,J,2)*(U1X(I,K,J)   + U1X(I,K-1,J))        &
                        + 0.5*CJOB_22(K,J+1,1)*(U1X(I,K,J) + U1X(I,K,J+1))        &
                        - 0.5*CJOB_22(K,J,1)*(U1X(I,K,J)   + U1X(I,K,J-1)))/      &
                        INT_JACOB(K,J)   )
          end do
        end do
      end do 
      else if (ij.eq.5) then
! Here, Sij is defined at GY points, interpolate C_DYN, etc 
! Use exact second order interpolation
! Sij(:,:,:,5)=0.5*(dU1/dz+dU3/dx) 
! But the dU1/dz term will be accounted for as an implicit eddy viscosity
! so subtract this term from Sij in the Smagorinsky part
      do j=1,NY+1
        do k=1,NZ+1
          do i=0,NXP
            temp_2(i,k,j)=                                                       &
               temp(i,k,j)-U_2BAR(i,k,j,U_index1(ij))*U_2BAR(i,k,j,U_index2(ij)) &
           -2.d0*C_DYN(k,j)*DELTA_Y(k,j)**2.d0*S2X(i,k,j)*Sij(i,k,j,ij)
 
            Sij(i,k,j,ij)=                                                       &
              temp(i,k,j)-U_2BAR(i,k,j,U_index1(ij))*U_2BAR(i,k,j,U_index2(ij))  &
           -2.d0*C_DYN(k,j)*DELTA_Y(k,j)**2.d0*S2X(i,k,j)*(Sij(i,k,j,ij)-0.5d0*  &
                        ( 0.5*CJOB_11(K+1,J,2)*(U1X(I,K,J)+ U1X(I,K+1,J))          &
                        - 0.5*CJOB_11(K,J,2)*(U1X(I,K,J)   + U1X(I,K-1,J))         &
                        + 0.5*CJOB_21(K,J+1,1)*(U1X(I,K,J) + U1X(I,K,J+1))         &
                        - 0.5*CJOB_21(K,J,1)*(U1X(I,K,J)   + U1X(I,K,J-1)) )/      &
                       INT_JACOB(K,J) )
          end do
        end do
      end do
      else if (ij.eq.6) then
! Here, Sij is defined at GY points, interpolate C_DYN, etc 
! Use exact second order interpolation
! Sij(:,:,:,6)=0.5*(dU3/dy + dU2/dz)
! But both  dU3/dy and dU2/dz  term will be accounted for as an implicit eddy viscosity
! so subtract this term from Sij in the Smagorinsky part
      do j=1,NY+1
        do k=1,NZ+1
          do i=0,NXP
            temp_3(i,k,j)=                                                         &
               temp(i,k,j) -U_2BAR(i,k,j,U_index1(ij))*U_2BAR(i,k,j,U_index2(ij))  &
           -2.d0*C_DYN(k,j)*DELTA_Y(k,j)**2.d0*S2X(i,k,j)*(Sij(i,k,j,ij)-          &
                       0.5* ( 0.5*CJOB_12(K+1,J,2)*(U3X(I,K,J)   + U3X(I,K+1,J))   &
                        - 0.5*CJOB_12(K,J,2)*(U3X(I,K,J)   + U3X(I,K-1,J))         &
                        + 0.5*CJOB_22(K,J+1,1)*(U3X(I,K,J) + U3X(I,K,J+1))         &
                        - 0.5*CJOB_22(K,J,1)*(U3X(I,K,J)   + U3X(I,K,J-1)))/       &
                        INT_JACOB(K,J) )

            Sij(i,k,j,ij)=                                                         &
              temp(i,k,j)  -U_2BAR(i,k,j,U_index1(ij))*U_2BAR(i,k,j,U_index2(ij))  &
           -2.d0*C_DYN(k,j)*DELTA_Y(k,j)**2.d0*S2X(i,k,j)*(Sij(i,k,j,ij)-          &
                        0.5*( 0.5*CJOB_11(K+1,J,2)*(U2X(I,K,J)+ U2X(I,K+1,J))      &
                        - 0.5*CJOB_11(K,J,2)*(U2X(I,K,J)   + U2X(I,K-1,J))         &
                        + 0.5*CJOB_21(K,J+1,1)*(U2X(I,K,J) + U2X(I,K,J+1))         &
                        - 0.5*CJOB_21(K,J,1)*(U2X(I,K,J)   + U2X(I,K,J-1)))/       &
                       INT_JACOB(K,J)  ) 
          end do
        end do
      end do 

! End if ij
       end if
        write(6,*)'I am here impli'

       ELSE
         do j=JSTART,JEND
         do k=ZSTART,ZEND
          do i=0,NXP
            Sij(i,k,j,ij)=temp(i,k,j)   &
              -U_2BAR(i,k,j,U_index1(ij))*U_2BAR(i,k,j,U_index2(ij)) &
              -2.d0*C_DYN(k,j)*DELTA_Y(k,j)**2.d0*S2X(i,k,j)*Sij(i,k,j,ij)
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

       ENDIF 

! End if Mixed Model
       end if

! Convert TAUij, now stored in Sij to Fourier space
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

!       call FFT_X_TO_FOURIER(Sij(0,0,0,ij),CSij(0,0,0,ij),0,NY+1,0,NZ+1)
! End do ij
       end do 
! Convert TAUij, now stored in Sij to Fourier space (for only u momentum equation)

      S1X=temp_1
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
      Ctemp_1=CS1X

      S1X=temp_2
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
      Ctemp_2=CS1X

      S1X=temp_3
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
      Ctemp_3=CS1X

!       call FFT_X_TO_FOURIER(temp_1,ctemp_1,0,NY+1,0,NZ+1)
!       call FFT_X_TO_FOURIER(temp_2,ctemp_2,0,NY+1,0,NZ+1) 
!       call FFT_X_TO_FOURIER(s2,cs2,0,NY+1,0,NZ+1)
! End if LES_MODEL_TYPE dynamic Smagorinsky or Dynamic mixed model
      else
        pause 'Error, unsupported LES_MODEL_TYPE chosen'
      end if

! Now, add the subgrid scale forcing to CFi
! (This includes the subgrid scale stress as an explicit R-K term

       

      IF (ADD_LES) THEN
  
        IF (LES_IMPLICIT) THEN
         DO J=JSTART,JEND
          DO K=ZSTART,ZEND
            DO I=0,NX2P
             CF1X(I,K,J)=CF1X(I,K,J)                                              &
                     -INT_JACOB(K,J)*CIKXP(I)*CSij(I,K,J,1)                       &
                      -( 0.5*CJOB_12(K+1,J,2)*(CSij(I,K,J,4)  + CSij(I,K+1,J,4))  &
                        - 0.5*CJOB_12(K,J,2)*(CSij(I,K,J,4)   + CSij(I,K-1,J,4))  &
                        + 0.5*CJOB_22(K,J+1,1)*(CSij(I,K,J,4) + CSij(I,K,J+1,4))  &
                        - 0.5*CJOB_22(K,J,1)*(CSij(I,K,J,4)   + CSij(I,K,J-1,4)) )& 

                      -(  0.5*CJOB_11(K+1,J,2)*(CSij(I,K,J,5) + CSij(I,K+1,J,5))  &
                        - 0.5*CJOB_11(K,J,2)*(CSij(I,K,J,5)   + CSij(I,K-1,J,5))  &
                        + 0.5*CJOB_21(K,J+1,1)*(CSij(I,K,J,5) + CSij(I,K,J+1,5))  &
                        - 0.5*CJOB_21(K,J,1)*(CSij(I,K,J,5)   + CSij(I,K,J-1,5)) )             
   

             CF2X(I,K,J)=CF2X(I,K,J)                                              &
                      -INT_JACOB(K,J)*CIKXP(I)*Ctemp_1(i,k,j)                     &
                      -( 0.5*CJOB_12(K+1,J,2)*(CSij(I,K,J,2)  + CSij(I,K+1,J,2))  &
                        - 0.5*CJOB_12(K,J,2)*(CSij(I,K,J,2)   + CSij(I,K-1,J,2))  &
                        + 0.5*CJOB_22(K,J+1,1)*(CSij(I,K,J,2) + CSij(I,K,J+1,2))  &
                        - 0.5*CJOB_22(K,J,1)*(CSij(I,K,J,2)   + CSij(I,K,J-1,2)) )&

                       -(  0.5*CJOB_11(K+1,J,2)*(CSij(I,K,J,6) + CSij(I,K+1,J,6)) &
                        - 0.5*CJOB_11(K,J,2)*(CSij(I,K,J,6)   + CSij(I,K-1,J,6))  &
                        + 0.5*CJOB_21(K,J+1,1)*(CSij(I,K,J,6) + CSij(I,K,J+1,6))  &
                        - 0.5*CJOB_21(K,J,1)*(CSij(I,K,J,6)   + CSij(I,K,J-1,6)) ) 

           
           CF3X(I,K,J)=CF3X(I,K,J)                                                        &
                     -INT_JACOB(K,J)*CIKXP(I)*Ctemp_2(i,k,j)                              &

                      -( 0.5*CJOB_12(K+1,J,2)*(Ctemp_3(I,K,J)  + Ctemp_3(I,K+1,J))        &
                        - 0.5*CJOB_12(K,J,2)*(Ctemp_3(I,K,J)   + Ctemp_3(I,K-1,J))        &
                        + 0.5*CJOB_22(K,J+1,1)*(Ctemp_3(I,K,J) + Ctemp_3(I,K,J+1))        &
                        - 0.5*CJOB_22(K,J,1)*(Ctemp_3(I,K,J)   + Ctemp_3(I,K,J-1)) )      &

                      -(  0.5*CJOB_11(K+1,J,2)*(CSij(I,K,J,3) + CSij(I,K+1,J,3))          &
                        - 0.5*CJOB_11(K,J,2)*(CSij(I,K,J,3)   + CSij(I,K-1,J,3))          &
                        + 0.5*CJOB_21(K,J+1,1)*(CSij(I,K,J,3) + CSij(I,K,J+1,3))          &
                        - 0.5*CJOB_21(K,J,1)*(CSij(I,K,J,3)   + CSij(I,K,J-1,3)) ) 
                           
            END DO
          END DO
        END DO

       ELSE
        DO J=JSTART,JEND
          DO K=ZSTART,ZEND
            DO I=0,NX2P
             CF1X(I,K,J)=CF1X(I,K,J)                                             &
                     -INT_JACOB(K,J)*CIKXP(I)*CSij(I,K,J,1)                       &
                      -( 0.5*CJOB_12(K+1,J,2)*(CSij(I,K,J,4)  + CSij(I,K+1,J,4)) &
                        - 0.5*CJOB_12(K,J,2)*(CSij(I,K,J,4)   + CSij(I,K-1,J,4)) &
                        + 0.5*CJOB_22(K,J+1,1)*(CSij(I,K,J,4) + CSij(I,K,J+1,4)) &
                        - 0.5*CJOB_22(K,J,1)*(CSij(I,K,J,4)   + CSij(I,K,J-1,4)) ) &

                      -(  0.5*CJOB_11(K+1,J,2)*(CSij(I,K,J,5) + CSij(I,K+1,J,5)) &
                        - 0.5*CJOB_11(K,J,2)*(CSij(I,K,J,5)   + CSij(I,K-1,J,5)) &
                        + 0.5*CJOB_21(K,J+1,1)*(CSij(I,K,J,5) + CSij(I,K,J+1,5)) &
                        - 0.5*CJOB_21(K,J,1)*(CSij(I,K,J,5)   + CSij(I,K,J-1,5)) )


             CF2X(I,K,J)=CF2X(I,K,J)                                                &
                      -INT_JACOB(K,J)*CIKXP(I)*CSij(I,K,J,4)                      &
                      -( 0.5*CJOB_12(K+1,J,2)*(CSij(I,K,J,2)  + CSij(I,K+1,J,2)) &
                        - 0.5*CJOB_12(K,J,2)*(CSij(I,K,J,2)   + CSij(I,K-1,J,2)) &
                        + 0.5*CJOB_22(K,J+1,1)*(CSij(I,K,J,2) + CSij(I,K,J+1,2)) &
                        - 0.5*CJOB_22(K,J,1)*(CSij(I,K,J,2)   + CSij(I,K,J-1,2)) ) &

                       -(  0.5*CJOB_11(K+1,J,2)*(CSij(I,K,J,6) + CSij(I,K+1,J,6))&
                        - 0.5*CJOB_11(K,J,2)*(CSij(I,K,J,6)   + CSij(I,K-1,J,6)) &
                        + 0.5*CJOB_21(K,J+1,1)*(CSij(I,K,J,6) + CSij(I,K,J+1,6)) &
                        - 0.5*CJOB_21(K,J,1)*(CSij(I,K,J,6)   + CSij(I,K,J-1,6)) )
      
             CF3X(I,K,J)=CF3X(I,K,J)                                             &
                     -INT_JACOB(K,J)*CIKXP(I)*CSij(I,K,J,5)                       &
                     -( 0.5*CJOB_12(K+1,J,2)*(CSij(I,K,J,6)  + CSij(I,K+1,J,6))  &
                        - 0.5*CJOB_12(K,J,2)*(CSij(I,K,J,6)   + CSij(I,K-1,J,6)) &
                        + 0.5*CJOB_22(K,J+1,1)*(CSij(I,K,J,6) + CSij(I,K,J+1,6)) &
                        - 0.5*CJOB_22(K,J,1)*(CSij(I,K,J,6)   + CSij(I,K,J-1,6)) ) &

                      -(  0.5*CJOB_11(K+1,J,2)*(CSij(I,K,J,3) + CSij(I,K+1,J,3)) &
                        - 0.5*CJOB_11(K,J,2)*(CSij(I,K,J,3)   + CSij(I,K-1,J,3)) &
                        + 0.5*CJOB_21(K,J+1,1)*(CSij(I,K,J,3) + CSij(I,K,J+1,3)) &
                        - 0.5*CJOB_21(K,J,1)*(CSij(I,K,J,3)   + CSij(I,K,J-1,3)) )

            END DO
          END DO
        END DO



       ENDIF
      ENDIF




! Periodically, output mean quantities
      IF ((MOD(TIME_STEP,SAVE_STATS_INT).EQ.0).AND.(RK_STEP.EQ.1)) THEN
! Get plane averages
!       allocate (tke_sgs_t(0:NZ+1,0:NY+1))
!       allocate (tke_sgs_p(0:NZ+1,0:NY+1))
!       allocate (tke_sgs_diss(0:NZ+1,0:NY+1))
!       allocate (Sij_mean(0:NZ+1,0:NY+1,1:6))
!       allocate (TAU_mean(0:NZ+1,0:NY+1,1:6))
!       allocate (NU_T_mean(0:NZ+1,1:NY+1))         



       do j=0,NY+1
        do k=0,NZ+1
         NU_T_MEAN(k,j)=SUM(NU_T(0:min(NXP,NXP_L),k,j))/dble(NX)
        end do
       end do

       CALL MPI_COMBINE_STATS(NU_T_MEAN,NZ+2,NY+2)

!     Output SGS contribution to the TKE equation
!     LES term in the tke equation: -<u_i'dtau_ij'/dx_j>

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!     sgs contribution to the transport: -dtau_ij'u_i'/dx_j 
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

     
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!     sgs contribution to the production: -<tau_ij><S_ij>
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC      

!     To get S_ij for strain rate tensor we need to store tau_ij 
  
      DO ij = 1,6
      
        do j=1,NY
         do k=1,NZ
          do i=0,NX2P
            ctemp(i,k,j) = CSij(i,k,j,ij)   
          enddo
         enddo
        enddo

      if ( (ij.eq. 1).or.(ij.eq.3).or.(ij.eq.5)) then

      CS1X=Ctemp
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
      temp=S1X

       
!        call fft_x_to_physical(ctemp,temp,0,NY+1,0,NZ+1)
        do j=1,NY
         do k=1,NZ
          do i=0,NXP
            Sij(i,k,j,ij) = temp(i,k,j)
          enddo
         enddo
        enddo
      else if (ij.eq.2) then
! we need to take into account Smagorinsky part du2/dy
      CS1X=Ctemp
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
      temp=S1X

!        call fft_x_to_physical(ctemp,temp,0,NY+1,0,NZ+1)
        do j=1,NY
         do k=1,NZ
          do i=0,NXP
            Sij(i,k,j,ij) = temp(i,k,j)-2.d0*C_DYN(k,j)*DELTA_YF(j)**2.d0 &
                           *S2X(i,k,j)*st_rij(i,k,j,ij)
          enddo
         enddo
        enddo
      else if (ij.eq.4) then

      CS1X=Ctemp_1
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
      temp_1=S1X

!        call fft_x_to_physical(ctemp_1,temp_1,0,NY+1,0,NZ+1)
        do j=1,NY
         do k=1,NZ
          do i=0,NXP
            Sij(i,k,j,ij) = temp_1(i,k,j)
          enddo
         enddo
        enddo
      else if (ij.eq.6) then

      CS1X=Ctemp_2
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
      temp_2=S1X


!        call fft_x_to_physical(ctemp_2,temp_2,0,NY+1,0,NZ+1)
        do j=1,NY
         do k=1,NZ
          do i=0,NXP
            Sij(i,k,j,ij) = temp_2(i,k,j)
          enddo
         enddo
        enddo
      endif
       
      ENDDO
      
! Need to calculate Strain Rate

      do ij=1,6
       do j=0,NY+1
        do k=0,NZ+1 
         TAU_mean(k,j,ij) = SUM(Sij(0:min(NXP,NXP_L),k,j,ij))/dble(NX)
        enddo
       enddo
        CALL MPI_COMBINE_STATS(TAU_mean(0,0,ij),NZ+2,NY+2)
      enddo
         
!c      do j=1,NY
!c          check_tau(j) = 0.d0
!c          check_tau(j)=SUM(temp_1(0:NXM,0:NZM,j))/dble(NX*NZ)
!c      enddo 



      do j=1,NY
       do k=1,NZ 
         tke_sgs_p(k,j) = 0.d0
         do ij=1,3 
       tke_sgs_p(k,j)=tke_sgs_p(k,j)-Sij_mean(k,j,ij)*TAU_mean(k,j,ij) 
         enddo 
      enddo
      enddo     

      
      do j=1,NY
       do k=1,NZ
       tke_sgs_p(k,j)= tke_sgs_p(k,j)- 2.0*Sij_mean(k,j,4)*TAU_mean(k,j,4)  &
                                - 2.0*Sij_mean(k,j,6)*TAU_mean(k,j,6)       &
                      - 2.0*Sij_mean(k,j,5)*TAU_mean(k,j,5)  
      enddo
      enddo
        
         


!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!     sgs contribution to the dissipation: -<tau_ij*S_ij>
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      do j=1,NY
       do k=1,NZ
         tke_sgs_diss(k,j) = 0.d0
         do ij=1,6
           if (ij.le.3) then
! We are computing a diagonal term
             cross=1.d0
           else
! We are computing a cross term
            cross=2.d0
           end if
           do i=0,NXP
             temp(i,k,j) =   cross*Sij(i,k,j,ij)*st_rij(i,k,j,ij)
           enddo
       tke_sgs_diss(k,j)=tke_sgs_diss(k,j)+SUM(temp(0:min(NXP,NXP_L),k,j)) &
                 /dble(NX) 
          enddo
       enddo
       enddo      
      CALL MPI_COMBINE_STATS(tke_sgs_diss,NZ+2,NY+2)


       if (rank .eq. 0) then
       k = time_step/SAVE_STATS_INT

      file_name = 'plane_les/data_les_'  &
             //CHAR(MOD(k,100000)/10000+48) &
             //CHAR(MOD(k,10000)/1000+48)   &
             //CHAR(MOD(k,1000)/100+48)     &
             //CHAR(MOD(k,100)/10+48)       &
             //CHAR(MOD(k,10)+48) //        &
             '.pln'

       call plot_les_tecplot(file_name)
       endif 
      

!      deallocate (tke_sgs_diss)
!      deallocate (tke_sgs_p)
!      deallocate (tke_sgs_t)
!      deallocate (Sij_mean)
!      deallocate (TAU_mean)
!      deallocate (NU_T_mean)

      write(6,*) 'Deallocating TKE_LES_vars'
 
      END IF


      RETURN
      END



      subroutine compute_strain_chan

!C This subroutine computes S_ij for the filtered velocity field
!C The input velocity field should be in fourier space in the periodic
!C directions.
!C For use in the LES model in channel flow (2 periodic directions)
use ntypes
use Domain
use Grid
use Fft_var
use TIME_STEP_VAR
use run_variable
use les_chan_var

implicit none       

      integer I,J,K,ij

       
      DO J=JSTART,JEND
        DO K=ZSTART,ZEND
          DO I=0,NX2P
            CSij(I,K,J,1)=CIKXP(I)*CU1X(I,K,J)

! du2/dy=> d[J-1*\zeta_y*u2]/d\zeta + d[J-1*\eta_y*u2]/d\eta

            CSij(I,K,J,2)= ( 0.5*CJOB_12(K+1,J,2)*(CU2X(I,K,J)  + CU2X(I,K+1,J))     &
                        - 0.5*CJOB_12(K,J,2)*(CU2X(I,K,J)       + CU2X(I,K-1,J))     &
                        + 0.5*CJOB_22(K,J+1,1)*(CU2X(I,K,J)     + CU2X(I,K,J+1))     &
                        - 0.5*CJOB_22(K,J,1)*(CU2X(I,K,J)       + CU2X(I,K,J-1)) )/  &
                           INT_JACOB(K,J) 

!  du3/dz => d[J-1*\zeta_x*p]/d\zeta + d[J-1*\eta_x*p]/d\eta
            CSij(I,K,J,3) = ( 0.5*CJOB_11(K+1,J,2)*(CU3X(I,K,J) + CU3X(I,K+1,J))     &
                        - 0.5*CJOB_11(K,J,2)*(CU3X(I,K,J)       + CU3X(I,K-1,J))     &
                        + 0.5*CJOB_21(K,J+1,1)*(CU3X(I,K,J)     + CU3X(I,K,J+1))     &
                        - 0.5*CJOB_21(K,J,1)*(CU3X(I,K,J)       + CU3X(I,K,J-1)) ) / &
                            INT_JACOB(K,J)

!  1/2(du1/dy+du2/dx) =>
            CSij(I,K,J,4)=0.5d0*( ( 0.5*CJOB_12(K+1,J,2)*(CU1X(I,K,J) + CU1X(I,K+1,J))  &
                        - 0.5*CJOB_12(K,J,2)*(CU1X(I,K,J)   + CU1X(I,K-1,J))            &
                        + 0.5*CJOB_22(K,J+1,1)*(CU1X(I,K,J) + CU1X(I,K,J+1))            &
                        - 0.5*CJOB_22(K,J,1)*(CU1X(I,K,J)   + CU1X(I,K,J-1)) ) /        &
                             INT_JACOB(K,J)                                             &
                        + CIKXP(I)*CU2X(I,K,J) )

!  1/2(du1/dz+du3/dx) =>
            CSij(I,K,J,5)=0.5d0*( ( 0.5*CJOB_11(K+1,J,2)*(CU1X(I,K,J) + CU1X(I,K+1,J)) &
                        - 0.5*CJOB_11(K,J,2)*(CU1X(I,K,J)   + CU1X(I,K-1,J))           &
                        + 0.5*CJOB_21(K,J+1,1)*(CU1X(I,K,J) + CU1X(I,K,J+1))           &
                        - 0.5*CJOB_21(K,J,1)*(CU1X(I,K,J)   + CU1X(I,K,J-1))  )/       &
                             INT_JACOB(K,J)                                            &
                        +CIKXP(I)*CU3X(I,K,J) )
    
!  1/2(du2/dz+du3/dy) =>
            CSij(I,K,J,6)=0.5d0*( 0.5*CJOB_11(K+1,J,2)*(CU2X(I,K,J) + CU2X(I,K+1,J)) &
                        - 0.5*CJOB_11(K,J,2)*(CU2X(I,K,J)   + CU2X(I,K-1,J))         &
                        + 0.5*CJOB_21(K,J+1,1)*(CU2X(I,K,J) + CU2X(I,K,J+1))         &
                        - 0.5*CJOB_21(K,J,1)*(CU2X(I,K,J)   + CU2X(I,K,J-1))         &
                                + 0.5*CJOB_12(K+1,J,2)*(CU3X(I,K,J) + CU3X(I,K+1,J)) &
                        - 0.5*CJOB_12(K,J,2)*(CU3X(I,K,J)   + CU3X(I,K-1,J))         &
                        + 0.5*CJOB_22(K,J+1,1)*(CU3X(I,K,J) + CU3X(I,K,J+1))         &
                        - 0.5*CJOB_22(K,J,1)*(CU3X(I,K,J)   + CU3X(I,K,J-1)) )/      &
                              INT_JACOB(K,J)

          END DO
        END DO
      END DO


! Convert rate of strain tensor to physical space
      do ij=1,6
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
! We now have S_ij in Physical space
      RETURN
      END

 
      subroutine les_filter_chan(A,zstart,zend,jstart,jend, filter_type)
! This subroutine applies the les filter to the input field
! The indices to the start and end of the array in the y-direction
! The array that is passed should be in physical space
use ntypes
use Domain
use Grid

      implicit none

      integer i,k,j,zstart,zend,jstart,jend,filter_type
 

      real*8 A(0:NXV-1,0:NZP,0:NY+1)
      real*8 B(0:NXV-1,0:NZP,0:NY+1)

      integer im2(0:NX-1),im1(0:NX-1),ip1(0:NX+1),ip2(0:NX+2)
      integer km2(0:NZ-1),km1(0:NZ-1),kp1(0:NZ+1),kp2(0:NZ+2)

! These are the weights for the filtering operation used
      real*8 W0,W1,W2,Wm1,Wm2,Wm1_j,W0_j,W1_j

! filter type = 1: grid filter operation
! filter type = 2: test filter operation
! The following is for the 3-point trapezoidal rule, alpha*beta=sqrt(6)
      
      If ( filter_type .eq. 1 ) then
       Wm2=0.d0
       Wm1=1.d0/8.d0
       W0=6.d0/8.d0
       W1=1.d0/8.d0
       W2=0.d0
      elseif ( filter_type .eq. 2 ) then
!     else 
!      Wm1_j=1.d0/4.d0  
!      W0_j=1.d0/2.d0
!      W1_j=1.d0/4.d0
! The following is for the 5-point trapezoidal rule, alpha*beta=9
!       Wm2=1.d0/8.d0
!       Wm1=1.d0/4.d0
!       W0=1.d0/4.d0
!       W1=1.d0/4.d0
!       W2=1.d0/8.d0

       Wm2=0.d0
       Wm1=1.d0/4.d0
       W0=1.d0/2.d0
       W1=1.d0/4.d0
       W2=0.d0 
      else
       pause 'Error, unsupported LES_FILTER_TYPE chosen'
      endif

!      NXM=NX-1
!      NZM=NZ-1

!      do j=0,NY+1
!        do k=0,NZM
!          do i=0,NXM
!            B(i,k,j)=A(i,k,j)
!          end do
!        end do
!      end do

! Filter in the periodic directions using cshift
! Note, cshift if not used since it appears to be slower
! Apply filter in the x-direction
!      B=Wm2*CSHIFT(B,-2,1)+Wm1*CSHIFT(B,-1,1)+W0*B+W1*CSHIFT(B,1,1)
!     &       +W2*CSHIFT(B,2,1)

! Filter using more efficient F77 syntax:
! Set up array to loop around periodic directions
      do i=2,NXM
        im2(i)=i-2
      end do
      im2(1)=NXM
      im2(0)=NX-2
      do i=1,NXM
        im1(i)=i-1
      end do
      im1(0)=NXM
      do i=0,NX-2
        ip1(i)=i+1
      end do
      ip1(NXM)=0
      do i=0,NX-3
        ip2(i)=i+2    
      end do
      ip2(NX-2)=0
      ip2(NXM)=1

      do j=jstart,jend
        do k=zstart,zend
          do i=0,NXM
            B(i,k,j)=Wm2*A(im2(i),k,j)+Wm1*A(im1(i),k,j)+W0*A(i,k,j) &
              +W1*A(ip1(i),k,j)+W2*A(ip2(i),k,j)
          end do
        end do  
      end do

      do j=jstart,jend
        do k=zstart,zend
          do i=0,NXM
            A(i,k,j) = B(i,k,j)
          end do
        end do
      end do

      return
      end



      subroutine les_filter_chan_fourier(A,jstart,jend)

      return
      end

! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! C *** Subroutine to write a .plt file (uses TecPlot tecio.a library)
! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
subroutine plot_les_tecplot(file_name)

use ntypes
use Domain
!use run_variable, only : CU1X, CU2X, CU3X
use les_chan_var

      implicit none

      integer mm,nk,k

      integer  i,j,imax,jmax,kmax
      integer  debug,ier,itot
      integer  tecini,tecdat,teczne,tecnod,tecfil,tecend
      integer  visdouble,disdouble
      character*1 nulchar
      CHARACTER*28 file_name

      real(r8),allocatable,dimension(:,:)   :: xpoint, ypoint, th_wh

      allocate (xpoint(0:NZ+1,0:NY+1), ypoint(0:NZ+1,0:NY+1))


      nulchar = char(0)
      debug   = 0
      visdouble = 0
      disdouble = 1
      imax = NZ+2
      jmax = NY+2
      kmax = 1

! c      if (mm .eq. 1) then
      open(21,file='GRID_XZ.dat',form='formatted',status='old')
      xpoint(:,:) =0.0d0
      ypoint(:,:) =0.0d0
      DO J=0,NY+1
       DO K=0,NZ+1
         READ(21,*)xpoint(K,J),ypoint(K,J)
       ENDDO
      ENDDO
      close(21)

      write(6,*)mm,file_name

      open(22,file=file_name,status='unknown',form='unformatted')
      write(22)xpoint,ypoint,U3_bar_les(:,:),U2_bar_les(:,:),&
              U1_bar_les(:,:),C_DYN(:,:),NU_T_mean(:,:)

      close(22)  


      write(6,*) 'end writing in plt format'
      deallocate (xpoint, ypoint)

      return
      end




