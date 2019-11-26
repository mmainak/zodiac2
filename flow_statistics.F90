! C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE SAVE_STATS_CURVI(FINAL)
! C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
use ntypes
use Domain
use Grid
use Fft_var, only :  pi
use TIME_STEP_VAR
use run_variable
use variable_stat
use mpi_var

use mg_vari, only : INIT_FLAG
      
implicit none

LOGICAL FINAL
CHARACTER*35 FNAME
integer i,j,k,n
real(r8) uc, ubulk, vbulk, area,RNUM1
CHARACTER*29 file_name
CHARACTER*3 PID
CHARACTER*5 file_num
LOGICAL     TKE_BUDGET, SAVE_3D, MEAN_ENERGY,XY_PLANE, XZ_PLANE 

INTEGER, DIMENSION(:), ALLOCATABLE :: seed

! C Initialize the random number generator
      CALL RANDOM_SEED(SIZE = K)
      Allocate (seed(1:K))
      seed(1:K)=10
      CALL RANDOM_SEED(PUT = seed)


TKE_BUDGET = .FALSE.
SAVE_3D    = .FALSE.
MEAN_ENERGY= .FALSE.
XY_PLANE   = .FALSE.        
XZ_PLANE   = .FALSE.

 IF ( rank .eq. 0) THEN
open(99,file='out_screen.txt',form='formatted', &
      status='unknown', position='append')
WRITE(99,*) 'Saving flow statistics.'


!WRITE(6,*) 'Saving flow statistics.'
write(*,*) 'Allocate all the tmp arrays'
ENDIF
WRITE(6,*) 'Saving flow statistics.', rank

!write(*,*) 'Allocate all the tmp arrays'


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! need to allocate temp variables
      call allocate_temps
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


! Compute and write out the centerline velocity
      if (int(float(NY)/2.) .eq. float(NY)/2.) then
! IF NY is even
        if (int(float(NZ)/2.) .eq. float(NZ)/2.) then
!          uc=dble(CU3X(0,NZ/2,int(float(NY)/2.)))
        else
!          uc=0.5*(dble(CU3X(0,int(float(NZ)/2.)-1,NY/2.)   &
!                +  CU3X(0,int(float(NZ)/2.),NY/2) ))
        endif  
      else
!        uc=0.5*(dble(CU3X(0,NZ/2,int(float(NY)/2.)-1)) &
!              +dble(CU3X(0,NZ/2,int(float(NY)/2.))))        
      end if
     
      IF (RANK .eq. 0) THEN
!      write(*,*) 'Centerline velocity = ', uc 
! Compute and write out bulk velocity
! Integrat the instantaneous mean profile numerically at GY points
      UBULK=0.
      VBULK=0.
      DO J=1,NY
       DO K=1,NZ
        UBULK=UBULK+0.25*(dble(CU3X(0,K,J))+dble(CU3X(0,K-1,J)) + &
                dble(CU3X(0,K-1,J-1)) + dble(CU3X(0,K,J-1))) 
        VBULK=VBULK+0.25*(dble(CU1X(0,K,J))+dble(CU1X(0,K-1,J)) + &
                dble(CU1X(0,K-1,J-1)) + dble(CU1X(0,K,J-1)))
       ENDDO
      ENDDO
      UBULK=UBULK/real(NY*NZ)  
      VBULK=VBULK/real(NY*NZ)    
! Write out UBULK
      write(*,*) 'UBULK: ',UBULK
      ENDIF


! Save CUi
      do k=0,NZ+1
        do i=0,NX2P
          do j=0,NY+1
            CR1X(i,k,j)=CU1X(i,k,j)
            CR2X(i,k,j)=CU2X(i,k,j)
            CR3X(i,k,j)=CU3X(i,k,j)
!   THIS STEPs ARE REQURIED WHEN DERIVATIVES w.r.t X IS REQURIED LATER
            CF1X(i,k,j)=CU1X(i,k,j)
            CF2X(i,k,j)=CU2X(i,k,j)
            CF3X(i,k,j)=CU3X(i,k,j)
          end do
        end do
      end do
!    TRANSFERRING MEAN TO ALL NODES 
     do k=0,NZ+1
         write(660,*) abs(CU2X(0,k,0)),k
     enddo 
      do k=0,NZ+1
          do j=0,NY+1
            p_mean(k,j)=CPX(0,k,j)
            CALL MPI_BCAST_COMPLEX(p_mean(k,j),1,1)       
          enddo
      enddo

      do k=0,NZV-1
          do j=0,NY+1
            CALL MPI_BCAST_COMPLEX(CR1X(0,k,j),1,1)
            CALL MPI_BCAST_COMPLEX(CR2X(0,k,j),1,1)
            CALL MPI_BCAST_COMPLEX(CR3X(0,k,j),1,1)       
          enddo
      enddo

!      CALL MPI_BCAST_COMPLEX(CR1X(0,:,:),NZV,NY+2)
!      CALL MPI_BCAST_COMPLEX(CR2X(0,:,:),NZV,NY+2)
!      CALL MPI_BCAST_COMPLEX(CR3X(0,:,:),NZV,NY+2)
!      CALL MPI_BCAST_COMPLEX(CF1X(0,:,:),NZV,NY+2) 
      
! Convert to physical space
      

      CALL REAL_FOURIER_TRANS_U1 (.false.)
      CALL REAL_FOURIER_TRANS_U2 (.false.)
      CALL REAL_FOURIER_TRANS_U3 (.false.)
      CALL REAL_FOURIER_TRANS_P (.false.)


!C       call fft_x_to_physical(CU1X,U1X,0,NY+1,0,NZ+1)
!C       call fft_x_to_physical(CU2X,U2X,0,NY+1,0,NZ+1)
!C       call fft_x_to_physical(CU3X,U3X,0,NY+1,0,NZ+1)
!C       call fft_x_to_physical(CPX,PX,0,NY+1,0,NZ+1)

! c        
! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! CC
! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! C        call save_tke_budget
! 
! 
! ! Get the turbulent kinetic energy at each level
     
      do k=0,NZ+1
        do j=0,NY+1
         TKE(k,j)=(dble(CR1X(0,k,j))**2.0+dble(CR2X(0,k,j))**2.0 &
                   + dble(CR3X(0,k,j))**2.0)/                  &
                   (AMP_OMEGA0/OMEGA0)**2.0
        end do
      end do

    

      do k=0,NZ+1 
        do j=0,NY+1
          urms(k,j)=0.d0
          vrms(k,j)=0.d0
          wrms(k,j)=0.d0
          tim(k,j)=TIME/(60.0d0*60.0d0)
      do i=0,min(NXP,NXP_L) 
        urms(k,j)=urms(k,j)+(U1X(i,k,j)-dble(CR1X(0,k,j)))**2
        vrms(k,j)=vrms(k,j)+(U2X(i,k,j)-dble(CR2X(0,k,j)))**2
        wrms(k,j)=wrms(k,j)+(U3X(i,k,j)-dble(CR3X(0,k,j)))**2
      end do
!        urms(k,j)=dsqrt(urms(k,j)/(dble(NX)))
!        vrms(k,j)=dsqrt(vrms(k,j)/(dble(NX)))
!        wrms(k,j)=dsqrt(wrms(k,j)/(dble(NX)))
      end do 
      end do
      

! Compute the Reynolds stress and mean velocity gradient
      do k=1,NZ
      do j=1,NY
        uv(k,j)=0. 
        uw(k,j)=0.
        wv(k,j)=0.
        pv(k,j)=0.
      do i=0,min(NXP,NXP_L)
        uv(k,j)=uv(k,j)+(U1X(i,k,j)-dble(CR1X(0,k,j))) &
         *(U2X(i,k,j)-dble(CR2X(0,k,j)))

        wv(k,j)=wv(k,j)+ &
           (U3X(i,k,j)-dble(CR3X(0,k,j)))* &
           (U2X(i,k,j)-dble(CR2X(0,k,j))) 
 
        uw(k,j)=uw(k,j)+(U1X(i,k,j)-dble(CR1X(0,k,j))) &
         *(U3X(i,k,j)-dble(CR3X(0,k,j))) 

        pu(k,j)=pv(k,j)+(U3X(i,k,j)-dble(CR3X(0,k,j))) &
         *(PX(i,k,j)-dble(p_mean(k,j)))

        pv(k,j)=pv(k,j)+(U2X(i,k,j)-dble(CR2X(0,k,j))) &
         *(PX(i,k,j)-dble(p_mean(k,j)))
      end do
        uv(k,j)=uv(k,j)/(float(NX))
        uw(k,j)=uw(k,j)/(float(NX))
        wv(k,j)=wv(k,j)/(float(NX))
        pu(k,j)=pv(k,j)/(float(NX))
        pv(k,j)=pv(k,j)/(float(NX))
      end do
      end do        

      
! Get the y-derivative of the mean velocity at GYF points
!      do j=1,NY
!        dudy(j)=dble(CR1X(0,0,j+1)-CR1X(0,0,j-1))/(2.*DYF(j))
!        dwdy(j)=dble(CR3X(0,0,j+1)-CR3X(0,0,j-1))/(2.*DYF(j))
!      end do
! Get the y-derivative of the mean velocity at GY points
      do k=1,NZ 
      do j=1,NY
        dudy(k,j)=(0.125*(CJOB_12(K,J,1)+CJOB_12(K,J+1,1)   &
                     + CJOB_12(K,J,2)+CJOB_12(K+1,J,2) )*   &
                 dble(CR2X(0,k+1,j)-CR2X(0,k-1,j))            &
               +  0.125*(CJOB_22(K,J,1)+CJOB_22(K,J+1,1)    &
                     + CJOB_22(K,J,2)+CJOB_22(K+1,J,2) )*   &
                  dble(CR2X(0,k,j+1)-CR2X(0,k,j-1)))/         &
                  INT_JACOB(K,J)
        dwdy(k,j)=(0.125*(CJOB_12(K,J,1)+CJOB_12(K,J+1,1)   &
                     + CJOB_12(K,J,2)+CJOB_12(K+1,J,2) )*  &
                 dble(CR3X(0,k+1,j)-CR3X(0,k-1,j))           &
               +  0.125*(CJOB_22(K,J,1)+CJOB_22(K,J+1,1)   &
                     + CJOB_22(K,J,2)+CJOB_22(K+1,J,2) )*  &
                  dble(CR3X(0,k,j+1)-CR3X(0,k,j-1)))/        &
                  INT_JACOB(K,J)
      end do
      end do

      do k=1,NZ
       dudy(k,0)=dudy(k,1)
       dwdy(k,0)=dwdy(k,1) 
       dudy(k,NY+1)=dudy(k,NY)
       dwdy(k,NY+1)=dwdy(k,NY)       
      enddo
      

      do k=1,NZ              
      do j=1,NY
        dudz(k,j)=(0.125*(CJOB_21(K,J,1)+CJOB_21(K,J+1,1)   &
                     + CJOB_21(K,J,2)+CJOB_21(K+1,J,2) )*  &
                 dble(CR2X(0,k,j+1)-CR2X(0,k,j-1))           &
               +  0.125*(CJOB_11(K,J,1)+CJOB_11(K,J+1,1)   &
                     + CJOB_11(K,J,2)+CJOB_11(K+1,J,2) )*  &
                  dble(CR2X(0,k+1,j)-CR2X(0,k-1,j)))/        &
                  INT_JACOB(K,J)
        dwdz(k,j)=(0.125*(CJOB_21(K,J,1)+CJOB_21(K,J+1,1)   &
                     + CJOB_21(K,J,2)+CJOB_21(K+1,J,2) )*  &
                 dble(CR3X(0,k,j+1)-CR3X(0,k,j-1))           &
               +  0.125*(CJOB_11(K,J,1)+CJOB_11(K,J+1,1)   &
                     + CJOB_11(K,J,2)+CJOB_11(K+1,J,2) )*  &
                  dble(CR3X(0,k+1,j)-CR3X(0,k-1,j)))/        &
                  INT_JACOB(K,J)
      end do
      enddo
      
      do k=1,NZ
       dudz(k,0)=dudz(k,1)
       dwdz(k,0)=dwdz(k,1)
       dudz(k,NY+1)=dudz(k,NY)
       dwdz(k,NY+1)=dwdz(k,NY)
      enddo

      do j=1,NY
       dwdz(1,j)=dwdz(2,j)
       dudz(1,j)=dudz(2,j)
      end do
        
 
     

! Calculate the mean square shear
      do k=1,NZ
      do j=1,NY
        shear(k,j)=0.d0
          do i=0,min(NXP,NXP_L)
            shear(k,j)=shear(k,j)         &
                 +((U1X(i,k,j+1)-U1X(i,k,j-1))/(2.d0*DYF(j)))**2.d0 &
                 +((U3X(i,k,j+1)-U3X(i,k,j-1))/(2.d0*DYF(j)))**2.d0 &
                 +((U1X(i,k+1,j)-U1X(i,k-1,j))/(2.d0*DZF(k)))**2.d0 &
                 +((U3X(i,k+1,j)-U3X(i,k-1,j))/(2.d0*DZF(k)))**2.d0
          end do
        end do
        shear(k,j)=shear(k,j)/dble(NX)
      end do


! First, get the x-component in fourier space
      do j=1,NY
      do k=1,NZ
        omega_x(k,j) =(0.125*(CJOB_12(K,J,1)+CJOB_12(K,J+1,1)  &
                     + CJOB_12(K,J,2)+CJOB_12(K+1,J,2) )*  &
                 dble(CR3X(0,k+1,j)-CR3X(0,k-1,j))           &
               +  0.125*(CJOB_22(K,J,1)+CJOB_22(K,J+1,1)   &
                     + CJOB_22(K,J,2)+CJOB_22(K+1,J,2) )*  &
                  dble(CR3X(0,k,j+1)-CR3X(0,k,j-1)))/        &
                  INT_JACOB(K,J)                           &
                 -(0.125*(CJOB_21(K,J,1)+CJOB_21(K,J+1,1)  &
                     + CJOB_21(K,J,2)+CJOB_21(K+1,J,2) )*  &
                 dble(CR2X(0,k,j+1)-CR2X(0,k,j-1))           &
               +  0.125*(CJOB_11(K,J,1)+CJOB_11(K,J+1,1)   &
                     + CJOB_11(K,J,2)+CJOB_11(K+1,J,2) )*  &
                  dble(CR2X(0,k+1,j)-CR2X(0,k-1,j)))/        &
                  INT_JACOB(K,J)
      end do
      end do
    
      do k=1,NZ
       omega_x(k,0)   =omega_x(k,1)
       omega_x(k,NY+1)=omega_x(k,NY)
      enddo
      
! Convert to physical space
!      call fft_xz_to_physical(CS1,S1X,0,NY+1)
! Get the rms value
!      do j=1,NY
!      omega_x(j)=0.d0
!      do k=1,NZM
!      do i=1,NXP
!        omega_x(j)=omega_x(j)+S1X(i,k,j)**2.d0
!      end do
!      end do
!      omega_x(j)=sqrt(omega_x(j)/(dble(NX-1)*dble(NZ-1)))
!      end do

! Now, get the y-component in fourier space
!      do j=1,NY
!      do k=0,TNKZ
!      do i=0,NX2P
!        CS1(i,k,j)=CIKZ(k)*CR1X(i,k,j)-CIKXP(i)*CR3X(i,k,j)
!      end do
!      end do
!      end do
! Convert to physical space
!      call fft_xz_to_physical(CS1,S1X,0,NY+1)
! Get the rms value
!      do j=1,NY
!      omega_y(j)=0.d0
!      do k=0,NZM
!      do i=0,NXP
!        omega_y(j)=omega_y(j)+S1X(i,k,j)**2.d0
!      end do
!      end do
!      omega_y(j)=sqrt(omega_y(j)/(dble(NX)*dble(NZ)))
!      end do

! Now, get the y-component in fourier space
!      do j=1,NY
!      do k=0,TNKZ
!      do i=0,NX2P
!        CS1(i,k,j)=CIKXP(i)*0.5d0*(CR2X(i,k,j+1)+CR2X(i,k,j))
!     &             -(CR1X(i,k,j+1)-CR1X(i,k,j-1))/(2.d0*DYF(j))
!      end do
!      end do
!        CS1(0,0,j)=CS1(0,0,j)-dudy(j)
!      end do
! Convert to physical space
!      call fft_xz_to_physical(CS1,S1X,0,NY+1)
! Get the rms value
!      do j=1,NY
!      omega_z(j)=0.d0
!      do k=0,NZM
!      do i=0,NXP
!        omega_z(j)=omega_z(j)+S1X(i,k,j)**2.d0
!      end do
!      end do
!      omega_z(j)=sqrt(omega_z(j)/(dble(NX)*dble(NZ)))
!      end do



! Write out the mean statistics at each time as a binary file
      IF ( rank .eq. 0) THEN
      open(66,file='plane_data/time_bulk.txt',form='formatted', &
      status='unknown', position='append' )
      write(6,*) TIME_STEP,TIME/60.0,DELTA_T
      write(66,565) TIME/(60.0d0*60.0d0),DELTA_T, dble(CR3X(0,1,int((NY+1)/2))), dble(CTHX(0,1,int((NY+1)/2),1))+theta_0, &
           dble(CTHX(0,1,int((NY+1)/2),2))+Sal_0
      close(66)
565   format(f12.5,5f13.8)
       
      do k=0,NZ
         write(667,565) xpoint(k,int((NY+1)/2)), dble(CR3X(0,k,int((NY+1)/2))),dble(CR2X(0,k,int((NY+1)/2))), dble(CTHX(0,k,int((NY+1)/2),1)), &
                        dble(CTHX(0,k,int((NY+1)/2),2))
       enddo
       close(667)

      write(99,*) TIME_STEP,TIME/60.0d0,DELTA_T
      ENDIF
 
      count_data = count_data + 1 
      open(60,file='count.txt',form='formatted',status='unknown')
      write(60,*) count_data
      close(60)
 


501   format(6(F25.9,' '))

 
! c      call tkebudget_chan_1   



! Do over the number of passive scalars
      do n=1,N_TH

! Save CTHX
      do k=0,NZ+1
        do i=0,NX2P
          do j=0,NY+1
            CRTHX(i,k,j,n)=CTHX(i,k,j,n)
          end do
        end do
      end do

       do k=0,NZV-1
          do j=0,NY+1
           CALL MPI_BCAST_COMPLEX(CRTHX(0,k,j,N),1,1)
          enddo
       enddo

      end do  ! done with passive scalar do loop

! Convert to physical space
 
      if (N_TH .gt. 0) then
       CALL REAL_FOURIER_TRANS_TH (.false.)
      endif

!C      call fft_x_to_physical(CTHX(0,0,0,n),THX(0,0,0,n),0,NY+1,0,NZ+1)

! Do over the number of passive scalars
      do n=1,N_TH
      
      do j=0,NY+1
      do k=0,NZ+1
        thrms(k,j,n)=0.
      do i=0,min(NXP,NXP_L)
        thrms(k,j,n)=thrms(k,j,n) + (abs(THX(i,k,j,n) &
                -dble(CRTHX(0,k,j,n))))**2.
      end do
!        thrms(k,j,n)=sqrt(thrms(k,j,n)/float(NX))
      end do
      end do
! Compute the Reynolds stress and mean velocity gradient
      do j=1,NY
      do k=0,NZ
        thv(k,j,n)=0.
        thw(k,j,n)=0.
      do i=0,min(NXP,NXP_L)
       thv(k,j,n)=thv(k,j,n)+(THX(i,k,j,n)-dble(CRTHX(0,k,j,n))) &
         *( U2X(i,k,j) - dble(CR2X(0,k,j)) )

       thw(k,j,n)=thw(k,j,n)+(THX(i,k,j,n)-dble(CRTHX(0,k,j,n))) &
         *( U3X(i,k,j) - dble(CR3X(0,k,j)) )
      end do
      thv(k,j,n)=thv(k,j,n)/(float(NX))
      thw(k,j,n)=thw(k,j,n)/(float(NX))
      end do
      end do

      
! Get the y-derivative of the mean scalar 
! C     dt/dy=zeta_y*(T_zeta) + eta_y*(T_eta)
      do k=1,NZ
      do j=1,NY
        dthdy(k,j,n)=(0.125*(CJOB_12(K,J,1)+CJOB_12(K,J+1,1)  &
                     + CJOB_12(K,J,2)+CJOB_12(K+1,J,2) )*     &
                 dble(CRTHX(0,k+1,j,n)-CRTHX(0,k-1,j,n))        &
               +  0.125*(CJOB_22(K,J,1)+CJOB_22(K,J+1,1)      &
                     + CJOB_22(K,J,2)+CJOB_22(K+1,J,2) )*     &
                  dble(CRTHX(0,k,j+1,n)-CRTHX(0,k,j-1,n)))/     &
                  INT_JACOB(K,J)
!        DELTA_N2(k,j,n)=-RI_TAU(N)*dthdy(k,j,n)
       DELTA_N2(k,j,n)=-Ratio_gr*dthdy(k,j,n)
       IF (CONT_STRAT) THEN
       ELSE
        dthdy(k,j,n)= dthdy(k,j,n)+                          &
                   (0.125*(CJOB_12(K,J,1)+CJOB_12(K,J+1,1)   &
                      + CJOB_12(K,J,2)+CJOB_12(K+1,J,2) )*   &
                   (TH_BAR(k+1,j,N)-TH_BAR(k-1,j,N))             &
                +  0.125*(CJOB_22(K,J,1)+CJOB_22(K,J+1,1)    &
                      + CJOB_22(K,J,2)+CJOB_22(K+1,J,2) )*   &
                   (TH_BAR(k,j+1,N)-TH_BAR(k,j-1,N)))/           &
                   INT_JACOB(K,J)
       ENDIF 
      end do
      end do

      do j=1,NY
       dthdy(0,j,n)=dthdy(1,j,n)
       dthdy(NZ+1,j,n)=dthdy(NZ,j,n)

       DELTA_N2(0,j,n)=DELTA_N2(1,j,n)
       DELTA_N2(NZ+1,j,n)=DELTA_N2(NZ,j,n)
      enddo

      do k=0,NZ+1
        dthdy(k,0,n)=2.0*dthdy(k,1,n)-dthdy(k,2,n)
        dthdy(k,NY+1,n)=2.0*dthdy(k,NY,n)-dthdy(k,NY-1,n)

        DELTA_N2(k,0,n)=2.0*DELTA_N2(k,1,n)-DELTA_N2(k,2,n)
        DELTA_N2(k,NY+1,n)=2.0*DELTA_N2(k,NY,n)-DELTA_N2(k,NY-1,n)
      enddo
      
! C     dt/dz=zeta_x*(T_zeta) + eta_x*(T_eta)
      do k=1,NZ
      do j=1,NY
        dthdz(k,j,n)=(0.125*(CJOB_21(K,J,1)+CJOB_21(K,J+1,1) &
                     + CJOB_21(K,J,2)+CJOB_21(K+1,J,2) )*    &
                 dble(CRTHX(0,k,j+1,n)-CRTHX(0,k,j-1,n))       &
               +  0.125*(CJOB_11(K,J,1)+CJOB_11(K,J+1,1)     &
                     + CJOB_11(K,J,2)+CJOB_11(K+1,J,2) )*    &
                  dble(CRTHX(0,k+1,j,n)-CRTHX(0,k-1,j,n)))/    &
                  INT_JACOB(K,J)
      end do
      end do     

      do j=1,NY
       dthdz(0,j,n)=dthdz(1,j,n)
       dthdz(NZ+1,j,n)=dthdz(NZ,j,n)
      enddo

      do k=0,NZ+1
        dthdz(k,0,n)=2.0*dthdz(k,1,n)-dthdz(k,2,n)
        dthdz(k,NY+1,n)=2.0*dthdz(k,NY,n)-dthdz(k,NY-1,n)
      enddo
      
      UBULK=0.
      DO J=1,NY
       DO K=1,NZ
        UBULK=UBULK+0.25*(dble(CRTHX(0,K,J,n))+dble(CRTHX(0,K-1,J,n)) + &
          dble(CRTHX(0,K-1,J-1,n)) + dble(CRTHX(0,K,J-1,n)) )
        IF (CONT_STRAT) THEN 
         IF( dwdy(k,j) .eq. 0.d0 ) THEN
          Rig(k,j) = -Ratio_gr*(dthdy(k,j,n)-1.d0)*(10.0d8)**2.0
         ELSE 
          Rig(k,j) = -Ratio_gr*(dthdy(k,j,n)-1.d0)/  &
                    (dwdy(k,j)**2.0 + dudz(k,j)**2.0)
         ENDIF
        ELSE
         IF( dwdy(k,j) .eq. 0.d0 ) THEN
          Rig(k,j) = -Ratio_gr*(dthdy(k,j,n))*(10.0d8)**2.0
         ELSE
          Rig(k,j) = -Ratio_gr*(dthdy(k,j,n))/         &
                    (dwdy(k,j)**2.0+ dudz(k,j)**2.0)
         ENDIF
        ENDIF

       ENDDO
      ENDDO
      UBULK=UBULK/real(NY*NZ) 

      If (rank .eq. 0) then
! Write out UBULK
      write(*,*) 'THBULK: ',UBULK, 'Ri', RI_TAU(N)
      write(99,*) 'THBULK: ',UBULK, 'Ri', RI_TAU(N)
      endif
! ! Compute the potential energy dissipation, grad(THX) \cdot grad(THX)
! c      do j=1,NY
! c        pe_diss(j,n)=0.d0
! c        do k=0,NZM
! c          do i=0,NXP
! c            pe_diss(j,n)=pe_diss(j,n)
! c     &          +R1X(i,k,j)**2.d0+R2X(i,k,j)**2.d0+R3X(i,k,j)**2.d0
! c          end do
! c        end do
! c        pe_diss(j,n)=pe_diss(j,n)/dble(NX*NZ)
! c      end do

!!!!!!!!!!!!SAVING FULL 3D DATA!!!!!!!!!!!!!!!!!!!!!!!!!!!

!      IF(N==1)THEN
!      IF (SAVE_3D)  THEN
!       IF (USE_MPI) THEN
!        I=1+RANK
!      ELSE
!        I=1
!      ENDIF
!       PID = CHAR(MOD(I,1000)/100+48)     &
!             //CHAR(MOD(I,100)/10+48)    &
!             //CHAR(MOD(I,10)+48)
!
!       k = time_step/SAVE_STATS_INT
!       file_num = CHAR(MOD(k,10000)/1000+48) &
!             //CHAR(MOD(k,1000)/100+48)     &
!             //CHAR(MOD(k,100)/10+48)       &
!             //CHAR(MOD(k,10)+48)//'_'
!
!       file_name = 'plane_3D/data_'//file_num//PID//'.pln'
!
!       write(6,*) file_name, RANK
!       call plane_3D_binary(file_name)
!
!      ENDIF
!      ENDIF 
!!!!!!!!!!END oF SAVING 3D DATA!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
 !     IF (XY_PLANE)THEN
!!!!!!!!!!!!SAVING SPAN-WISE PLANE DATA!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !     IF (USE_MPI) THEN
 !      	I=1+RANK
 !     ELSE
 !     	I=1
 !     ENDIF
!      PID = CHAR(MOD(I,1000)/100+48)     &
!             //CHAR(MOD(I,100)/10+48)    &
!             //CHAR(MOD(I,10)+48)
       
!      k = time_step/SAVE_STATS_INT
!      file_num = CHAR(MOD(k,10000)/1000+48) &
!             //CHAR(MOD(k,1000)/100+48)     &
!             //CHAR(MOD(k,100)/10+48)       &
!             //CHAR(MOD(k,10)+48)//'_'
       
!          k = 308 !(NZ+1)/2
!      file_name = 'yz_plane/span1_'//file_num//PID//'.pln'
!      call plane_XY_binary(file_name,k)     
!	  k = 347 !290!615
!      file_name = 'yz_plane/span2_'//file_num//PID//'.pln'
!      call plane_XY_binary(file_name,k)     
!	  k = 385 !323!736
!      file_name = 'yz_plane/span3_'//file_num//PID//'.pln'
!      call plane_XY_binary(file_name,k)     
!	  k = 485 !362!795
!      file_name = 'yz_plane/span4_'//file_num//PID//'.pln'
!      call plane_XY_binary(file_name,k)     
!      ENDIF 
!!!!!!!!!!END oF SAVING SPAN-WISE PLANE DATA!!!!!!!!!!!!!!!!!!!!!

!C      call FFT_X_TO_FOURIER(THX(0,0,0,n),CTHX(0,0,0,n),0,NY+1,0,NZ+1)

! End do over number of passive scalars, n
      end do

      IF (SAVE_3D)  THEN
       IF (USE_MPI) THEN
        I=1+RANK
      ELSE
        I=1
      ENDIF
       PID = CHAR(MOD(I,1000)/100+48)     &
             //CHAR(MOD(I,100)/10+48)    &
             //CHAR(MOD(I,10)+48)

       k = time_step/SAVE_STATS_INT
       file_num = CHAR(MOD(k,10000)/1000+48) &
             //CHAR(MOD(k,1000)/100+48)     &
             //CHAR(MOD(k,100)/10+48)       &
             //CHAR(MOD(k,10)+48)//'_'

       file_name = 'plane_3D/data_'//file_num//PID//'.pln'

!       write(6,*) file_name, RANK
       call plane_3D_binary(file_name)

      ENDIF
     


      if (rank .eq. 0 ) then
      if(XZ_PLANE) THEN
       k = time_step/SAVE_STATS_INT
       call plane_parav_vel_xy(k)
      endif
      endif 
      if (XY_PLANE) then
       k = time_step/SAVE_STATS_INT
       call plane_parav_zx(k)
      endif
 ! Convert back to Fourier space
     if (N_TH .gt. 0) then 
       CALL REAL_FOURIER_TRANS_TH (.true.)
     endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   calling energy budget 
     IF ( MEAN_ENERGY) THEN
      call energy_budget_curvi_duct
     ENDIF
!   calling tke  budget
      

      CALL MPI_COMBINE_STATS(urms,NZ+2,NY+2)
      CALL MPI_COMBINE_STATS(vrms,NZ+2,NY+2)
      CALL MPI_COMBINE_STATS(wrms,NZ+2,NY+2)
      CALL MPI_COMBINE_STATS(uv,NZ+2,NY+2)
      CALL MPI_COMBINE_STATS(uw,NZ+2,NY+2)
      CALL MPI_COMBINE_STATS(wv,NZ+2,NY+2)
      CALL MPI_COMBINE_STATS(pu,NZ+2,NY+2)
      CALL MPI_COMBINE_STATS(pv,NZ+2,NY+2)
!      CALL MPI_COMBINE_STATS(dudz,NZ+2,NY+2)
!      CALL MPI_COMBINE_STATS(dudy,NZ+2,NY+2)
!      CALL MPI_COMBINE_STATS(dwdz,NZ+2,NY+2)
!      CALL MPI_COMBINE_STATS(dwdy,NZ+2,NY+2)
!      CALL MPI_COMBINE_STATS(omega_x,NZ+2,NY+2)


! Get the bulk rms value

      if(rank .eq. 0) then

      do k=0,NZ+1
      do j=0,NY+1
        urms(k,j)=dsqrt(urms(k,j)/(dble(NX)))
        vrms(k,j)=dsqrt(vrms(k,j)/(dble(NX)))
        wrms(k,j)=dsqrt(wrms(k,j)/(dble(NX)))
      end do
      end do      
  
      urms_b=0.d0
      vrms_b=0.d0
      wrms_b=0.d0
      area  = 0.d0
      do k=1,NZ
      do j=1,NY
        area  = area + INT_JACOB(K,J)
        urms_b=urms_b+0.25d0*(urms(k,j)+urms(k,j-1) &
            +  urms(k-1,j)+urms(k-1,j-1))*INT_JACOB(K,J)
        vrms_b=vrms_b+0.25d0*(vrms(k,j)+vrms(k,j-1) &
            +  vrms(k-1,j)+vrms(k-1,j-1))*INT_JACOB(K,J)
        wrms_b=wrms_b+0.25d0*(wrms(k,j)+wrms(k,j-1) &
            +  wrms(k-1,j)+wrms(k-1,j-1))*INT_JACOB(K,J)
      end do
      enddo
      urms_b=urms_b/(area)
      vrms_b=vrms_b/(area)
      wrms_b=wrms_b/(area)
 

! Write out the bulk rms velocity
      write(*,*) '<U_rms>: ',urms_b
      write(*,*) '<V_rms>: ',vrms_b
      write(*,*) '<W_rms>: ',wrms_b

      write(99,*) '<U_rms>: ',urms_b
      write(99,*) '<V_rms>: ',vrms_b
      write(99,*) '<W_rms>: ',wrms_b

      endif 

     
      

      Do n =1,N_th
      CALL MPI_COMBINE_STATS(thrms(0,0,n),NZ+2,NY+2)
      CALL MPI_COMBINE_STATS(thv(0,0,n),NZ+2,NY+2)
      CALL MPI_COMBINE_STATS(thw(0,0,n),NZ+2,NY+2)
!      CALL MPI_COMBINE_STATS(dthdz(0,0,n),NZ+2,NY+2)
!      CALL MPI_COMBINE_STATS(dthdy(0,0,n),NZ+2,NY+2)
!      CALL MPI_COMBINE_STATS(Rig,NZ+2,NY+2)

      if (rank .eq. 0) then
      do j=0,NY+1
      do k=0,NZ+1
        thrms(k,j,n)=sqrt(thrms(k,j,n)/float(NX))
      end do
      end do
      endif

      enddo

!!!!!!!!!!!!!TKE BUDGET!!!!!!!!!!!!!!!!!!!!!
     IF (TKE_BUDGET) THEN
      call tkebudget_curvi_duct
     ENDIF
!!!!!!!!!!!!END OF TKE BUDGET!!!!!!!!!!!!!!!

! C Convert velocity back to Fourier space
      CALL REAL_FOURIER_TRANS_U1 (.true.) 
      CALL REAL_FOURIER_TRANS_U2 (.true.)
      CALL REAL_FOURIER_TRANS_U3 (.true.)
      CALL REAL_FOURIER_TRANS_P (.true.)

      CF1X = (0.d0,0.d0)
      CF2X = (0.d0,0.d0)
      CF3X = (0.d0,0.d0)     
      
! 
!C       call fft_x_to_fourier(U1X,CU1X,0,NY+1,0,NZ+1)
!C       call fft_x_to_fourier(U2X,CU2X,0,NY+1,0,NZ+1)
!C       call fft_x_to_fourier(U3X,CU3X,0,NY+1,0,NZ+1)
!C       call fft_x_to_fourier(PX,CPX,0,NY+1,0,NZ+1)

!    combine all statistics and send to root==rank-0
      

      if (rank .eq. 0 ) then
      k = time_step/SAVE_STATS_INT
      file_name = 'plane_data/data_tec_'    &
             //CHAR(MOD(k,100000)/10000+48) &
             //CHAR(MOD(k,10000)/1000+48)   &
             //CHAR(MOD(k,1000)/100+48)     &
             //CHAR(MOD(k,100)/10+48)       &
             //CHAR(MOD(k,10)+48) //        &
             '.pln'

       call plane_parav_vel(k)
!      call  plot_tecplot(file_name)  
       
      write(99,*) 'done save_stats chan'
      write(*,*) 'done save_stats chan' 

      close(99)
      endif

!     need to deallocate all the arrays
      call deallocate_temps
!     CCCCCCCCCCCCCCCCCCCCCCCCCC
      if (rank .eq. 0 ) then
      write(*,*) 'Deallocate all the tmp arrays'
      endif

      RETURN
      END
      
      
! C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE SAVE_STATS_DUCT(FINAL)

      RETURN
      END
  
!-------------------------------------C---
subroutine energy_budget_curvi_duct
!-------------------------------------C--



! NOte, it is important to only run this routine after complete R-K
!  time advancement since F1 is overwritten which is needed between R-K steps
! This subroutine should be called in SAVE_STATS_CHAN after computing
! plane averaged statistics, with the velocity in physical space, and 
! CRi containing the velocity in Fourier space
use ntypes
use Domain
use Grid
use Fft_var, only :  CIKXP
use TIME_STEP_VAR
use run_variable
use mg_vari, only : INIT_FLAG
use variable_stat
use mpi_var
implicit none

      integer i,j,k,n
      CHARACTER*31 file_energy


      do j=0,NY+1
       do k=0,NZ+1
! tke_mean defined at GY points
        energy_mean_old(k,j)=energy_mean(k,j)
        energy_mean(k,j)=0.5d0*( dble(CR1X(0,k,j))**2.d0  &
             + dble(CR2X(0,k,j))**2.d0                 &
             + dble(CR3X(0,k,j))**2.d0 )
        tke_1(k,j)=(energy_mean(k,j)-energy_mean_old(k,j)) &
             /(TIME-TIME_old_eng)
       end do
      end do

      time_old_eng = TIME

!     U2*d<E>/dy

      do j=1,NY
       do k=1,NZ
        tke_2_1(k,j)=-dble(CR2X(0,K,J))*  &
                 (0.125*(CJOB_12(K,J,1)+CJOB_12(K,J+1,1)      &
               + CJOB_12(K,J,2)+CJOB_12(K+1,J,2) )*           &
                 (energy_mean(k+1,j)-energy_mean(k-1,j))      &
               +  0.125*(CJOB_22(K,J,1)+CJOB_22(K,J+1,1)      &
               + CJOB_22(K,J,2)+CJOB_22(K+1,J,2) )*           &
                 (energy_mean(k,j+1)-energy_mean(k,j-1)))/    &
                  INT_JACOB(K,J)
       end do
      end do

!     U3*d<E>/dz

      do j=1,NY
      do k=1,NZ
       tke_2_2(k,j)=-dble(CR3X(0,K,J))*                       &
                  (0.125*(CJOB_21(K,J,1)+CJOB_21(K,J+1,1)    &
               + CJOB_21(K,J,2)+CJOB_21(K+1,J,2) )*          &
                  (energy_mean(k,j+1)-energy_mean(k,j-1))    &
               +  0.125*(CJOB_11(K,J,1)+CJOB_11(K,J+1,1)     &
               + CJOB_11(K,J,2)+CJOB_11(K+1,J,2) )*          &
                  (energy_mean(k+1,j)-energy_mean(k-1,j)))/  &
                  INT_JACOB(K,J)
       tke_2(k,j) = tke_2_1(k,j) + tke_2_2(k,j)
      end do
      end do


    


! Get the production 

      do j=1,NY
      do k=1,NZ
!        tke_3(k,j)=uv(k,j)*dUdy(k,j)+wv(k,j)*dWdy(k,j)+vv(k,j)*dVdy(k,j)  
!                 +vw(k,j)*dVdz(k,j)+uw(k,j)*dUdz(k,j)+ww(k,j)*dWdy(k,j)
       tke_3(k,j)=+uv(k,j)*(0.125*(CJOB_12(K,J,1)+CJOB_12(K,J+1,1)     &
               + CJOB_12(K,J,2)+CJOB_12(K+1,J,2) )*                   &
                  dble(CR1X(0,k+1,j)-CR1X(0,k-1,j))                   &
               +  0.125*(CJOB_22(K,J,1)+CJOB_22(K,J+1,1)            &
               + CJOB_22(K,J,2)+CJOB_22(K+1,J,2) )*                 &
                  dble(CR1X(0,k,j+1)-CR1X(0,k,j-1)))/                 &
                  INT_JACOB(K,J)                                    &
            +uw(k,j)*(0.125*(CJOB_21(K,J,1)+CJOB_21(K,J+1,1)        &
                + CJOB_21(K,J,2)+CJOB_21(K+1,J,2) )*                &
                  dble(CR1X(0,k,j+1)-CR1X(0,k,j-1))                   &
                + 0.125*(CJOB_11(K,J,1)+CJOB_11(K,J+1,1)            &
                + CJOB_11(K,J,2)+CJOB_11(K+1,J,2) )*                &
                  dble(CR1X(0,k+1,j)-CR1X(0,k-1,j)))/                 &
                  INT_JACOB(K,J)                                    &
            +vrms(k,j)**2.0*(0.125*(CJOB_12(K,J,1)+CJOB_12(K,J+1,1) &
                + CJOB_12(K,J,2)+CJOB_12(K+1,J,2) )*              &
                  dble(CR2X(0,k+1,j)-CR2X(0,k-1,j))                 &
                +  0.125*(CJOB_22(K,J,1)+CJOB_22(K,J+1,1)         &
                + CJOB_22(K,J,2)+CJOB_22(K+1,J,2) )*              &
                  dble(CR2X(0,k,j+1)-CR2X(0,k,j-1)))/               &
                  INT_JACOB(K,J)                                  &
            +wv(k,j)*(0.125*(CJOB_21(K,J,1)+CJOB_21(K,J+1,1)      &
                + CJOB_21(K,J,2)+CJOB_21(K+1,J,2) )*              &
                  dble(CR2X(0,k,j+1)-CR2X(0,k,j-1))                 &
                + 0.125*(CJOB_11(K,J,1)+CJOB_11(K,J+1,1)          &
                + CJOB_11(K,J,2)+CJOB_11(K+1,J,2) )*              &
                  dble(CR2X(0,k+1,j)-CR2X(0,k-1,j)))/               &
                  INT_JACOB(K,J)                                   &
            +wv(k,j)*(0.125*(CJOB_12(K,J,1)+CJOB_12(K,J+1,1)       &
                + CJOB_12(K,J,2)+CJOB_12(K+1,J,2) )*               &
                  dble(CR3X(0,k+1,j)-CR3X(0,k-1,j))                  &
                +  0.125*(CJOB_22(K,J,1)+CJOB_22(K,J+1,1)          &
                + CJOB_22(K,J,2)+CJOB_22(K+1,J,2) )*               &
                  dble(CR3X(0,k,j+1)-CR3X(0,k,j-1)))/                &
                  INT_JACOB(K,J)                                    &
            +wrms(k,j)**2.0*(0.125*(CJOB_21(K,J,1)+CJOB_21(K,J+1,1) &
                + CJOB_21(K,J,2)+CJOB_21(K+1,J,2) )*              &
                  dble(CR3X(0,k,j+1)-CR3X(0,k,j-1))                &
                + 0.125*(CJOB_11(K,J,1)+CJOB_11(K,J+1,1)          &
                + CJOB_11(K,J,2)+CJOB_11(K+1,J,2) )*              &
                  dble(CR3X(0,k+1,j)-CR3X(0,k-1,j)))/               &
                  INT_JACOB(K,J)   
      end do
      end do


!    viscous diffusion NU*d2tke/dxidxi

      do j=1,NY
       do k=1,NZ
        tke_4(k,j)= NU*(                                       &
            GMAT_22(K,J+1,1)*(energy_mean(K,J+1)-energy_mean(K,J))     &
            - GMAT_22(K,J,1)*(energy_mean(K,J) -energy_mean(K,J-1))   &
            + GMAT_11(K+1,J,2)*(energy_mean(K+1,J)-energy_mean(K,J)) &
            - GMAT_11(K,J,2)*(energy_mean(K,J) - energy_mean(K-1,J))  &
            +  0.25*(GMAT_12(K+1,J,2)*                 &
               (energy_mean(K+1,J+1) + energy_mean(K,J+1)     &
              - energy_mean(K,J-1) - energy_mean(K+1,J-1))    &
              - GMAT_12(K,J,2)*(energy_mean(K-1,J+1) + energy_mean(K,J+1)  &
              - energy_mean(K,J-1) - energy_mean(K-1,J-1)))  &
            + 0.25*(GMAT_12(K,J+1,1)*                  &
                (energy_mean(K+1,J+1) + energy_mean(K+1,J)   &
               - energy_mean(K-1,J+1) - energy_mean(K-1,J))   &
               - GMAT_12(K,J,1)*(energy_mean(K+1,J) + energy_mean(K+1,J-1)  &
               - energy_mean(K-1,J) - energy_mean(K-1,J-1))) )/   &
              INT_JACOB(K,J)
      end do
      end do

      do n=1,N_TH
        do j=1,NY
         do k=1,NZ
!          tke_5(k,j,n)=-RI_TAU(n)*dble(CTHX(0,k,j,n))*dble(CR2X(0,k,j))
           tke_5(k,j,n)=-Ratio_gr*dble(CTHX(0,k,j,n))*dble(CR2X(0,k,j))
         end do
        end do
      end do



! Construct the resolved turbulent transport terms 
! Convert the pressure to physical space for computation of the
! pressure transport
! Get the mean of the pressure
!      do j=0,NY+1
!       do k=0,NZ+1
!        p_mean(k,j)=dble(CF1X(0,k,j))
!       end do
!      end do

      do j=0,NY+1      
      do k=0,NZ+1
        transport(k,j)=0.d0
! Vertical Pressure Transport term:
        transport(k,j)=dble(CR2X(0,k,j))*p_mean(k,J)        
      end do
      end do

      do j=1,NY
       do k=1,NZ
        tke_6_1_1(k,j)=-(0.125*(CJOB_12(K,J,1)+CJOB_12(K,J+1,1) &
                     + CJOB_12(K,J,2)+CJOB_12(K+1,J,2) )* &
                 (transport(k+1,j)-transport(k-1,j))      &
               +  0.125*(CJOB_22(K,J,1)+CJOB_22(K,J+1,1)  &
                     + CJOB_22(K,J,2)+CJOB_22(K+1,J,2) )* &
                 (transport(k,j+1)-transport(k,j-1)))/    &
                  INT_JACOB(K,J)
       end do 
      end do

      do j=0,NY+1      
      do k=0,NZ+1
        transport(k,j)=0.d0
! Horizontal (streamwise) Pressure Transport term:
          transport(k,j)=dble(CR3X(0,k,j))*p_mean(k,J)
      end do
      end do

      do j=1,NY
       do k=1,NZ
        tke_6_1_2(k,j)=-(0.125*(CJOB_21(K,J,1)+CJOB_21(K,J+1,1) &
                     + CJOB_21(K,J,2)+CJOB_21(K+1,J,2) )*       &
                 (transport(k,j+1)-transport(k,j-1))            &
               +  0.125*(CJOB_11(K,J,1)+CJOB_11(K,J+1,1)        &
                     + CJOB_11(K,J,2)+CJOB_11(K+1,J,2) )*       &
                 (transport(k+1,j)-transport(k-1,j)))/          &
                  INT_JACOB(K,J)
        tke_6_1(k,j)=tke_6_1_1(k,j) + tke_6_1_2(k,j)
       end do 
      end do

! Turbulent transport terms:
!   d(0.5*u_i^2.0*v)dy


      do j=0,NY+1
       do k=0,NZ+1
        transport(k,j)=0.d0
        transport(k,j)=transport(k,j)                  &
! <U1*u2>*U1
          + uv(k,j)*dble(CR1X(0,K,J))                   &
! <u3*u2>U3
          + wv(k,j)*dble(CR3X(0,K,J))                   &
! <u2*u2>*U2
          + vrms(k,j)**2.0*dble(CR2X(0,K,J))
       end do
      end do

! Now, the vertical derivative of the transport term:

      do j=1,NY
      do k=1,NZ
        tke_6_2_1(k,j)=-(0.125*(CJOB_12(K,J,1)+CJOB_12(K,J+1,1) &
               + CJOB_12(K,J,2)+CJOB_12(K+1,J,2) )*            &
                 (transport(k+1,j)-transport(k-1,j))           &
               + 0.125*(CJOB_22(K,J,1)+CJOB_22(K,J+1,1)        &
               + CJOB_22(K,J,2)+CJOB_22(K+1,J,2) )*            &
                 (transport(k,j+1)-transport(k,j-1)))/         &
                 INT_JACOB(K,J) 
      end do
      end do


!   d(0.5*u_i^2.0*w)dz

      do j=0,NY+1
       do k=0,NZ+1
         transport(k,j)=0.d0
         transport(k,j)=transport(k,j)   &
! <u1*u3>*U1
          + uw(k,j)*dble(CR1X(0,K,J))     &
! <u2*u3>*U2
          + wv(k,j)*dble(CR2X(0,K,J))     &
! <u3*u3>*U3
          + wrms(k,j)**2.0*dble(CR3X(0,K,J)) 
        end do
      end do

! Now, the horizontal derivative of the transport term:

      do j=1,NY
       do k=1,NZ
        tke_6_2_2(k,j)=-(0.125*(CJOB_21(K,J,1)+CJOB_21(K,J+1,1) &
               + CJOB_21(K,J,2)+CJOB_21(K+1,J,2) )*            &
                 (transport(k,j+1)-transport(k,j-1))           &
               +  0.125*(CJOB_11(K,J,1)+CJOB_11(K,J+1,1)       &
               + CJOB_11(K,J,2)+CJOB_11(K+1,J,2) )*            &
                 (transport(k+1,j)-transport(k-1,j)))/         &
                  INT_JACOB(K,J)
        tke_6_2(k,j) = tke_6_2_1(k,j) + tke_6_2_2(k,j)
       end do
      end do
       
 





 
! Compute the turbulent dissipation rate, epsilon=nu*<du_i/dx_j du_i/dx_j>
      do j=0,NY+1
       do k=0,NZ+1
        epsilon(k,j)=0.d0
       end do
      end do

      
      
! Compute du/dy note remove mean

      do j=1,NY
      do k=1,NZ
       tke_7(k,j) = (0.125*(CJOB_12(K,J,1)+CJOB_12(K,J+1,1) &
               +  CJOB_12(K,J,2)+CJOB_12(K+1,J,2) )*       & 
                  (dble(CR1X(0,k+1,j))-dble(CR1X(0,k-1,j)))  &
               +  0.125*(CJOB_22(K,J,1)+CJOB_22(K,J+1,1)   &
               +  CJOB_22(K,J,2)+CJOB_22(K+1,J,2) )*       &
                  (dble(CR1X(0,k,j+1))-dble(CR1X(0,k,j-1))))/&
                  INT_JACOB(K,J)

       epsilon(k,j)=epsilon(k,j)+(tke_7(k,j)**2.0)
      end do
      end do

      
      
! Compute du/dz 
      do j=1,NY
      do k=1,NZ
        tke_7(k,j) = (0.125*(CJOB_21(K,J,1)+CJOB_21(K,J+1,1)  &
               + CJOB_21(K,J,2)+CJOB_21(K+1,J,2) )*         &
                 (dble(CR1X(0,k,j+1))-dble(CR1X(0,k,j-1)))    &
               + 0.125*(CJOB_11(K,J,1)+CJOB_11(K,J+1,1)     & 
               + CJOB_11(K,J,2)+CJOB_11(K+1,J,2) )*         &
                 (dble(CR1X(0,k+1,j))-dble(CR1X(0,k-1,j))))/  &
                  INT_JACOB(K,J)

        epsilon(k,j)=epsilon(k,j)+ (tke_7(k,j)**2.0)
      end do
      end do
    

! Compute dv/dy  note remove mean

      do j=1,NY
      do k=1,NZ
       tke_7(k,j) = (0.125*(CJOB_12(K,J,1)+CJOB_12(K,J+1,1) &
               +  CJOB_12(K,J,2)+CJOB_12(K+1,J,2) )*      &
                  (dble(CR2X(0,k+1,j))-dble(CR2X(0,k-1,j))) &
               +  0.125*(CJOB_22(K,J,1)+CJOB_22(K,J+1,1)  &
               +  CJOB_22(K,J,2)+CJOB_22(K+1,J,2) )*      &
                  (dble(CR2X(0,k,j+1))-dble(CR2X(0,k,j-1))))/             &
                  INT_JACOB(K,J)
 
       epsilon(k,j)=epsilon(k,j)+ (tke_7(k,j)**2.0)
      end do
      end do

! Compute dw/dy

      do j=1,NY
      do k=1,NZ
       tke_7(k,j)=(0.125*(CJOB_12(K,J,1)+CJOB_12(K,J+1,1)   &
                     + CJOB_12(K,J,2)+CJOB_12(K+1,J,2) )*  &
                  (dble(CR3X(0,k+1,j))-dble(CR3X(0,k-1,j)))  &  
               +  0.125*(CJOB_22(K,J,1)+CJOB_22(K,J+1,1)   &
                     + CJOB_22(K,J,2)+CJOB_22(K+1,J,2) )*  &
                  (dble(CR3X(0,k,j+1))-dble(CR3X(0,k,j-1))))/ &
                  INT_JACOB(K,J)

       epsilon(k,j)=epsilon(k,j)+ (tke_7(k,j)**2.0)
      end do
      end do
     
! Store dv/dz 
      

      do j=1,NY
      do k=1,NZ
        tke_7(k,j) = (0.125*(CJOB_21(K,J,1)+CJOB_21(K,J+1,1) &
               +  CJOB_21(K,J,2)+CJOB_21(K+1,J,2) )*         &
                  (dble(CR2X(0,k,j+1))-dble(CR2X(0,k,j-1)))    & 
               +  0.125*(CJOB_11(K,J,1)+CJOB_11(K,J+1,1)     &
               +  CJOB_11(K,J,2)+CJOB_11(K+1,J,2) )*         &
                  (dble(CR2X(0,k+1,j))-dble(CR2X(0,k-1,j))))/   &
                  INT_JACOB(K,J)

         epsilon(k,j)=epsilon(k,j)+ (tke_7(k,j)**2.0)
      end do
      end do 


! Store dw/dz 
      

      do j=1,NY
      do k=1,NZ
        tke_7(k,j) = (0.125*(CJOB_21(K,J,1)+CJOB_21(K,J+1,1) &
               + CJOB_21(K,J,2)+CJOB_21(K+1,J,2) )*          &
                  (dble(CR3X(0,k,j+1))-dble(CR3X(0,k,j-1)))    &
               +  0.125*(CJOB_11(K,J,1)+CJOB_11(K,J+1,1)     &
                 + CJOB_11(K,J,2)+CJOB_11(K+1,J,2) )*        &
                  (dble(CR3X(0,k+1,j))-dble(CR3X(0,k-1,j))))/   &
                  INT_JACOB(K,J)
 
        epsilon(k,j)=epsilon(k,j)+ (tke_7(k,j)**2.0)
      end do
      end do

      do j=1,NY
      do k=1,NZ
        tke_7(k,j)=-NU*epsilon(k,j)
      end do
      end do

      
      do k=1,NZ
        tke_7(k,0)=tke_7(k,1)
      end do

!      if (rank .eq. 0 ) then
!       k = time_step/SAVE_STATS_INT
!       file_energy = 'plane_energy/data_eng_' &
!             //CHAR(MOD(k,100000)/10000+48)  &
!             //CHAR(MOD(k,10000)/1000+48)    &
!             //CHAR(MOD(k,1000)/100+48)      &
!             //CHAR(MOD(k,100)/10+48)        &
!             //CHAR(MOD(k,10)+48) //         &
!             '.pln'

!       call plot_tec_eng(file_energy )
!      endif




      return
      end


subroutine tkebudget_curvi_duct
! NOte, it is important to only run this routine after complete R-K
!  time advancement since F1 is overwritten which is needed between R-K steps
! This subroutine should be called in SAVE_STATS_CHAN after computing
! plane averaged statistics, with the velocity in physical space, and 
! CRi containing the velocity in Fourier space

use ntypes
use Domain
use Grid
use Fft_var, only : CIKXP
use TIME_STEP_VAR
use run_variable
use mg_vari, only : INIT_FLAG
use variable_stat
use mpi_var      
implicit none

      integer i,j,k,n
      CHARACTER*28 file_tke



! Define working arrays
!       real*8 tke_1(0:NZ+1,0:NY+1)
!       real*8 tke_2(0:NZ+1,0:NY+1)
!       real*8 tke_2_1(0:NZ+1,0:NY+1)
!       real*8 tke_2_2(0:NZ+1,0:NY+1)
!       real*8 tke_3(0:NZ+1,0:NY+1)
!       real*8 tke_3_1(0:NZ+1,0:NY+1)
!       real*8 tke_3_2(0:NZ+1,0:NY+1)
!       real*8 tke_3_3(0:NZ+1,0:NY+1)
!       real*8 tke_4(0:NZ+1,0:NY+1)
!       real*8 tke_5(0:NZ+1,0:NY+1,1:N_TH)
!       real*8 tke_6_1(0:NZ+1,0:NY+1)
!       real*8 tke_6_1_1(0:NZ+1,0:NY+1)
!       real*8 tke_6_1_2(0:NZ+1,0:NY+1)
!       real*8 tke_6_1_3(0:NZ+1,0:NY+1)
!       real*8 tke_6_2(0:NZ+1,0:NY+1)
!       real*8 tke_6_2_1(0:NZ+1,0:NY+1)
!       real*8 tke_6_2_2(0:NZ+1,0:NY+1)
!       real*8 tke_6_2_3(0:NZ+1,0:NY+1)
!       real*8 tke_7(0:NZ+1,0:NY+1)
!       real*8 S1_mean(0:NZ+1,0:NY+1)
!       real*8 p_mean(0:NZ+1,0:NY+1)
!       real*8 transport(0:NZ+1,0:NY+1)

      


      do j=0,NY
       do k=0,NZ
! tke_mean defined at GY points
        tke_mean_old(k,j)=tke_mean(k,j)
        tke_mean(k,j)=0.5d0*( urms(k,j)**2.d0  &
             +vrms(k,j)**2.d0                 &
             + wrms(k,j)**2.d0 ) 
        tke_1(k,j)=(tke_mean(k,j)-tke_mean_old(k,j))/(TIME-TIME_old)
       end do
      end do

      time_old=TIME

!     U2*dtke/dy

      do j=1,NY
       do k=1,NZ
        tke_2_1(k,j)=-dble(CR2X(0,K,J))*  &
                 (0.125*(CJOB_12(K,J,1)+CJOB_12(K,J+1,1)  &
               + CJOB_12(K,J,2)+CJOB_12(K+1,J,2) )*     &
                 (tke_mean(k+1,j)-tke_mean(k-1,j))      &
               +  0.125*(CJOB_22(K,J,1)+CJOB_22(K,J+1,1) &
               + CJOB_22(K,J,2)+CJOB_22(K+1,J,2) )*   &
                 (tke_mean(k,j+1)-tke_mean(k,j-1)))/  &
                  INT_JACOB(K,J)
       end do
      end do

!     U3*dtke/dz

      do j=1,NY
      do k=1,NZ
       tke_2_2(k,j)=-dble(CR3X(0,K,J))*                     &
                  (0.125*(CJOB_21(K,J,1)+CJOB_21(K,J+1,1)  &
               + CJOB_21(K,J,2)+CJOB_21(K+1,J,2) )*        &
                  (tke_mean(k,j+1)-tke_mean(k,j-1))        &
               +  0.125*(CJOB_11(K,J,1)+CJOB_11(K,J+1,1)   &
               + CJOB_11(K,J,2)+CJOB_11(K+1,J,2) )*        & 
                  (tke_mean(k+1,j)-tke_mean(k-1,j)))/      &
                  INT_JACOB(K,J)
       tke_2(k,j) = tke_2_1(k,j) + tke_2_2(k,j)
      end do
      end do
       


! Get the production at GY points
      do j=1,NY
      do k=1,NZ
!        tke_3(k,j)=-uv(k,j)*dUdy(k,j)-wv(k,j)*dWdy(k,j)-vv(k,j)*dVdy(k,j)  
!                 -vw(k,j)*dVdz(k,j)-uw(k,j)*dUdz(k,j)-ww(k,j)*dWdy(k,j)
       tke_3(k,j)=-uv(k,j)*(0.125*(CJOB_12(K,J,1)+CJOB_12(K,J+1,1)     &
               + CJOB_12(K,J,2)+CJOB_12(K+1,J,2) )*                   &
                  dble(CR1X(0,k+1,j)-CR1X(0,k-1,j))                   &
               +  0.125*(CJOB_22(K,J,1)+CJOB_22(K,J+1,1)            &
               + CJOB_22(K,J,2)+CJOB_22(K+1,J,2) )*                 &
                  dble(CR1X(0,k,j+1)-CR1X(0,k,j-1)))/                 &
                  INT_JACOB(K,J)                                    &
            -uw(k,j)*(0.125*(CJOB_21(K,J,1)+CJOB_21(K,J+1,1)         &
                + CJOB_21(K,J,2)+CJOB_21(K+1,J,2) )*                &
                  dble(CR1X(0,k,j+1)-CR1X(0,k,j-1))                    &
                + 0.125*(CJOB_11(K,J,1)+CJOB_11(K,J+1,1)            &
                + CJOB_11(K,J,2)+CJOB_11(K+1,J,2) )*                &
                  dble(CR1X(0,k+1,j)-CR1X(0,k-1,j)))/                 &
                  INT_JACOB(K,J)                                    &
            -vrms(k,j)**2.0*(0.125*(CJOB_12(K,J,1)+CJOB_12(K,J+1,1) &
                + CJOB_12(K,J,2)+CJOB_12(K+1,J,2) )*              &
                  dble(CR2X(0,k+1,j)-CR2X(0,k-1,j))                 &
                +  0.125*(CJOB_22(K,J,1)+CJOB_22(K,J+1,1)         &
                + CJOB_22(K,J,2)+CJOB_22(K+1,J,2) )*              &
                  dble(CR2X(0,k,j+1)-CR2X(0,k,j-1)))/               &
                  INT_JACOB(K,J)                                  &
            -wv(k,j)*(0.125*(CJOB_21(K,J,1)+CJOB_21(K,J+1,1)      &
                + CJOB_21(K,J,2)+CJOB_21(K+1,J,2) )*              &
                  dble(CR2X(0,k,j+1)-CR2X(0,k,j-1))                 &
                + 0.125*(CJOB_11(K,J,1)+CJOB_11(K,J+1,1)          &
                + CJOB_11(K,J,2)+CJOB_11(K+1,J,2) )*              &
                  dble(CR2X(0,k+1,j)-CR2X(0,k-1,j)))/               &
                  INT_JACOB(K,J)                                   &
            -wv(k,j)*(0.125*(CJOB_12(K,J,1)+CJOB_12(K,J+1,1)       &
                + CJOB_12(K,J,2)+CJOB_12(K+1,J,2) )*               &
                  dble(CR3X(0,k+1,j)-CR3X(0,k-1,j))                  &
                +  0.125*(CJOB_22(K,J,1)+CJOB_22(K,J+1,1)          &
                + CJOB_22(K,J,2)+CJOB_22(K+1,J,2) )*               &
                  dble(CR3X(0,k,j+1)-CR3X(0,k,j-1)))/                &
                  INT_JACOB(K,J)                                    &
            -wrms(k,j)**2.0*(0.125*(CJOB_21(K,J,1)+CJOB_21(K,J+1,1) &
                + CJOB_21(K,J,2)+CJOB_21(K+1,J,2) )*              &
                  dble(CR3X(0,k,j+1)-CR3X(0,k,j-1))                &
                + 0.125*(CJOB_11(K,J,1)+CJOB_11(K,J+1,1)          &
                + CJOB_11(K,J,2)+CJOB_11(K+1,J,2) )*              &
                  dble(CR3X(0,k+1,j)-CR3X(0,k-1,j)))/               &
                  INT_JACOB(K,J)   
      end do
      end do

! Get the components of the production
! Use S1X as a working variable
! <u1*u1*dU1/dz>==0.0

! u1*u2*dU1/dy
! 
! c      do j=1,NY
! c        do k=1,NZ
! c          do i=0,NXP
! c            S1X(i,k,j)=(0.125*(CJOB_12(K,J,1)+CJOB_12(K,J+1,1)
! c     &          +  CJOB_12(K,J,2)+CJOB_12(K+1,J,2) )*
! c     &             dble(CR1X(0,k+1,j)-CR1X(0,k-1,j))
! c     &          +  0.125*(CJOB_22(K,J,1)+CJOB_22(K,J+1,1)
! c     &          +  CJOB_22(K,J,2)+CJOB_22(K+1,J,2) )*
! c     &             dble(CR1X(0,k,j+1)-CR1X(0,k,j-1)))/
! c     &             INT_JACOB(K,J)
! c
! c            S1X(i,k,j)=S1X(i,k,j)*(U1(i,k,j)-dble(CR1X(0,k,j)))
! c     &                *(U2(i,k,j)-dble(CR2X(0,k,j)))
! c          end do
! c        end do
! c      end do
! c      do j=1,NY
! c       do k=1,NZ
! c        tke_3_1(k,j)=-SUM(S1X(0:NXP,k,j))/dble(NX)
! c       end do
! c      end do
! c      do k=0,NZ+1 
! c       tke_3_1(k,0)=0.d0
! c       tke_3_1(k,NY+1)=tke_3_1(k,NY)
! c      end do
! 
! ! u1*u3*dU1/dz
! c      do j=1,NY
! c       do k=1,NZ
! c        do i=0,NXP
! c         S1X(i,k,j)=(0.125*(CJOB_21(K,J,1)+CJOB_21(K,J+1,1)
! c     &          + CJOB_21(K,J,2)+CJOB_21(K+1,J,2) )*
! c     &            dble(CR1X(0,k,j+1)-CR1X(0,k,j-1))
! c     &          +  0.125*(CJOB_11(K,J,1)+CJOB_11(K,J+1,1)
! c     &          + CJOB_11(K,J,2)+CJOB_11(K+1,J,2) )*
! c     &            dble(CR1X(0,k+1,j)-CR1X(0,k-1,j)))/
! c     &             INT_JACOB(K,J)
! c         S1X(i,k,j)=S1X(i,k,j)*(U1(i,k,j)-dble(CR1X(0,k,j)))
! c     &              *(U3(i,k,j)-dble(CR3X(0,k,j)))
! c        end do
! c       end do
! c      end do
! 
! C      do j=1,NY
! C       do k=1,NZ
! C        tke_3_1(k,j)=tke_3_1(k,j)-SUM(S1X(0:NXP,k,j))/dble(NX)
! C       end do
! C      end do
     
! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
       do j=1,NY
       do k=1,NZ
        tke_3_1(k,j)=  &
          -wrms(k,j)**2.0*(0.125*(CJOB_21(K,J,1)+CJOB_21(K,J+1,1) &
                + CJOB_21(K,J,2)+CJOB_21(K+1,J,2) )*  &
                  dble(CR3X(0,k,j+1)-CR3X(0,k,j-1))     &
                + 0.125*(CJOB_11(K,J,1)+CJOB_11(K,J+1,1) &
                + CJOB_11(K,J,2)+CJOB_11(K+1,J,2) )*  &
                  dble(CR3X(0,k+1,j)-CR3X(0,k-1,j)))/  &
                  INT_JACOB(K,J)              &
           -vrms(k,j)**2.0*(0.125*(CJOB_12(K,J,1)+CJOB_12(K,J+1,1) &
                + CJOB_12(K,J,2)+CJOB_12(K+1,J,2) )*      &
                  dble(CR2X(0,k+1,j)-CR2X(0,k-1,j))         &
                +  0.125*(CJOB_22(K,J,1)+CJOB_22(K,J+1,1) &
                + CJOB_22(K,J,2)+CJOB_22(K+1,J,2) )*  &
                  dble(CR2X(0,k,j+1)-CR2X(0,k,j-1)))/   &
                  INT_JACOB(K,J)
         
       end do
       end do
! <u2*u1*dU1dx> == 0 

! <u2*u2*dU2/dy>
! c      do j=1,NY
! c        do k=1,NZ
! c          do i=0,NXP
! c            S1X(i,k,j)=(U2(i,k,j)-dble(CR2X(0,k,j)))**2.0
! c     &              *(0.125*(CJOB_12(K,J,1)+CJOB_12(K,J+1,1)
! c     &          +  CJOB_12(K,J,2)+CJOB_12(K+1,J,2) )*
! c     &             dble(CR2X(0,k+1,j)-CR2X(0,k-1,j))
! c     &          +  0.125*(CJOB_22(K,J,1)+CJOB_22(K,J+1,1)
! c     &          +  CJOB_22(K,J,2)+CJOB_22(K+1,J,2) )*
! c     &             dble(CR2X(0,k,j+1)-CR2X(0,k,j-1)))/
! c     &             INT_JACOB(K,J)
! c          end do
! c        end do
! c      end do
! 
! c      do j=1,NY
! c       do k=1,NZ
! c        tke_3_2(k,j)=-SUM(S1X(0:NXP,k,j))/dble(NX)
! c       end do
! c      end do
! c      do k=1,NZ+1 
! c       tke_3_2(k,0)=0.d0
! c       tke_3_2(k,NY+1)=tke_3_2(k,NY)
! c      end do
! 
! ! <u2*u3*dU2dz
! c      do j=1,NY
! c       do k=1,NZ
! c        do i=0,NXP
! c         S1X(i,k,j)=(0.125*(CJOB_21(K,J,1)+CJOB_21(K,J+1,1)
! c     &                + CJOB_21(K,J,2)+CJOB_21(K+1,J,2) )*
! c     &            dble(CR2X(0,k,j+1)-CR2X(0,k,j-1))
! c     &          +  0.125*(CJOB_11(K,J,1)+CJOB_11(K,J+1,1)
! c     &                + CJOB_11(K,J,2)+CJOB_11(K+1,J,2) )*
! c     &             dble(CR2X(0,k+1,j)-CR2X(0,k-1,j)))/
! c     &             INT_JACOB(K,J)
! c         S1X(i,k,j)=S1X(i,k,j)*(U2(i,k,j)-dble(CR2X(0,k,j)))
! c     &                *(U3(i,k,j)-dble(CR3X(0,k,j)))
! c         end do
! c        end do
! c      end do
! 
! c      do j=1,NY
! c       do k=1,NZ
! c        tke_3_2(k,j)=tke_3_2(k,j)-SUM(S1X(0:NXP,k,j))/dble(NX)
! c       end do
! c      end do

       do j=1,NY
       do k=1,NZ
        tke_3_2(k,j)=tke_3(k,j)-tke_3_1(k,j)
       end do
       end do


! <u3*u1*dU3/dx>==0.0

      
! <u3*u2*dU3/dy>

! c      do j=1,NY
! c        do k=1,NZ
! c          do i=0,NXP
! c            S1X(i,k,j)=(0.125*(CJOB_12(K,J,1)+CJOB_12(K,J+1,1)
! c     &          +  CJOB_12(K,J,2)+CJOB_12(K+1,J,2) )*
! c     &             dble(CR3X(0,k+1,j)-CR3X(0,k-1,j))
! c     &          +  0.125*(CJOB_22(K,J,1)+CJOB_22(K,J+1,1)
! c     &          +  CJOB_22(K,J,2)+CJOB_22(K+1,J,2) )*
! c     &             dble(CR3X(0,k,j+1)-CR3X(0,k,j-1)))/
! c     &             INT_JACOB(K,J)
! c
! c            S1X(i,k,j)=S1X(i,k,j)*(U3(i,k,j)-dble(CR3X(0,k,j)))
! c     &                *(U2(i,k,j)-dble(CR2X(0,k,j)))
! c          end do
! c        end do
! c      end do
! 
! c      do j=1,NY
! c       do k=1,NZ
! c        tke_3_3(k,j)=-SUM(S1X(0:NXP,k,j))/dble(NX)
! c       end do
! c      end do
! 
! c      do k=0,NZ+1 
! c       tke_3_3(k,0)=0.d0
! c       tke_3_3(k,NY+1)=tke_3_3(k,NY)
! c      end do
! 
! ! <u3*u3*dU3dz
! 
! c      do j=1,NY
! c       do k=1,NZ
! c        do i=0,NXP
! c         S1X(i,k,j)=(0.125*(CJOB_21(K,J,1)+CJOB_21(K,J+1,1)
! c     &                + CJOB_21(K,J,2)+CJOB_21(K+1,J,2) )*
! c     &            dble(CR3X(0,k,j+1)-CR3X(0,k,j-1))
! c     &          +  0.125*(CJOB_11(K,J,1)+CJOB_11(K,J+1,1)
! c     &                + CJOB_11(K,J,2)+CJOB_11(K+1,J,2) )*
! c     &             dble(CR3X(0,k+1,j)-CR3X(0,k-1,j)))/
! c     &             INT_JACOB(K,J)
! c         S1X(i,k,j)=S1X(i,k,j)*(U3(i,k,j)-dble(CR3X(0,k,j)))
! c     &                *(U3(i,k,j)-dble(CR3X(0,k,j)))
! c         end do
! c        end do
! c      end do
! 
! c      do j=1,NY
! c       do k=1,NZ
! c        tke_3_3(k,j)=tke_3_3(k,j)-SUM(S1X(0:NXP,k,j))/dble(NX)
! c       end do
! c      end do

      do n=1,N_TH
      do j=0,NY+1
       do k=0,NZ+1
        tke_3_3(k,j)=dble(CRTHX(0,k,j,n))+TH_BAR(K,J,N) 
       end do
      end do     
      end do

!    viscous diffusion NU*d2tke/dxidxi

      do j=1,NY
       do k=1,NZ
        tke_4(k,j)= NU*(                                       &
            GMAT_22(K,J+1,1)*(tke_mean(K,J+1)-tke_mean(K,J))     &
            - GMAT_22(K,J,1)*(tke_mean(K,J) - tke_mean(K,J-1))   &
            + GMAT_11(K+1,J,2)*(tke_mean(K+1,J) - tke_mean(K,J)) &
            - GMAT_11(K,J,2)*(tke_mean(K,J) - tke_mean(K-1,J))  &
            +  0.25*(GMAT_12(K+1,J,2)*                 &
               (tke_mean(K+1,J+1) + tke_mean(K,J+1)     &
              - tke_mean(K,J-1) - tke_mean(K+1,J-1))    &
              - GMAT_12(K,J,2)*(tke_mean(K-1,J+1) + tke_mean(K,J+1)  &
              - tke_mean(K,J-1) - tke_mean(K-1,J-1)))  &
            + 0.25*(GMAT_12(K,J+1,1)*                  &
                (tke_mean(K+1,J+1) + tke_mean(K+1,J)   &
               - tke_mean(K-1,J+1) - tke_mean(K-1,J))   &
               - GMAT_12(K,J,1)*(tke_mean(K+1,J) + tke_mean(K+1,J-1)  &
               - tke_mean(K-1,J) - tke_mean(K-1,J-1))) )/   &
              INT_JACOB(K,J)              
      end do 
      end do

      
      do n=1,N_TH

       if (n==N_TH)then
        do j=1,NY
         do k=1,NZ
           tke_5(k,j,n)=(-Ratio_gr_a*thv(k,j,1) +      &
                Ratio_gr_g*thv(k,j,2))*Sin(ANG_BETA)
         end do
        end do
       else
        do j=1,NY
         do k=1,NZ
           tke_5(k,j,n)=-Ratio_gr*thv(k,j,n)
           tke_5_1(k,j)=(-Ratio_gr_a*thw(k,j,1)          &
                          +Ratio_gr_g*thw(k,j,2))/Ratio_gr
         end do
        end do
       endif
      end do
  
! Construct the resolved turbulent transport terms 
! Convert the pressure to physical space for computation of the
! pressure transport
! Get the mean of the pressure
!      do j=0,NY+1
!       do k=0,NZ+1
!        p_mean(k,j)=dble(CF1X(0,k,j))
!       end do
!      end do

      do j=0,NY+1      
      do k=0,NZ+1
        transport(k,j)=0.d0
        do i=0,min(NXP,NXP_L)
          transport(k,j)=transport(k,j)     &
! Vertical Pressure Transport term:
             +(U2X(I,K,J)-dble(CR2X(0,k,j)))*(PX(I,K,J)-p_mean(k,J)) 
          end do        
        transport(k,j)=transport(k,j)/dble(NX)
      end do
      end do

     CALL MPI_COMBINE_STATS(transport,NZ+2,NY+2)

      do j=1,NY
       do k=1,NZ
        tke_6_1_1(k,j)=-(0.125*(CJOB_12(K,J,1)+CJOB_12(K,J+1,1) &
                     + CJOB_12(K,J,2)+CJOB_12(K+1,J,2) )* &
                 (transport(k+1,j)-transport(k-1,j))      &
               +  0.125*(CJOB_22(K,J,1)+CJOB_22(K,J+1,1)  &
                     + CJOB_22(K,J,2)+CJOB_22(K+1,J,2) )* &
                 (transport(k,j+1)-transport(k,j-1)))/    &
                  INT_JACOB(K,J)
       end do 
      end do

      do j=0,NY+1      
      do k=0,NZ+1
        transport(k,j)=0.d0
        do i=0,min(NXP,NXP_L)
          transport(k,j)=transport(k,j)  &
! Horizontal Pressure Transport term:
             +(U3X(I,K,J)-dble(CR3X(0,k,j)))*(PX(I,K,J)-p_mean(k,J))
          end do        
        transport(k,j)=transport(k,j)/dble(NX)
      end do
      end do

      CALL MPI_COMBINE_STATS(transport,NZ+2,NY+2)

      do j=1,NY
       do k=1,NZ
        tke_6_1_2(k,j)=-(0.125*(CJOB_21(K,J,1)+CJOB_21(K,J+1,1) &
                     + CJOB_21(K,J,2)+CJOB_21(K+1,J,2) )*       &
                 (transport(k,j+1)-transport(k,j-1))            &
               +  0.125*(CJOB_11(K,J,1)+CJOB_11(K,J+1,1)        &
                     + CJOB_11(K,J,2)+CJOB_11(K+1,J,2) )*       &
                 (transport(k+1,j)-transport(k-1,j)))/          &
                  INT_JACOB(K,J)
        tke_6_1(k,j)=tke_6_1_1(k,j) + tke_6_1_2(k,j)
       end do 
      end do

! Turbulent transport terms:
!   d(0.5*u_i^2.0*v)dy


      do j=0,NY+1
      do k=0,NZ+1
      transport(k,j)=0.d0
          do i=0,min(NXP,NXP_L)
            transport(k,j)=transport(k,j)             &
! u1^2*u2
          + 0.5d0*(U1X(I,K,J)-dble(CR1X(0,k,J)))**2.d0   &
                 *(U2X(I,K,J)- dble(CR2X(0,k,J)) )       &
! U2^3
          + 0.5d0*(U2X(I,K,J)- dble(CR2X(0,k,J)))**3.d0  &
! U3^2*U2
          + 0.5d0*(U3X(I,K,J)-dble(CR3X(0,k,J)))**2.d0    &
                 *(U2X(I,K,J)- dble(CR2X(0,k,J)))
         end do
        transport(k,j)=transport(k,j)/dble(NX)
      end do
      end do

      CALL MPI_COMBINE_STATS(transport,NZ+2,NY+2)
! Now, the vertical derivative of the transport term:
      do j=1,NY
      do k=1,NZ
        tke_6_2_1(k,j)=-(0.125*(CJOB_12(K,J,1)+CJOB_12(K,J+1,1) &
               + CJOB_12(K,J,2)+CJOB_12(K+1,J,2) )*            &
                 (transport(k+1,j)-transport(k-1,j))           &
               + 0.125*(CJOB_22(K,J,1)+CJOB_22(K,J+1,1)        &
               + CJOB_22(K,J,2)+CJOB_22(K+1,J,2) )*            &
                 (transport(k,j+1)-transport(k,j-1)))/         &
                 INT_JACOB(K,J) 
      end do
      end do


!   d(0.5*u_i^2.0*w)dz

      do j=0,NY+1
      do k=0,NZ+1
      transport(k,j)=0.d0
          do i=0,min(NXP,NXP_L)
            transport(k,j)=transport(k,j)   &
! U1^2*U3
          + 0.5d0*(U1X(I,K,J)-dble(CR1X(0,K,J)))**2.d0 &
                 *(U3X(I,K,J)- dble(CR3X(0,K,J)) ) &
! U3^3
          + 0.5d0*(U3X(I,K,J)- dble(CR3X(0,K,J)))**3.d0 &
! U2^2.0*U3
          + 0.5d0*(U2X(I,K,J)-dble(CR2X(0,K,J)))**2.d0 &
                 *(U3X(I,K,J)- dble(CR3X(0,K,J)))
         end do
        transport(k,j)=transport(k,j)/dble(NX)
        end do
      end do

      CALL MPI_COMBINE_STATS(transport,NZ+2,NY+2)
! Now, the horizontal derivative of the transport term:

      do j=1,NY
       do k=1,NZ
        tke_6_2_2(k,j)=-(0.125*(CJOB_21(K,J,1)+CJOB_21(K,J+1,1) &
               + CJOB_21(K,J,2)+CJOB_21(K+1,J,2) )*            &
                 (transport(k,j+1)-transport(k,j-1))           &
               +  0.125*(CJOB_11(K,J,1)+CJOB_11(K,J+1,1)       &
               + CJOB_11(K,J,2)+CJOB_11(K+1,J,2) )*            &
                 (transport(k+1,j)-transport(k-1,j)))/         &
                  INT_JACOB(K,J)
        tke_6_2(k,j) = tke_6_2_1(k,j) + tke_6_2_2(k,j)
       end do
      end do
       
 
! Compute the turbulent dissipation rate, epsilon=nu*<du_i/dx_j du_i/dx_j>
      do j=0,NY+1
       do k=0,NZ+1
        epsilon(k,j)=0.d0
       end do
      end do

! Store du/dx in CS1
      do j=1,NY
      do k=1,NZ
      do i=0,NX2P
        CS1X(i,k,j)=CIKXP(i)*CF1X(i,k,j) ! WE HAVE ALSO STORED CU1 IN CF1X
!        CS1(i,k,j)=CIKXP(i)*CR1X(i,k,j) 
      end do
      end do
      end do
! Convert to physical space

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

!C      call fft_x_to_physical(CS1,S1X,0,NY+1,0,NZ+1)


      do j=1,NY
      do k=1,NZ
      do i=0,min(NXP,NXP_L)
!        epsilon(k,j)=epsilon(k,j)+(S1X(i,k,j)**2.0)
         epsilon(k,j)=epsilon(k,j) + S1X(i,k,j)**2.0
      end do
      end do
      end do
! Store dv/dx in CS1
      do j=0,NY
      do k=0,NZ
      do i=0,NX2P
        CS1X(i,k,j)=CIKXP(i)*CF2X(i,k,j)   ! WE HAVE ALSO STORED CU2 IN CF2X
!        CS1(i,k,j)=CIKXP(i)*CR2X(i,k,j)
      end do
      end do
      end do

! Convert to physical space
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

!C       call fft_x_to_physical(CS1,S1X,0,NY+1,0,NZ+1)
      do j=0,NY
      do k=0,NZ
      do i=0,min(NXP,NXP_L)
        epsilon(k,j)=epsilon(k,j)+(S1X(i,k,j)**2.0)
      end do
      end do
      end do
      
      
      
! Compute du/dy note remove mean
      do j=0,NY+1
      do k=0,NZ+1
      do i=0,min(NXP,NXP_L)
         S2X(i,k,j)= U1X(i,k,j)-dble(CR1X(0,k,j))
      end do
      end do
      end do

      do j=1,NY
      do k=1,NZ
      do i=0,min(NXP,NXP_L)
       S1X(i,k,j) = (0.125*(CJOB_12(K,J,1)+CJOB_12(K,J+1,1) &
               +  CJOB_12(K,J,2)+CJOB_12(K+1,J,2) )*      & 
                  (S2X(i,k+1,j)-S2X(i,k-1,j))               &
               +  0.125*(CJOB_22(K,J,1)+CJOB_22(K,J+1,1)  &
               +  CJOB_22(K,J,2)+CJOB_22(K+1,J,2) )*      &
                  (S2X(i,k,j+1)-S2X(i,k,j-1)))/             &
                  INT_JACOB(K,J)
      end do
      end do
      end do


      do j=1,NY
      do k=1,NZ
      do i=0,min(NXP,NXP_L)
        epsilon(k,j)=epsilon(k,j)+(S1X(i,k,j)**2.0)
      end do
      end do
      end do
      
      
      
! Store dw/dx in CS1
      do j=1,NY
      do k=1,NZ
      do i=0,NX2P
        CS1X(i,k,j)=CIKXP(i)*CF3X(i,k,j)  ! WE HAVE ALSO STORED CU3 IN CF3X
      end do
      end do
      end do

! Convert to physical space
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

!C      call fft_x_to_physical(CS1,S1X,0,NY+1,0,NZ+1)
      do j=0,NY+1
      do k=0,NZ+1
      do i=0,min(NXP,NXP_L)
         epsilon(k,j)=epsilon(k,j)+ (S1X(i,k,j)**2.0)
      end do
      end do
      end do

      
! Compute du/dz 
      do j=1,NY
      do k=1,NZ
      do i=0,min(NXP,NXP_L)
         S2X(i,k,j)= U1X(i,k,j)-dble(CR1X(0,k,j))
      end do
      end do
      end do

      do j=1,NY
      do k=1,NZ
      do i=0,min(NXP,NXP_L)
        S1X(i,k,j) = (0.125*(CJOB_21(K,J,1)+CJOB_21(K,J+1,1)  &
               + CJOB_21(K,J,2)+CJOB_21(K+1,J,2) )*         &
                 (S2X(i,k,j+1)-S2X(i,k,j-1))                  &
               + 0.125*(CJOB_11(K,J,1)+CJOB_11(K,J+1,1)     & 
               + CJOB_11(K,J,2)+CJOB_11(K+1,J,2) )*         &
                 (S2X(i,k+1,j)-S2X(i,k-1,j)))/                &
                  INT_JACOB(K,J)

      end do
      end do
      end do    

      do j=1,NY
      do k=1,NZ
      do i=0,min(NXP,NXP_L)
         epsilon(k,j)=epsilon(k,j)+ (S1X(i,k,j)**2.0)
      end do
      end do
      end do

! Compute dv/dy  note remove mean

      do j=0,NY+1
      do k=0,NZ+1
      do i=0,min(NXP,NXP_L)
       S2X(i,k,j)= U2X(i,k,j)-dble(CR2X(0,k,j))
      end do
      end do
      end do

      do j=1,NY
      do k=1,NZ
      do i=0,min(NXP,NXP_L)
       S1X(i,k,j) = (0.125*(CJOB_12(K,J,1)+CJOB_12(K,J+1,1) &
               +  CJOB_12(K,J,2)+CJOB_12(K+1,J,2) )*      &
                  (S2X(i,k+1,j)-S2X(i,k-1,j))               &
               +  0.125*(CJOB_22(K,J,1)+CJOB_22(K,J+1,1)  &
               +  CJOB_22(K,J,2)+CJOB_22(K+1,J,2) )*      &
                  (S2X(i,k,j+1)-S2X(i,k,j-1)))/             &
                  INT_JACOB(K,J)
      end do
      end do
      end do

      do j=1,NY
      do k=1,NZ
      do i=0,min(NXP,NXP_L)
        epsilon(k,j)=epsilon(k,j)+(S1X(i,k,j)**2.0)
      end do
      end do
      end do

! Compute dw/dy, note remove mean
      do j=0,NY+1
      do k=0,NZ+1
      do i=0,min(NXP,NXP_L)
       S2X(i,k,j)= U3X(i,k,j)-dble(CR3X(0,k,j))
      end do
      end do
      end do

      do j=1,NY
      do k=1,NZ
      do i=0,min(NXP,NXP_L)
       S1X(i,k,j)=(0.125*(CJOB_12(K,J,1)+CJOB_12(K,J+1,1)   &
                     + CJOB_12(K,J,2)+CJOB_12(K+1,J,2) )* &
                  (S2X(i,k+1,j)-S2X(i,k-1,j))               &  
               +  0.125*(CJOB_22(K,J,1)+CJOB_22(K,J+1,1)  &
                     + CJOB_22(K,J,2)+CJOB_22(K+1,J,2) )* &
                  (S2X(i,k,j+1)-S2X(i,k,j-1)))/             &
                  INT_JACOB(K,J)
      end do
      end do
      end do
     

      do j=1,NY
      do k=1,NZ
      do i=0,min(NXP,NXP_L)
        epsilon(k,j)=epsilon(k,j)+(S1X(i,k,j)**2.0)
      end do
      end do
      end do

! Store dv/dz in CF1X
      do j=0,NY+1
      do k=0,NZ+1
      do i=0,min(NXP,NXP_L)
         S2X(i,k,j)= U2X(i,k,j)-dble(CR2X(0,k,j))
      end do
      end do
      end do

      do j=1,NY
      do k=1,NZ
      do i=0,min(NXP,NXP_L)
        S1X(i,k,j) = (0.125*(CJOB_21(K,J,1)+CJOB_21(K,J+1,1) &
               +  CJOB_21(K,J,2)+CJOB_21(K+1,J,2) )*       &
                  (S2X(i,k,j+1)-S2X(i,k,j-1))                & 
               +  0.125*(CJOB_11(K,J,1)+CJOB_11(K,J+1,1)   &
               +  CJOB_11(K,J,2)+CJOB_11(K+1,J,2) )*       &
                  (S2X(i,k+1,j)-S2X(i,k-1,j)))/              &
                  INT_JACOB(K,J)

      end do
      end do
      end do 

      do j=1,NY
      do k=1,NZ
      do i=0,min(NXP,NXP_L)
        epsilon(k,j)=epsilon(k,j)+(S1X(i,k,j)**2.0)
      end do
      end do
      end do

! Store dw/dz in CS1
      do j=0,NY+1
      do k=0,NZ+1
      do i=0,min(NXP,NXP_L)
         S2X(i,k,j)= U3X(i,k,j)-dble(CR3X(0,k,j))
      end do
      end do
      end do

      do j=1,NY
      do k=1,NZ
      do i=0,min(NXP,NXP_L)
        S1X(i,k,j) = (0.125*(CJOB_21(K,J,1)+CJOB_21(K,J+1,1) &
               + CJOB_21(K,J,2)+CJOB_21(K+1,J,2) )*        &
                  (S2X(i,k,j+1)-S2X(i,k,j-1))                 &
               +  0.125*(CJOB_11(K,J,1)+CJOB_11(K,J+1,1)    &
                 + CJOB_11(K,J,2)+CJOB_11(K+1,J,2) )*       &
                  (S2X(i,k+1,j)-S2X(i,k-1,j)))/              &
                  INT_JACOB(K,J)

      end do
      end do
      end do

      do j=1,NY
      do k=1,NZ
      do i=0,min(NXP,NXP_L)
         epsilon(k,j)=epsilon(k,j)+(S1X(i,k,j)**2.0)
      end do
      end do
      end do

!!!!!! NEED ADD FOR ALL PROCESSORS!!!!!!!!!!!!
      CALL MPI_COMBINE_STATS(epsilon,NZ+2,NY+2)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      do j=1,NY
      do k=1,NZ
        tke_7(k,j)=-NU*epsilon(k,j)/float(NX)
      end do
      end do

      
      do k=1,NZ
        tke_7(k,0)=tke_7(k,1)
      end do

      

      do j=0,NY+1
      do k=0,NZ+1
      do i=0,NX2P
        CF1X(I,K,J)=0.
      end do
      end do
      end do
      

      IF (RANK .eq. 0 ) THEN
      k = time_step/SAVE_STATS_INT
!      file_tke = 'plane_tke/data_tke_'        &
!             //CHAR(MOD(k,100000)/10000+48)  &
!             //CHAR(MOD(k,10000)/1000+48)    &
!             //CHAR(MOD(k,1000)/100+48)      &
!             //CHAR(MOD(k,100)/10+48)        &
!             //CHAR(MOD(k,10)+48) //         &
!             '.plt'

!      call plot_tec_tke(file_tke )  
       call plot_para_tke(k)      
      ENDIF
      return 
       end
