                                             
!!!!!!!!!!!! Nonlinear equation of state

subroutine density_TC(Back_Density,Stat_CAL)

use ntypes
use Domain
use Grid
use run_variable
use variable_stat, only : th_wh
use mpi_var, only : RANK
implicit none

       integer i,j,k,kk,N
       real(r8) :: v01,v02,v03,v04,v05,v06,v07,v08,v09,v10, &
                   v11,v12,v13,v14,v15,v16,v17,v18,v19,     &
                   v20,v21,v22,v23,v24,v25,v26,v27,v28,     &
                   v29,v30,v31,v32,v33,v34,v35,v36
       real(r8) :: SA,CT, sigma0,sqrtSA,v_hat_denominator,v_hat_numerator
       logical  :: Back_Density,Stat_CAL 

v01 =  9.998420897506056e+2;
v02 =  2.839940833161907;
v03 = -3.147759265588511e-2;
v04 =  1.181805545074306e-3;
v05 = -6.698001071123802;
v06 = -2.986498947203215e-2;
v07 =  2.327859407479162e-4;
v08 = -3.988822378968490e-2;
v09 =  5.095422573880500e-4;
v10 = -1.426984671633621e-5;
v11 =  1.645039373682922e-7;


v21 =  1.0;
v22 =  2.775927747785646e-3;
v23 = -2.349607444135925e-5;
v24 =  1.119513357486743e-6;
v25 =  6.743689325042773e-10;
v26 = -7.521448093615448e-3;
v27 = -2.764306979894411e-5;
v28 =  1.262937315098546e-7;
v29 =  9.527875081696435e-10;
v30 = -1.811147201949891e-11;
v31 = -3.303308871386421e-5;
v32 =  3.801564588876298e-7;
v33 = -7.672876869259043e-9;
v34 = -4.634182341116144e-11;
v35 =  2.681097235569143e-12;
v36 =  5.419326551148740e-6;



!CT = 0 ;
!SA = 35 ;

IF (Back_Density) THEN

IF (Stat_CAL) THEN
DO J=0,NY+1
DO K=0,NZ+1
CT =  th_wh(K,J,1)
SA =  th_wh(K,J,2)

sqrtSA = sqrt(SA);
v_hat_denominator = v01 + CT*(v02 + CT*(v03 + v04*CT))  &
             + SA*(v05 + CT*(v06 + v07*CT) &
         + sqrtSA*(v08 + CT*(v09 + CT*(v10 + v11*CT))));

v_hat_numerator = v21 + CT*(v22 + CT*(v23 + CT*(v24 + v25*CT))) &
           + SA*(v26 + CT*(v27 + CT*(v28 + CT*(v29 + v30*CT))) + v36*SA &
       + sqrtSA*(v31 + CT*(v32 + CT*(v33 + CT*(v34 + v35*CT)))));

sigma0 = v_hat_denominator/(v_hat_numerator*rho_0) ;

th_wh(K,J,3)= sigma0

ENDDO
ENDDO

ELSE
DO J=0,NY+1
DO K=0,NZ+1

CT =  TH_BAR(K,J,1) 
SA =  TH_BAR(K,J,2)

sqrtSA = sqrt(SA);
v_hat_denominator = v01 + CT*(v02 + CT*(v03 + v04*CT))  &
             + SA*(v05 + CT*(v06 + v07*CT) &
         + sqrtSA*(v08 + CT*(v09 + CT*(v10 + v11*CT))));

v_hat_numerator = v21 + CT*(v22 + CT*(v23 + CT*(v24 + v25*CT))) &
           + SA*(v26 + CT*(v27 + CT*(v28 + CT*(v29 + v30*CT))) + v36*SA &
       + sqrtSA*(v31 + CT*(v32 + CT*(v33 + CT*(v34 + v35*CT)))));

sigma0 = v_hat_denominator/(v_hat_numerator*rho_0) ;

TH_BAR(K,J,3)= sigma0 

ENDDO
ENDDO
ENDIF

ELSE
DO J=JSTART,JEND
DO K=ZSTART,ZEND
DO I=0,NXP

CT =  TH_BAR(K,J,1) + THX(I,K,J,1) ;
SA =  TH_BAR(K,J,2) + THX(I,K,J,2) ;

sqrtSA = sqrt(SA);
v_hat_denominator = v01 + CT*(v02 + CT*(v03 + v04*CT))  &
             + SA*(v05 + CT*(v06 + v07*CT) &
         + sqrtSA*(v08 + CT*(v09 + CT*(v10 + v11*CT))));

v_hat_numerator = v21 + CT*(v22 + CT*(v23 + CT*(v24 + v25*CT))) &
           + SA*(v26 + CT*(v27 + CT*(v28 + CT*(v29 + v30*CT))) + v36*SA &
       + sqrtSA*(v31 + CT*(v32 + CT*(v33 + CT*(v34 + v35*CT)))));

sigma0 = v_hat_denominator/(v_hat_numerator*rho_0) ;


S1X(I,K,J)= INT_JACOB(K,J)*(sigma0-TH_BAR(K,J,3))*Gravity_g 
     
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

DO J=JSTART,JEND
DO K=ZSTART,ZEND
DO I=0,NX2P
CF2X(I,K,J)=CF2X(I,K,J) - CS1X(I,K,J)
END DO
END DO
END DO

ENDIF

!if(rank .eq. 0) then
!      call plane_parav_sponge
!endif

!write(6,*)'Density is ', sigma0
!stop
end subroutine density_TC
  
! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! C *** Subroutine to write a .plt file (uses TecPlot tecio.a library)
! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!      subroutine plot_tec_tke(file_tke)
!	use ntypes
!	use Domain
!        use Grid
!	use run_variable, only : CR2X, CR3X, TH_BAR
!	use variable_stat

!      implicit none

!      integer  mm,nk,k
!      integer  i,j,imax,jmax,kmax
!      integer  debug,ier,itot
!      integer  tecini,tecdat,teczne,tecnod,tecfil,tecend
!      integer  visdouble,disdouble
!      character*1 nulchar      
!      CHARACTER*28 file_tke

!      real(r8),allocatable,dimension(:,:)   :: xpoint, ypoint
!      allocate (xpoint(0:NZ+1,0:NY+1), ypoint(0:NZ+1,0:NY+1))

!      nulchar = char(0)
!      debug   = 0
!      visdouble = 0
!      disdouble = 1
!      imax = NZ+2
!      jmax = NY+2
!      kmax = 1

! c      if (mm .eq. 1) then
!      open(21,file='GRID_XZ.dat',form='formatted',status='old') 
!      xpoint(:,:) =0.0d0 
!      ypoint(:,:) =0.0d0 
!      DO J=0,NY+1
!       DO K=0,NZ+1
!         READ(21,*)xpoint(K,J),ypoint(K,J)
!       ENDDO
!      ENDDO
!      close(21)

      !
      ! Open the file and write the tecplot datafile header information.
      !
!      write(6,*)mm,file_tke
!      ier = tecini('flow_data'//nulchar,'x,y,ume,wme,tke_m,dkdt,advec,& 
!              &advec_y,advec_z,Prod,Prod_u,&
!              &Prod_v, Prod_w,vis_diff,&
!              &buoy_f,Pr_trns,Pr_trns_y,&
!              &Pr_trns_z,tur_trns,&
!              &tur_trns_y,tur_trns_z,dissip'&
!                   &//nulchar,&
!                  &file_tke//nulchar,&
!                  &'.'//nulchar,&
!                  &debug,visdouble)
! c      endif
 
!     write(6,*) 'strat writing in plt format'
      !
      ! Write the zone header information.
      !
!      ier = teczne('flow_data'//nulchar,  &
!                  imax,jmax,kmax,         &
!                  'BLOCK'//nulchar,nulchar)

      !
      ! Write out the field data.

       
!       open(23,file=file_tke,status='unknown',form='unformatted')
!       write(23)xpoint,ypoint,dble(CR3X(0,0:NZ+1,0:NY+1)),dble(CR2X(0,0:NZ+1,0:NY+1)), &
!                tke_mean,tke_1,tke_2,tke_2_1,tke_2_2,tke_3,tke_3_1, &
!                tke_3_2,tke_3_3,tke_4,tke_5(:,:,1),tke_6_1,tke_6_1_1, &
!                tke_6_1_2,tke_6_2,tke_6_2_1,tke_6_2_2,tke_7


!       close(23)


!      itot = imax*jmax*kmax
!      ier = tecdat(itot,xpoint,disdouble)
!      ier = tecdat(itot,ypoint,disdouble)
!      ier = tecdat(itot,dble(CR3X(0,:,:)),disdouble)
!      ier = tecdat(itot,dble(CR2X(0,:,:)),disdouble)
!      ier = tecdat(itot,tke_mean,disdouble)
!      ier = tecdat(itot,tke_1,disdouble)
!      ier = tecdat(itot,tke_2,disdouble)
!      ier = tecdat(itot,tke_2_1,disdouble)
!      ier = tecdat(itot,tke_2_2,disdouble)
!      ier = tecdat(itot,tke_3,disdouble)
!      ier = tecdat(itot,tke_3_1,disdouble)
!      ier = tecdat(itot,tke_3_2,disdouble)
!      ier = tecdat(itot,tke_3_3,disdouble)
!      ier = tecdat(itot,tke_4,disdouble)
!      ier = tecdat(itot,tke_5(:,:,1),disdouble)
!      ier = tecdat(itot,tke_6_1,disdouble)
!      ier = tecdat(itot,tke_6_1_1,disdouble)
!      ier = tecdat(itot,tke_6_1_2,disdouble)
!      ier = tecdat(itot,tke_6_2,disdouble)
!      ier = tecdat(itot,tke_6_2_1,disdouble)
!      ier = tecdat(itot,tke_6_2_2,disdouble)
!      ier = tecdat(itot,tke_7,disdouble)

! c      if (mm .eq.nk) then 
!       !
!       ! close file
      !
!      ier = tecend()
! c      endif 
!      deallocate (xpoint, ypoint!)
!      return
!      end

! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! C *** Subroutine to write a 3D.PLN file of the instantaneous data
! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

!subroutine plane_3D_binary(file_name)

!use ntypes
!use Domain
!use Grid, only : dx
!use run_variable, only : U1X,U2X,U3X,THX,TIME,TIME_STEP
!use mg_vari, only : INIT_FLAG
!use mpi_var, only : RANK
!use TIME_STEP_VAR, only: DELTA_T      
!use mpi_stitch
!implicit none

!      integer  ind
!      integer  i,j,k, NI,kmax
!      CHARACTER*29 file_name
!      real(r8),allocatable,dimension(:,:)   :: xpoint, ypoint
!      allocate (xpoint(0:NZ+1,0:NY+1), ypoint(0:NZ+1,0:NY+1))

!       write(6,*) file_name, RANK

!      open(21,file='GRID_XZ.dat',form='formatted',status='old')
!      xpoint(:,:) =0.0d0 
!      ypoint(:,:) =0.0d0 
!      DO J=0,NY+1
!       DO I=0,NZ+1
!         READ(21,*)xpoint(I,J),ypoint(I,J)
!       ENDDO
!      ENDDO
!      close(21)


!        NI = min(NXP,NXP_L)
 
!       write(6,*) file_name, RANK,(NI+1)*(NZ+2)*(NY+2) 
!        write(6,*) 'Saving span-wise plane information'
!       open(22,file=file_name,form='unformatted',status='unknown')
!       write(22) TIME,DELTA_T,TIME_STEP,NI+1,NZ+2,NY+2
!       write(22) dx(1),U3X(0:NI,0:NZ+1,0:NY+1),               &
!            U2X(0:NI,0:NZ+1,0:NY+1), U1X(0:NI,0:NZ+1,0:NY+1), &
!            THX(0:NI,0:NZ+1,0:NY+1,1)
!       close(22)

!       write(6,*) 'Done 3D data saving in = ', file_name, rank
!      deallocate (xpoint, ypoint)

!      return
!      end subroutine plane_3D_binary

! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! C *** Subroutine to write a .PLN file of the instantaneous data
! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

!subroutine plane_XY_binary(file_name,ind)

!use ntypes
!use Domain
!use Grid, only : dx,xpoint,ypoint
!use run_variable, only : U1X,U2X,U3X,THX,PX,TIME,TIME_STEP
!use mg_vari, only : INIT_FLAG
!use TIME_STEP_VAR, only: DELTA_T      
!implicit none

!      integer  ind
!      integer  i,j,k, NI
!      CHARACTER*29 file_name
!      real(r8),allocatable,dimension(:,:)   :: xpoint, ypoint
!      allocate (xpoint(0:NZ+1,0:NY+1), ypoint(0:NZ+1,0:NY+1))

!       write(6,*) file_name

!      open(21,file='GRID_XZ.dat',form='formatted',status='old')
!      xpoint(:,:) =0.0d0 
!      ypoint(:,:) =0.0d0 
!      DO J=0,NY+1
!       DO I=0,NZ+1
!         READ(21,*)xpoint(I,J),ypoint(I,J)
!       ENDDO
!      ENDDO
!      close(21)

!	   k=ind
!	   NI = min(NXP,NXP_L)
	   
!        write(6,*) 'Saving span-wise plane information'
!       open(22,file=file_name,form='unformatted',status='unknown')
!       write(22) TIME,DELTA_T,TIME_STEP,xpoint(k,0), NI+1,NY+2
!       write(22) ypoint(k,0:NY+1),dx(1),U3X(0:NI,k,0:NY+1), &
!       			U2X(0:NI,k,0:NY+1), U1X(0:NI,k,0:NY+1), &
!       			PX(0:NI,k,0:NY+1),THX(0:NI,k,0:NY+1,1)
!       close(22)
 
!      deallocate (xpoint, ypoint)

!      return
!      end subroutine plane_XY_binary


	subroutine plane_parav_sponge

use ntypes
use Domain
use Grid, only : dx,xpoint,ypoint
use run_variable
implicit none
             
       integer i,j,k,np1,np2
       !FILENAMES
       character(len=100) :: InFileName, OutFileName
       character(len=100) :: dataDIR, outputDIR
       character(len=80) :: basename
       !OTHER STRINGS
       character(len=25) :: ss
       real(r4),allocatable,dimension(:) :: g1vtk,g2vtk
       real(r4),allocatable,dimension(:,:) :: var1
       !STATUS VARIABLES
       integer :: s1

       allocate(g1vtk(0:NZ+1), g2vtk(0:NY+1))
       allocate(var1(0:NZ+1,0:NY+1))

       np1 = NZ+2
       np2 = NY+2
!       allocate(var1(1:3*np1*np2)) 

!       do j=0,NY+1
!         g1vtk(j) = xpoint(k,j)
!       enddo

       do j=0,NY+1
       do k=0,NZ+1
         g1vtk(k) = xpoint(k,j)
         g2vtk(j) = ypoint(k,j)
         var1(k,j)= SPONGE_SIGMA_OUT(k,j)
       enddo
       enddo
       outputDIR= 'plane_data/'
       basename = 'sponde_density'

       write(OutFileName,'(a,a,i5.5,a4)') trim(outputDIR),trim(basename)//".vtk"

        open(unit=13,file=OutFileName,access='stream',form='unformatted',status='unknown',&
        convert='big_endian',iostat=s1)
!HEADER: note termination with char(10)
        
        write(13) "# vtk DataFile Version 3.0"//char(10)
        write(13) trim(BaseName)//char(10)
        write(13) "BINARY"//char(10)
        write(13) "DATASET RECTILINEAR_GRID"//char(10)

        write(ss,fmt='(A10,3I5)') "DIMENSIONS",1,np1,np2
        write(13) ss//char(10)
!Xgrid
        write(ss,fmt='(A13,I6,A6)') "X_COORDINATES",1," float"
        write(13) char(10)//ss//char(10)
        write(13) 1.0
!Ygrid
        write(ss,fmt='(A13,I6,A6)') "Y_COORDINATES",np1," float"
        write(13) char(10)//ss//char(10)
        do k = 0, np1-1
        write(13) g1vtk(k)
        enddo
!Zgrid
        write(ss,fmt='(A13,I6,A6)') "Z_COORDINATES",np2," float"
        write(13) char(10)//ss//char(10)
        do j = 0, np2-1
          write(13) g2vtk(j)
        enddo

!Field
        write(ss,fmt='(A10,I15)') "POINT_DATA",np1*np2
        write(13) char(10)//ss//char(10)
        write(13) "SCALARS Sponge float 1"//char(10)
        write(13) "LOOKUP_TABLE default"//char(10)
        write(13) real(SPONGE_SIGMA_OUT)
        write(13) "SCALARS Sponge_extra_1 float 1"//char(10)
        write(13) "LOOKUP_TABLE default"//char(10)
        write(13) real(SPONGE_SIGMA)
        write(13) "SCALARS Sponge_extra_2 float 1"//char(10)
        write(13) "LOOKUP_TABLE default"//char(10)
        write(13) real(SPONGE_TEMP)
        write(13) "SCALARS Temerature float 1"//char(10)
        write(13) "LOOKUP_TABLE default"//char(10)
        write(13) real(TH_BAR(0:NZ+1,0:NY+1,1))
        write(13) "SCALARS Salinity float 1"//char(10)
        write(13) "LOOKUP_TABLE default"//char(10)
        write(13) real(TH_BAR(0:NZ+1,0:NY+1,2))
        write(13) "SCALARS Density float 1"//char(10)
        write(13) "LOOKUP_TABLE default"//char(10)
        write(13) real(TH_BAR(0:NZ+1,0:NY+1,3))
!Close VTK File
        close(13)

deallocate(g1vtk,g2vtk,var1 )
!if (allocated(SPplane) ) deallocate(SPplane)
!if (allocated(DPplane) ) deallocate(DPplane)

return
end subroutine plane_parav_sponge 


subroutine plane_parav_vel_xy(kk)

use ntypes
use Domain
use Grid, only : dx,xpoint,ypoint
use run_variable
use variable_stat
implicit none
             
       integer i,j,k,np1,np2,kk,N
       !FILENAMES
       character(len=100) :: InFileName, OutFileName
       character(len=100) :: dataDIR, outputDIR
       character(len=80) :: basename
       !OTHER STRINGS
       character(len=25) :: ss
       real(r4),allocatable,dimension(:) :: g1vtk,g2vtk
       real(r4),allocatable,dimension(:) :: var1 
!       real(r8),allocatable,dimension(:,:,:) :: th_wh
       !STATUS VARIABLES
       integer :: s1

       np1 = NZ+2
       np2 = NY+2
!       allocate (th_wh(0:NZ+1,0:NY+1,1:N_TH+1))
       allocate (var1(1:3*np1*np2))
       allocate(g1vtk(0:NZ+1), g2vtk(0:NY+1))


!       do j=0,NY+1
!         g1vtk(j) = xpoint(k,j)
!       enddo


       do j=0,NY+1
       do k=0,NZ+1
         g1vtk(k) = xpoint(k,j)
         g2vtk(j) = ypoint(k,j)
         DO N=1,N_TH
         th_wh(k,j,N)=(THX(int(NXP/2),k,j,N)) +  TH_BAR(K,j,N) !
!         tim(k,j,N)=TIME
        ENDDO
       enddo
       enddo

       IF (Non_linear_ST) THEN
          call density_TC(.true.,.true.)
       ELSE 
       do j=0,NY+1
       do k=0,NZ+1
         th_wh(k,j,N_TH+1)= -alpha_w*(THX(int(NXP/2),K,J,1)) & 
               + gamma_w*(THX(int(NXP/2),K,J,2)) +  TH_BAR(K,j,N_TH+1) !
       enddo
       enddo
       ENDIF


        k = 1
        do j = 1, np2
        do i = 1, np1
        var1(k) =  real(U1X(int(NXP/2),i-1,j-1))
        k = k+1
        var1(k) =  real(U3X(int(NXP/2),i-1,j-1))
        k = k+1
        var1(k) =  real(U2X(int(NXP/2),i-1,j-1))
        k = k+1
        enddo
        enddo

       outputDIR='plane_data_xy/'
       basename = 'data_parav_xy_'


       write(OutFileName,'(a,a,i5.5,a4)') trim(outputDIR),trim(basename)//"_n",kk,".vtk"

      	open(unit=13,file=OutFileName,access='stream',form='unformatted',status='unknown',&
	convert='big_endian',iostat=s1)


        WRITE(6,*) 'Wrting in ', OutFileName,' for paraview'
!HEADER: note termination with char(10)
	
	write(13) "# vtk DataFile Version 3.0"//char(10)
	write(13) trim(BaseName)//char(10)
	write(13) "BINARY"//char(10)
	write(13) "DATASET RECTILINEAR_GRID"//char(10)

	write(ss,fmt='(A10,3I5)') "DIMENSIONS",1,np1,np2
	write(13) ss//char(10)
!Xgrid
	write(ss,fmt='(A13,I6,A6)') "X_COORDINATES",1," float"
	write(13) char(10)//ss//char(10)
	write(13) 1.0
!Ygrid
	write(ss,fmt='(A13,I6,A6)') "Y_COORDINATES",np1," float"
	write(13) char(10)//ss//char(10)
	do k = 0, np1-1
	write(13) g1vtk(k)
	enddo
!Zgrid
	write(ss,fmt='(A13,I6,A6)') "Z_COORDINATES",np2," float"
	write(13) char(10)//ss//char(10)
	do j = 0, np2-1
	  write(13) g2vtk(j)
	enddo

!Field
        write(ss,fmt='(A10,I15)') "POINT_DATA",np1*np2
        write(13) char(10)//ss//char(10)
        write(13) "VECTORS velocity_vectors float 1"//char(10)
        write(13) var1
!	write(ss,fmt='(A10,I15)') "POINT_DATA",np1*np2
!	write(13) char(10)//ss//char(10)
!	write(13) "SCALARS U1 float 1"//char(10)
!	write(13) "LOOKUP_TABLE default"//char(10)
!	write(13) real(dble(CU3X(0,0:NZ+1,0:NY+1)))
!        write(13) "SCALARS U2 float 1"//char(10)
!	write(13) "LOOKUP_TABLE default"//char(10)
!	write(13) real(dble(CU2X(0,0:NZ+1,0:NY+1)))
!        write(13) "SCALARS U3 float 1"//char(10)
!        write(13) "LOOKUP_TABLE default"//char(10)
!        write(13) real(dble(CU1X(0,0:NZ+1,0:NY+1))) 
        write(13) "SCALARS Temp_d float 1"//char(10)
        write(13) "LOOKUP_TABLE default"//char(10)
        write(13) real(THX(int(NXP/2),0:NZ+1,0:NY+1,1))
        write(13) "SCALARS Salinity_d float 1"//char(10)
        write(13) "LOOKUP_TABLE default"//char(10)
        write(13) real(THX(int(NXP/2),0:NZ+1,0:NY+1,2))
        write(13) "SCALARS Temp_w float 1"//char(10)
        write(13) "LOOKUP_TABLE default"//char(10)
        write(13) real(th_wh(0:NZ+1,0:NY+1,1))
        write(13) "SCALARS Salinity_w float 1"//char(10)
        write(13) "LOOKUP_TABLE default"//char(10)
        write(13) real(th_wh(0:NZ+1,0:NY+1,2))
        write(13) "SCALARS Density_w float 1"//char(10)
        write(13) "LOOKUP_TABLE default"//char(10)
        write(13) real(th_wh(0:NZ+1,0:NY+1,3))

!Close VTK File
	close(13)

deallocate(g1vtk,g2vtk,var1 )
!if (allocated(SPplane) ) deallocate(SPplane)
!if (allocated(DPplane) ) deallocate(DPplane)

return
end subroutine plane_parav_vel_xy 



	subroutine plane_parav_vel(kk)

use ntypes
use Domain
use Grid, only : dx,xpoint,ypoint
use run_variable
use variable_stat
implicit none
             
       integer i,j,k,np1,np2,kk,N
       !FILENAMES
       character(len=100) :: InFileName, OutFileName
       character(len=100) :: dataDIR, outputDIR
       character(len=80) :: basename
       !OTHER STRINGS
       character(len=25) :: ss
       real(r4),allocatable,dimension(:) :: g1vtk,g2vtk
       real(r4),allocatable,dimension(:) :: var1 
!       real(r8),allocatable,dimension(:,:,:) :: tim
       !STATUS VARIABLES
       integer :: s1

       np1 = NZ+2
       np2 = NY+2
!       allocate (tim(0:NZ+1,0:NY+1))
       allocate (var1(1:3*np1*np2))
       allocate(g1vtk(0:NZ+1), g2vtk(0:NY+1))


!       do j=0,NY+1
!         g1vtk(j) = xpoint(k,j!)
!       enddo
!        open(21,file='tim.dat',position='append',form='unformatted',status='unknown')

       do j=0,NY+1
       do k=0,NZ+1
         g1vtk(k) = xpoint(k,j)
         g2vtk(j) = ypoint(k,j)
         DO N=1,N_TH
         th_wh(k,j,N)=dble(CTHX(0,k,j,N)) +  TH_BAR(K,j,N) !
         ENDDO
       enddo
       enddo

       IF (Non_linear_ST) THEN
          call density_TC(.true.,.true.)
       ELSE 
       do j=0,NY+1
       do k=0,NZ+1
         th_wh(k,j,N_TH+1)= -alpha_w*dble(CTHX(0,K,J,1)) + gamma_w*dble(CTHX(0,K,J,2)) +  TH_BAR(K,j,N_TH+1) !
       enddo
       enddo
       ENDIF


        k = 1
        do j = 1, np2
        do i = 1, np1
        var1(k) =  real(dble(CU1X(0,i-1,j-1)))
        k = k+1
        var1(k) =  real(dble(CU3X(0,i-1,j-1)))
        k = k+1
        var1(k) =  real(dble(CU2X(0,i-1,j-1)))
        k = k+1
        enddo
        enddo

       outputDIR='plane_data/'
       basename = 'data_parav'

       write(OutFileName,'(a,a,i5.5,a4)') trim(outputDIR),trim(basename)//"_n",kk,".vtk"

      	open(unit=13,file=OutFileName,access='stream',form='unformatted',status='unknown',&
	convert='big_endian',iostat=s1)


        WRITE(6,*) 'Wrting in ', OutFileName,' for paraview'
!HEADER: note termination with char(10)
	
	write(13) "# vtk DataFile Version 3.0"//char(10)
	write(13) trim(BaseName)//char(10)
	write(13) "BINARY"//char(10)
	write(13) "DATASET RECTILINEAR_GRID"//char(10)

	write(ss,fmt='(A10,3I5)') "DIMENSIONS",1,np1,np2
	write(13) ss//char(10)
!Xgrid
	write(ss,fmt='(A13,I6,A6)') "X_COORDINATES",1," float"
	write(13) char(10)//ss//char(10)
	write(13) 1.0
!Ygrid
	write(ss,fmt='(A13,I6,A6)') "Y_COORDINATES",np1," float"
	write(13) char(10)//ss//char(10)
	do k = 0, np1-1
	write(13) g1vtk(k)
	enddo
!Zgrid
	write(ss,fmt='(A13,I6,A6)') "Z_COORDINATES",np2," float"
	write(13) char(10)//ss//char(10)
	do j = 0, np2-1
	  write(13) g2vtk(j)
	enddo
        
!        write(ss,fmt='(A4,I4,A6)') "TIME",1, " double"
!        write(13) char(10)//ss//char(10)
!        write(13) real(TIME/(60.00*60.00))

!        write(ss,fmt='(A4,I4)') "TIME",TIME/(60.00*60.00)
        write(ss,fmt='(A10,I15)') "POINT_DATA",np1*np2
        write(13) char(10)//ss//char(10)
        write(13) "VECTORS velocity_vectors float 1"//char(10)
        write(13) var1
!	write(ss,fmt='(A10,I15)') "POINT_DATA",np1*np2
!	write(13) char(10)//ss//char(10)
!	write(13) "SCALARS U1 float 1"//char(10)
!	write(13) "LOOKUP_TABLE default"//char(10)
!	write(13) real(dble(CU3X(0,0:NZ+1,0:NY+1)))
!        write(13) "SCALARS U2 float 1"//char(10)
!	write(13) "LOOKUP_TABLE default"//char(10)
!	write(13) real(dble(CU2X(0,0:NZ+1,0:NY+1)))
!        write(13) "SCALARS U3 float 1"//char(10)
!        write(13) "LOOKUP_TABLE default"//char(10)
!        write(13) real(dble(CU1X(0,0:NZ+1,0:NY+1))) 
        write(13) "SCALARS Temp_d float 1"//char(10)
        write(13) "LOOKUP_TABLE default"//char(10)
        write(13) real(dble(CTHX(0,0:NZ+1,0:NY+1,1)))
        write(13) "SCALARS Salinity_d float 1"//char(10)
        write(13) "LOOKUP_TABLE default"//char(10)
        write(13) real(dble(CTHX(0,0:NZ+1,0:NY+1,2)))
        write(13) "SCALARS Temp_w float 1"//char(10)
        write(13) "LOOKUP_TABLE default"//char(10)
        write(13) real(th_wh(0:NZ+1,0:NY+1,1))
        write(13) "SCALARS Salinity_w float 1"//char(10)
        write(13) "LOOKUP_TABLE default"//char(10)
        write(13) real(th_wh(0:NZ+1,0:NY+1,2))
        write(13) "SCALARS Density_w float 1"//char(10)
        write(13) "LOOKUP_TABLE default"//char(10)
        write(13) real(th_wh(0:NZ+1,0:NY+1,3))
        write(13) "SCALARS Pressure float 1"//char(10)
        write(13) "LOOKUP_TABLE default"//char(10)
        write(13) real(dble(CPX(0,0:NZ+1,0:NY+1)))
        write(13) "SCALARS Urms float 1"//char(10)
        write(13) "LOOKUP_TABLE default"//char(10)
        write(13) real(wrms(0:NZ+1,0:NY+1))
        write(13) "SCALARS Wrms float 1"//char(10)
        write(13) "LOOKUP_TABLE default"//char(10)
        write(13) real(vrms(0:NZ+1,0:NY+1))
        write(13) "SCALARS time float 1"//char(10)
        write(13) "LOOKUP_TABLE default"//char(10)
        write(13) real(tim(0:NZ+1,0:NY+1))

!Close VTK File
	close(13)

deallocate(g1vtk,g2vtk,var1 )
!if (allocated(SPplane) ) deallocate(SPplane)
!if (allocated(DPplane) ) deallocate(DPplane)

return
end subroutine plane_parav_vel 



subroutine plot_para_tke(kk)

!      INCLUDE 'header_duct'
 use ntypes
 use Domain
 use Grid, only : dx,xpoint,ypoint
 use run_variable
 use variable_stat
! use les_chan_var
 implicit none
             
       integer i,j,k,np1,np2,kk,N
       !FILENAMES
       character(len=100) :: InFileName, OutFileName
       character(len=100):: OutFileName2,outputDIR2
       character(len=100) :: dataDIR, outputDIR
       character(len=80) :: basename,basename2
       !OTHER STRINGS
       character(len=25) :: ss


      real(r4),allocatable,dimension(:) :: g1vtk,g2vtk
      real(r4),allocatable,dimension(:) :: var_1
      integer :: s11


       np1 = NZ+2
       np2 = NY+2
       allocate (var_1(1:3*np1*np2))
       allocate(g1vtk(0:NZ+1), g2vtk(0:NY+1))


!       do j=0,NY+1
!         g1vtk(j) = xpoint(k,j)
!       enddo

       do j=0,NY+1
       do k=0,NZ+1
         g1vtk(k) = xpoint(k,j)
         g2vtk(j) = ypoint(k,j)
         DO N=1,N_TH
         th_wh(k,j,N)=dble(CTHX(0,k,j,N)) +  TH_BAR(K,j,N) !
         ENDDO
       enddo
       enddo 

       IF (Non_linear_ST) THEN
          call density_TC(.true.,.true.)
       ELSE 
       do j=0,NY+1
       do k=0,NZ+1
         th_wh(k,j,N_TH+1)= -alpha_w*dble(CTHX(0,K,J,1)) + gamma_w*dble(CTHX(0,K,J,2)) +  TH_BAR(K,j,N_TH+1) !
       enddo
       enddo
       ENDIF

 
        k = 1
        do j = 1, np2
        do i = 1, np1
        var_1(k) =  real(dble(CR1X(0,i-1,j-1)))
        k = k+1
        var_1(k) =  real(dble(CR3X(0,i-1,j-1)))
        k = k+1
        var_1(k) =  real(dble(CR2X(0,i-1,j-1)))
        k = k+1
        enddo
        enddo
 

       outputDIR='plane_tke/'
       basename = 'data_tke_parav'


       write(OutFileName,'(a,a,i5.5,a4)') trim(outputDIR),  &
      trim(basename)//"_n",kk,".vtk"

       open(unit=13,file=OutFileName,access='stream', &
      form='unformatted',status='unknown',            &
      convert='big_endian',iostat=s11)


        WRITE(6,*) 'Wrting in ', OutFileName,' for paraview'
!HEADER: note termination with char(10)
        
        write(13) "# vtk DataFile Version 3.0"//char(10)
        write(13) trim(BaseName)//char(10)
        write(13) "BINARY"//char(10)
        write(13) "DATASET RECTILINEAR_GRID"//char(10)

!        write(ss,fmt='(A4,2I4,A6)') "TIME",1,1, " double"
!        write(13) char(10)//ss//char(10)
!        write(13) real(TIME/(60.00*60.00))

        write(ss,fmt='(A10,3I5)') "DIMENSIONS",1,np1,np2
        write(13) ss//char(10)
        
!Xgrid
        write(ss,fmt='(A13,I6,A6)') "X_COORDINATES",1," float"
        write(13) char(10)//ss//char(10)
        write(13) 1.0
!Ygrid
        write(ss,fmt='(A13,I6,A6)') "Y_COORDINATES",np1," float"
        write(13) char(10)//ss//char(10)
        do k = 0,np1-1
        write(13) g1vtk(k)
        enddo
!Zgrid
        write(ss,fmt='(A13,I6,A6)') "Z_COORDINATES",np2," float"
        write(13) char(10)//ss//char(10)
        do j = 0,np2-1
          write(13) g2vtk(j)
        enddo

     
        
!Field
        write(ss,fmt='(A10,I15)') "POINT_DATA",np1*np2
        write(13) char(10)//ss//char(10)
        write(13) "VECTORS velocity_vectors float 1"//char(10)
        write(13) var_1
        write(13) "SCALARS tke float 1"//char(10)
        write(13) "LOOKUP_TABLE default"//char(10)
        write(13) real(tke_mean(0:NZ+1,0:NY+1))
        write(13) "SCALARS dtkedt float 1"//char(10)
        write(13) "LOOKUP_TABLE default"//char(10)
        write(13) real(tke_1(0:NZ+1,0:NY+1))
        write(13) "SCALARS Tur_advection float 1"//char(10)
        write(13) "LOOKUP_TABLE default"//char(10)
        write(13) real(tke_2(0:NZ+1,0:NY+1))
        write(13) "SCALARS tur_production float 1"//char(10)
        write(13) "LOOKUP_TABLE default"//char(10)
        write(13) real(tke_3(0:NZ+1,0:NY+1))
        write(13) "SCALARS viscous_dissipation float 1"//char(10)
        write(13) "LOOKUP_TABLE default"//char(10)
        write(13) real(tke_4(0:NZ+1,0:NY+1))
        write(13) "SCALARS Buoyancy_flux float 1"//char(10)
        write(13) "LOOKUP_TABLE default"//char(10)
        write(13) real(tke_5(0:NZ+1,0:NY+1,1))
        write(13) "SCALARS pressure_transport float 1"//char(10)
        write(13) "LOOKUP_TABLE default"//char(10)
        write(13) real(tke_5(0:NZ+1,0:NY+1,2))
        write(13) "SCALARS trubulent_transport float 1"//char(10)
        write(13) "LOOKUP_TABLE default"//char(10)
        write(13) real(tke_6_2(0:NZ+1,0:NY+1))
        write(13) "SCALARS urho float 1"//char(10)
        write(13) "LOOKUP_TABLE default"//char(10)
!        write(13) real(tke_7(0:NZ+1,0:NY+1))
         write(13) real(tke_5_1(0:NZ+1,0:NY+1))

!Close VTK File
        close(13)

        deallocate(g1vtk,g2vtk,var_1)
        return
        end subroutine plot_para_tke


    subroutine plane_parav_zx(kk)

!      INCLUDE 'header_duct'
 use ntypes
 use Domain
 use Grid, only : dx,GZF, GYF,xpoint
 use run_variable
 use variable_stat
 use mpi_var, only: rank
 implicit none

       integer i,j,k,np1,np3,kk,N,jpn1
       !FILENAMES
       parameter(np1 = NX, np3 = NZ)
       character(len=100) :: InFileName, OutFileName
       character(len=100) :: dataDIR, outputDIR
       character(len=80) :: basename
       !OTHER STRINGS
       character(len=25) :: ss

       real(r8),allocatable,dimension(:,:) :: U1_ZX,U2_ZX,U3_ZX,TH_ZX,K_E_ZX
       real(r4) :: g1vtk(1:np1),g3vtk(1:np3), var_2(1:3*np1*np3)

!       th_wh(1:np1,1:np2,N_TH+1)
!       real(r4),allocatable,dimension(:,:) :: var1
!       real(r8),allocatable,dimension(:,:,:) :: th_wh
       !STATUS VARIABLES
       integer :: s11


       allocate (U1_ZX(1:(NXP+1)*NP,1:NZ+2))
       allocate (U2_ZX(1:(NXP+1)*NP,1:NZ+2))
       allocate (U3_ZX(1:(NXP+1)*NP,1:NZ+2))
       allocate (TH_ZX(1:(NXP+1)*NP,1:NZ+2))
       allocate (K_E_ZX(1:(NXP+1)*NP,1:NZ+2))

       U1_ZX(:,:)=0.0d0
       U2_ZX(:,:)=0.0d0
       U3_ZX(:,:)=0.0d0
       TH_ZX(:,:)=0.0d0
      K_E_ZX(:,:)=0.0d0

       jpn1 = int(NY/2)
       CALL  MPI_COMBINE (U1X,U1_ZX,NXP+1,NZ+2)
       CALL  MPI_COMBINE (U2X,U2_ZX,NXP+1,NZ+2)
       CALL  MPI_COMBINE (U3X,U3_ZX,NXP+1,NZ+2)
       CALL  MPI_COMBINE (THX,TH_ZX,NXP+1,NZ+2)
        do k=1,np1
         g1vtk(k) = dx(1)*(k-1)
        enddo

       do j=1,np3
         g3vtk(j) = xpoint(j,1)
       enddo
       
       
       if (rank.eq.0)then
       outputDIR='plane_data/'
       basename = 'data_parav_zx'

        write(OutFileName,'(a,a,i4.4,a2,i5.5,a4)') trim(outputDIR), &
        trim(basename)//"_p",jpn1,"_n",kk,".vtk"

       write(OutFileName,'(a,a,i5.5,a4)') trim(outputDIR), &
       trim(basename)//"_n",kk,".vtk"  

       open(unit=13,file=OutFileName,access='stream',  &
       form='unformatted',status='unknown',             &
       convert='big_endian',iostat=s11)
        WRITE(6,*) 'Wrting in ', OutFileName,' for paraview for plane jpn1=',jpn1
!HEADER: note termination with char(10)

        write(13) "# vtk DataFile Version 3.0"//char(10)
       write(13) trim(BaseName)//char(10)
        write(13) "BINARY"//char(10)
        write(13) "DATASET RECTILINEAR_GRID"//char(10)

        write(ss,fmt='(A10,3I5)') "DIMENSIONS",np1,np3,1
        write(13) ss//char(10)
!Xgrid
        write(ss,fmt='(A13,I6,A6)') "X_COORDINATES",np1," float"
        write(13) char(10)//ss//char(10)
        do k = 1,np1
        write(13) g1vtk(k)
        enddo
!Ygrid
        write(ss,fmt='(A13,I6,A6)') "Y_COORDINATES",np3," float"
        write(13) char(10)//ss//char(10)
        do j = 1,np3
          write(13) g3vtk(j)
        enddo

!Zgrid
        write(ss,fmt='(A13,I6,A6)') "Z_COORDINATES",1," float"
        write(13) char(10)//ss//char(10)
        write(13) real(GYF(jpn1))


!Field
        j = 1
        do k = 1, np3
        do i = 1, np1
        var_2(j) =  real(U1_ZX(i,k))
        j = j+1
        var_2(j) =  real(U3_ZX(i,k))
        j = j+1
        var_2(j) =  real(U2_ZX(i,k))
        j = j+1
        K_E_ZX(i,k)=0.5d0*((real(U1_ZX(i,k)))**2+(real(U2_ZX(i,k)))**2+(real(U3_ZX(i,k)))**2)
        enddo
        enddo

!Field
        write(ss,fmt='(A10,I15)') "POINT_DATA",np1*np3
        write(13) char(10)//ss//char(10)
        write(13) "VECTORS velocity_vectors float 1"//char(10)
        write(13) var_2
        write(13) "SCALARS Temp_div float 1"//char(10)
        write(13) "LOOKUP_TABLE default"//char(10)
        write(13) real(TH_ZX(1:np1,1:np3))
        write(13) "SCALARS Kinetic_Energy float 1"//char(10)
        write(13) "LOOKUP_TABLE default"//char(10)
        write(13) real(K_E_ZX(1:np1,1:np3))

!Close VTK File
        close(13)
        endif

        deallocate(U1_ZX,U2_ZX,U3_ZX,TH_ZX,K_E_ZX)

        return
        end subroutine plane_parav_zx

