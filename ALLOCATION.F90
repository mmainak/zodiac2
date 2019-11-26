subroutine allocate_var
use ntypes
use Domain
use Grid
use Fft_var
use TIME_STEP_VAR
use run_variable
use mg_vari
use ADI_var
use mpi_var
use forcing
use IO,     only: I_OUT


implicit none

!Passed Variables
 integer        :: ok

!Local Variables
 integer           :: s_1
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! grid allocation
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

allocate (GX(0:NX+1), stat=s_1) 
allocate (GY(0:NY+1), stat=s_1)
allocate (GZ(0:NZ+1), stat=s_1)
allocate (DX(0:NX+1), stat=s_1)
allocate (DY(0:NY+1), stat=s_1)
allocate (DZ(0:NZ+1), stat=s_1)
allocate (GXF(0:NX), stat=s_1)
allocate (GYF(0:NY+1), stat=s_1) 
allocate (GZF(0:NZ), stat=s_1)
allocate (DXF(0:NX), stat=s_1)
allocate (DYF(0:NY), stat=s_1)
allocate (DZF(0:NZ), stat=s_1)

allocate (xpoint(0:NZ+1,0:NY+1), ypoint(0:NZ+1,0:NY+1))
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
if (s_1.NE.0) then
  write(I_OUT,*) "Error Allocating Grid Variables" 
  goto 1000
 endif


!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! fft allocation
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
allocate  ( KX (0:NX/3), stat=s_1)
allocate  ( KY(0:2*(NY/3)), stat=s_1)
allocate  ( KZ  (0:2*(NZ/3)), stat=s_1)
allocate  ( KX2 (0:NX/3), stat=s_1) 
allocate  ( KY2 (0:2*(NY/3)), stat=s_1)
allocate  ( KZ2 (0:2*(NZ/3)), stat=s_1)
allocate  (CIKX(0:NX/3), stat=s_1)
allocate  (CIKY(0:2*(NY/3)), stat=s_1)
allocate  ( CIKZ(0:2*(NZ/3)), stat=s_1)
allocate  ( CZX_PLANE(0:NZ,0:NX2P), stat=s_1)
allocate  ( CYZ_PLANE(0:NY,0:2*(NZ/3)), stat=s_1)

allocate  ( KXP(0:NX2P), stat=s_1)
allocate  ( KX2P(0:NX2P), stat=s_1)
allocate  (CIKXP(0:NX2P), stat=s_1)
if (s_1.NE.0) then
  write(I_OUT,*) "Error Allocating FFT Variables" 
  goto 1000
 endif
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! Input parameters and runtime variables
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

allocate (SPONGE_SIGMA(0:NZ+1,0:NY+1), stat=s_1)
allocate (SPONGE_SIGMA_OUT(0:NZ+1,0:NY+1), stat=s_1)
allocate (SPONGE_TEMP(0:NZ+1,0:NY+1), stat=s_1)

allocate (C_int(0:NY+1),C_int_le(0:NY+1), stat=s_1)


allocate (U_BC_LOWER_NWM(0:NXP,0:NZ+1), stat=s_1)
allocate (W_BC_LOWER_NWM(0:NXP,0:NZ+1), stat=s_1)
allocate (U_BC_UPPER_NWM(0:NXP,0:NZ+1), stat=s_1)
allocate (W_BC_UPPER_NWM(0:NXP,0:NZ+1), stat=s_1)

allocate (CU_BC_LOWER_NWM(0:NX2P,0:NZ+1), stat=s_1)
allocate (CU_BC_UPPER_NWM(0:NX2P,0:NZ+1), stat=s_1)
allocate (CW_BC_LOWER_NWM(0:NX2P,0:NZ+1), stat=s_1)
allocate (CW_BC_UPPER_NWM(0:NX2P,0:NZ+1), stat=s_1)

allocate (U1_bar(0:NY+2), stat=s_1)
allocate (U2_bar(0:NZV-1,0:NY+1), stat=s_1)
allocate (U3_bar(0:NZV-1,0:NY+1), stat=s_1)
allocate (TH_BAR(0:NZ+2,0:NY+2,1:N_TH+1), stat=s_1)
allocate (temp_mean(0:NZV-1,0:NY+1), stat=s_1) 

allocate (f(0:NZ+1,0:NY+1), stat=s_1)
allocate (g(0:NZ+1,0:NY+1), stat=s_1)
allocate (f_prime(0:NZ+1,0:NY+1),g_prime(0:NZ+1,0:NY+1), stat=s_1)

allocate (u_prime(0:NY+1),v_prime(0:NY+1),beam_th(0:NY+1),stat=s_1)

if (s_1.NE.0) then
  write(I_OUT,*) "Error Allocating VEL_BAR Variables" 
  goto 1000
endif

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! ! Global variables
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! allocate      (U1 (0:NX+1,0:NZ+1,0:NY+1), U2 (0:NX+1,0:NZ+1,0:NY+1), stat=s_1)
! allocate      (U3 (0:NX+1,0:NZ+1,0:NY+1), P  (0:NX+1,0:NZ+1,0:NY+1), stat=s_1)
! allocate      (R1 (0:NX+1,0:NZ+1,0:NY+1), R2 (0:NX+1,0:NZ+1,0:NY+1), stat=s_1)
! allocate      (R3 (0:NX+1,0:NZ+1,0:NY+1), F1 (0:NX+1,0:NZ+1,0:NY+1), stat=s_1)
! allocate      (F2 (0:NX+1,0:NZ+1,0:NY+1), F3 (0:NX+1,0:NZ+1,0:NY+1), stat=s_1)
! allocate      (S1 (0:NX+1,0:NZ+1,0:NY+1), U1b(0:NX+1,0:NZ+1,0:NY+1), stat=s_1)
! allocate      (U2b(0:NX+1,0:NZ+1,0:NY+1), U3b(0:NX+1,0:NZ+1,0:NY+1), stat=s_1)
! allocate      (S2 (0:NX+1,0:NZ+1,0:NY+1)) 
! allocate      (TH (0:NX+1,0:NZ+1,0:NY+1,1:N_TH), stat=s_1)
! allocate      (TH_BACK(0:NY+1,1:N_TH), stat=s_1)
! allocate      (FTH (0:NX+1,0:NZ+1,0:NY+1,1:N_TH), stat=s_1) 
! allocate      (RTH (0:NX+1,0:NZ+1,0:NY+1,1:N_TH), stat=s_1)
! if (s_1.NE.0) then
!   write(I_OUT,*) "Error Allocating VEL_GOL_REL Variables" 
!   goto 1000
! endif
! 
! allocate     (CU1(0:NX/2,0:NZ+1,0:NY+1), CU2(0:NX/2,0:NZ+1,0:NY+1), stat=s_1) 
! allocate     (CU3(0:NX/2,0:NZ+1,0:NY+1), CP (0:NX/2,0:NZ+1,0:NY+1), stat=s_1)
! allocate     (CR1(0:NX/2,0:NZ+1,0:NY+1), CR2(0:NX/2,0:NZ+1,0:NY+1), stat=s_1)
! allocate     (CR3(0:NX/2,0:NZ+1,0:NY+1), CF1(0:NX/2,0:NZ+1,0:NY+1), stat=s_1)
! allocate     (CF2(0:NX/2,0:NZ+1,0:NY+1), CF3(0:NX/2,0:NZ+1,0:NY+1), stat=s_1)
! allocate     (CS1(0:NX/2,0:NZ+1,0:NY+1), CU1b(0:NX/2,0:NZ+1,0:NY+1), stat=s_1)
! allocate     (CU2b(0:NX/2,0:NZ+1,0:NY+1), CU3b(0:NX/2,0:NZ+1,0:NY+1), stat=s_1) 
! allocate     (CS2(0:NX/2,0:NZ+1,0:NY+1), stat=s_1)
! allocate     (CTH(0:NX/2,0:NZ+1,0:NY+1,1:N_TH), stat=s_1)
! allocate     (CFTH(0:NX/2,0:NZ+1,0:NY+1,1:N_TH), stat=s_1)
! allocate     (CRTH(0:NX/2,0:NZ+1,0:NY+1,1:N_TH), stat=s_1)
! 
! if (s_1.NE.0) then
!   write(I_OUT,*) "Error Allocating VEL_GOL_COMPLX Variables" 
!   goto 1000
! endif

allocate      (TH_BACK(0:NY+1,1:N_TH), stat=s_1)
! EQUIVALENCE (U1,CU1), (U2,CU2), (U3,CU3), (R1,CR1), (R2,CR2) &
!         , (R3,CR3) , (F1,CF1), (F2,CF2), (F3,CF3), (P,CP), (S1,CS1) &
!         , (S2,CS2), (U1b,CU1b), (U2b,CU2b), (U3b,CU3b), (RTH,CRTH)  &
!         , (TH, CTH), (FTH, CFTH)

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! ! ADI variables
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

!allocate     (MATL (0:NX-1,0:NY+1), MATD(0:NX-1,0:NY+1), stat=s_1)
!allocate     (MATU(0:NX-1,0:NY+1), VEC(0:NX-1,0:NY+1), stat=s_1)
!allocate     (MATL_Z (0:NX-1,0:NZ+1), MATD_Z(0:NX-1,0:NZ+1), stat=s_1)
!allocate     (MATU_Z(0:NX-1,0:NZ+1), VEC_Z(0:NX-1,0:NZ+1), stat=s_1)

allocate     ( MATLX(0:NXP,0:NZ+1),MATDX(0:NXP,0:NZ+1), stat=s_1)
allocate     ( MATUX(0:NXP,0:NZ+1),VECX(0:NXP,0:NZ+1), stat=s_1)

allocate     ( MATLY(0:NXP,0:NY+1),MATDY(0:NXP,0:NY+1), stat=s_1)
allocate     ( MATUY(0:NXP,0:NY+1),VECY(0:NXP,0:NY+1), stat=s_1)

if (s_1.NE.0) then
  write(I_OUT,*) "Error Allocating ADI Variables" 
  goto 1000
endif

! JACOBIAN VAR
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
allocate     (INT_JACOB(0:NZ+2,0:NY+2), stat=s_1)
allocate     (GMAT_11(0:NZ+2,0:NY+2,2))
allocate     (GMAT_12(0:NZ+2,0:NY+2,2), GMAT_22(0:NZ+2,0:NY+2,2), stat=s_1)
allocate     (CJOB_11(0:NZ+2,0:NY+2,2), CJOB_12(0:NZ+2,0:NY+2,2), stat=s_1)
allocate     (CJOB_21(0:NZ+2,0:NY+2,2), CJOB_22(0:NZ+2,0:NY+2,2), stat=s_1)
if (s_1.NE.0) then
  write(I_OUT,*) "Error Allocating JACOBIAN Variables for theta" 
  goto 1000
endif
!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
! Multigrid variables - ALL ARRAYS 1D
!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
allocate      (V(1:(NY+2)*(NZ+2),0:NX2P), stat=s_1 )
allocate      (VC(1:NM,0:NX2P), stat=s_1)
allocate      (VB(1:(NY+2)*(NZ+2),0:NX2P), stat=s_1) 
allocate      (VBC(1:NM,0:NX2P), stat=s_1)
allocate      (RHS(1:(NY+2)*(NZ+2),0:NX2P), stat=s_1)
allocate      (RHSC(1:NM,0:NX2P), stat=s_1)
allocate      (A(1:(NY+2)*(NZ+2)*9,0:NX2P), stat=s_1)
allocate      (AC(1:9*NM,0:NX2P), stat=s_1)
allocate      (LDU(1:3*(NY+2)*(NZ+2),0:NX2P), stat=s_1)
allocate      (LDUC(1:3*NM,0:NX2P), stat=s_1)
allocate      (WORK(1:(NZ+2)*12,0:NX2P), stat=s_1)
allocate      (WA(1:(NY+2)*(NZ+2),0:NX2P), stat=s_1)
allocate      (WAC(NM,0:NX2P), stat=s_1)
allocate      (WB(1:(NY+2)*(NZ+2),0:NX2P), stat=s_1)
allocate      (WBC(1:NM,0:NX2P), stat=s_1)

!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
! MPI variables - ALL ARRAYS 1D
!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
!allocate       (OUT_CX(0:NX2P,0:NZV-1,0:NY+1),IN_CZ(0:NX2V-1,0:NZP,0:NY+1), stat=s_1)
!allocate       (IN_CX(0:NX2P,0:NZV-1,0:NY+1),OUT_CZ(0:NX2V-1,0:NZP,0:NY+1), stat=s_1)

allocate      (M_IN_C(1:NX2V*(NY+2)*NZV/NP),M_OUT_C(1:NX2V*(NY+2)*NZV/NP), stat=s_1)
allocate      (M_IN(1:NXV*(NY+2)*NZV/NP), M_OUT(1:NXV*(NY+2)*NZV/NP), stat=s_1)


if (s_1.NE.0) then
  write(I_OUT,*) "Error Allocating MG Variables for theta" 
  goto 1000
endif

write(I_OUT,'(a)') "ALLOCATION COMPLETED"
return
1000 continue 
 ok = s_1
 write(I_OUT,'(a)') "ALLOCATION FAILED"

return
end

subroutine allocation_u (FLAG)
use ntypes
use Domain
use run_variable
use IO,     only: I_OUT

implicit none
logical FLAG
integer           :: s_1

if (FLAG)then

if (.not. allocated(U1X) ) allocate( U1X(0:NXP,0:NZV-1,0:NY+1) )
if ( allocated(CU1X) ) then
deallocate( CU1X, stat=s_1)
endif
if (s_1.NE.0) then
  write(I_OUT,*) "Error deallocating cu1"
endif

else
if (.not. allocated(CU1X) ) allocate( CU1X(0:NX2P,0:NZV-1,0:NY+1))
if ( allocated(U1X) ) then
deallocate(U1X, stat=s_1)
endif
if (s_1.NE.0) then
  write(I_OUT,*) "Error deallocating u1"
endif

endif



return
end

subroutine allocation_v (FLAG)
use ntypes
use Domain
use run_variable
use IO,     only: I_OUT

implicit none
logical FLAG
integer           :: s_1


if (FLAG)then

if (.not. allocated(U2X) ) allocate( U2X(0:NXP,0:NZV-1,0:NY+1) )
if ( allocated(CU2X) ) then
deallocate( CU2X, stat=s_1)
endif
if (s_1.NE.0) then
  write(I_OUT,*) "Error deallocating cu1"
endif

else
if (.not. allocated(CU2X) ) allocate( CU2X(0:NX2P,0:NZV-1,0:NY+1))
if ( allocated(U2X) ) then
deallocate(U2X, stat=s_1)
endif
if (s_1.NE.0) then
  write(I_OUT,*) "Error deallocating u2"
endif

endif


return
end

subroutine allocation_w (FLAG)
use ntypes
use Domain
use run_variable
use IO,     only: I_OUT

implicit none
logical FLAG
integer           :: s_1

if (FLAG)then

if (.not. allocated(U3X) ) allocate( U3X(0:NXP,0:NZV-1,0:NY+1) )
if ( allocated(CU3X) ) then
deallocate( CU3X, stat=s_1)
endif
if (s_1.NE.0) then
  write(I_OUT,*) "Error deallocating cu1"
endif

else
if (.not. allocated(CU3X) ) allocate( CU3X(0:NX2P,0:NZV-1,0:NY+1))
if ( allocated(U3X) ) then
deallocate(U3X, stat=s_1)
endif
if (s_1.NE.0) then
  write(I_OUT,*) "Error deallocating u3"
endif

endif


return
end

subroutine allocation_p (FLAG)
use ntypes
use Domain
use run_variable
use IO,     only: I_OUT

implicit none
logical FLAG
integer           :: s_1

if (FLAG)then

if (.not. allocated(PX) ) allocate( PX(0:NXP,0:NZV-1,0:NY+1) )
if ( allocated(CPX) ) then
deallocate( CPX, stat=s_1)
endif
if (s_1.NE.0) then
  write(I_OUT,*) "Error deallocating CPX"
endif

else
if (.not. allocated(CPX) ) allocate( CPX(0:NX2P,0:NZV-1,0:NY+1))
if ( allocated(PX) ) then
deallocate(PX, stat=s_1)
endif
if (s_1.NE.0) then
  write(I_OUT,*) "Error deallocating PX"
endif

endif

return
end

subroutine allocation_th (FLAG)
use ntypes
use Domain
use run_variable
use IO,     only: I_OUT

implicit none
logical FLAG
integer           :: s_1,n

if (FLAG)then
 
if (.not. allocated(THX) ) allocate( THX(0:NXP,0:NZV-1,0:NY+1,1:N_TH) ) 
if ( allocated(CTHX ) ) then
deallocate( CTHX, stat=s_1)
endif
if (s_1.NE.0) then
  write(I_OUT,*) "Error deallocating CTHX"
endif
 
else
if (.not. allocated(CTHX ) ) allocate( CTHX(0:NX2P,0:NZV-1,0:NY+1,1:N_TH) ) 
if ( allocated(THX) ) then
deallocate(THX, stat=s_1)
endif
if (s_1.NE.0) then
  write(I_OUT,*) "Error deallocating THX"
endif
 
endif

return
end


subroutine allocation_ub (FLAG)
use ntypes
use Domain
use run_variable
use IO,     only: I_OUT

implicit none
logical FLAG
integer           :: s_1


if (FLAG)then

if (allocated(CU2bx) ) then
 deallocate( CU2bX )
 deallocate( CU3bX )
endif

if (.not. allocated(U2bX) ) then
 allocate( U2bX(0:NXP,0:NZV-1,0:NY+1) )
 allocate( U3bX(0:NXP,0:NZV-1,0:NY+1) )
endif

else

if (allocated(U2bX) ) then
 deallocate( U2bX )
 deallocate( U3bX )
endif

if (.not. allocated(CU2bX) ) then
 allocate( CU2bX(0:NX2P,0:NZV-1,0:NY+1) )
 allocate( CU3bX(0:NX2P,0:NZV-1,0:NY+1) )
endif

endif


return
end


subroutine allocation_R1(FLAG)
use ntypes
use Domain
use run_variable
use IO,     only: I_OUT

implicit none
logical FLAG
integer           :: s_1

if (FLAG)then

if (.not. allocated(R1X) ) allocate( R1X(0:NXP,0:NZV-1,0:NY+1) )
if ( allocated(CR1X) ) then
deallocate( CR1X, stat=s_1)
endif
if (s_1.NE.0) then
  write(I_OUT,*) "Error deallocating CR1X"
endif

else
if (.not. allocated(CR1X) ) allocate( CR1X(0:NX2P,0:NZV-1,0:NY+1))
if ( allocated(R1X) ) then
deallocate(R1X, stat=s_1)
endif
if (s_1.NE.0) then
  write(I_OUT,*) "Error deallocating R1"
endif

endif

return
end

subroutine allocation_R2(FLAG)
use ntypes
use Domain
use run_variable
use IO,     only: I_OUT

implicit none
logical FLAG
integer           :: s_1

if (FLAG)then

if (.not. allocated(R2X) ) allocate( R2X(0:NXP,0:NZV-1,0:NY+1) )
if ( allocated(CR2X) ) then
deallocate( CR2X, stat=s_1)
endif
if (s_1.NE.0) then
  write(I_OUT,*) "Error deallocating CR2X"
endif

else
if (.not. allocated(CR2X) ) allocate( CR2X(0:NX2P,0:NZV-1,0:NY+1))
if ( allocated(R2X) ) then
deallocate(R2X, stat=s_1)
endif
if (s_1.NE.0) then
  write(I_OUT,*) "Error deallocating R2"
endif

endif

return
end

subroutine allocation_R3(FLAG)
use ntypes
use Domain
use run_variable
use IO,     only: I_OUT
 
implicit none
logical FLAG
integer           :: s_1
 
if (FLAG)then
 
if (.not. allocated(R3X) ) allocate( R3X(0:NXP,0:NZV-1,0:NY+1) )
if ( allocated(CR3X) ) then
deallocate( CR3X, stat=s_1)
endif
if (s_1.NE.0) then
  write(I_OUT,*) "Error deallocating CR3X"
endif
 
else
if (.not. allocated(CR3X) ) allocate( CR3X(0:NX2P,0:NZV-1,0:NY+1))
if ( allocated(R3X) ) then
deallocate(R3X, stat=s_1)
endif
if (s_1.NE.0) then
  write(I_OUT,*) "Error deallocating R3"
endif
 
endif

return
end

subroutine allocation_Rth (FLAG)
use ntypes
use Domain
use run_variable
use IO,     only: I_OUT

implicit none
logical FLAG
integer           :: s_1,n

if (FLAG)then
 
if (.not. allocated(RTHX) ) allocate( RTHX(0:NXP,0:NZV-1,0:NY+1,1:N_TH) ) 
if ( allocated(CRTHX ) ) then
deallocate( CRTHX, stat=s_1)
endif
if (s_1.NE.0) then
  write(I_OUT,*) "Error deallocating CRTHX"
endif
 
else
if (.not. allocated(CRTHX ) ) allocate( CRTHX(0:NX2P,0:NZV-1,0:NY+1,1:N_TH) ) 
if ( allocated(RTHX) ) then
deallocate(RTHX, stat=s_1)
endif
if (s_1.NE.0) then
  write(I_OUT,*) "Error deallocating RTHX"
endif
 
endif

return
end


subroutine allocation_F1(FLAG)
use ntypes
use Domain
use run_variable
use IO,     only: I_OUT

implicit none
logical FLAG
integer           :: s_1

if (FLAG)then

if (.not. allocated(F1X) ) allocate( F1X(0:NXP,0:NZV-1,0:NY+1) )
if ( allocated(CF1X) ) then
deallocate( CF1X, stat=s_1)
endif
if (s_1.NE.0) then
  write(I_OUT,*) "Error deallocating CF1X"
endif

else
if (.not. allocated(CF1X) ) allocate( CF1X(0:NX2P,0:NZV-1,0:NY+1))
if ( allocated(F1X) ) then
deallocate(F1X, stat=s_1)
endif
if (s_1.NE.0) then
  write(I_OUT,*) "Error deallocating F1X"
endif

endif

return
end

subroutine allocation_F2(FLAG)
use ntypes
use Domain
use run_variable
use IO,     only: I_OUT

implicit none
logical FLAG
integer           :: s_1

if (FLAG)then

if (.not. allocated(F2X) ) allocate( F2X(0:NXP,0:NZV-1,0:NY+1) )
if ( allocated(CF2X) ) then
deallocate( CF2X, stat=s_1)
endif
if (s_1.NE.0) then
  write(I_OUT,*) "Error deallocating CF2X"
endif

else
if (.not. allocated(CF2X) ) allocate( CF2X(0:NX2P,0:NZV-1,0:NY+1))
if ( allocated(F2X) ) then
deallocate(F2X, stat=s_1)
endif
if (s_1.NE.0) then
  write(I_OUT,*) "Error deallocating F2X"
endif

endif

return
end

subroutine allocation_F3(FLAG)
use ntypes
use Domain 
use run_variable 
use IO,     only: I_OUT
 
implicit none 
logical FLAG 
integer           :: s_1 
 
if (FLAG)then
 
if (.not. allocated(F3X) ) allocate( F3X(0:NXP,0:NZV-1,0:NY+1) )
if ( allocated(CF3X) ) then
deallocate( CF3X, stat=s_1)
endif
if (s_1.NE.0) then
  write(I_OUT,*) "Error deallocating CF3X"
endif
 
else
if (.not. allocated(CF3X) ) allocate( CF3X(0:NX2P,0:NZV-1,0:NY+1))
if ( allocated(F3X) ) then
deallocate(F3X, stat=s_1)
endif
if (s_1.NE.0) then
  write(I_OUT,*) "Error deallocating F3X"
endif

endif 
 
return 
end

subroutine allocation_Fth (FLAG)
use ntypes
use Domain
use run_variable
use IO,     only: I_OUT

implicit none
logical FLAG
integer           :: s_1,n

if (FLAG)then
 
if (.not. allocated(FTHX) ) allocate( FTHX(0:NXP,0:NZV-1,0:NY+1,1:N_TH) ) 
if ( allocated(CFTHX ) ) then
deallocate( CFTHX, stat=s_1)
endif
if (s_1.NE.0) then
  write(I_OUT,*) "Error deallocating CRTHX"
endif
 
else
if (.not. allocated(CFTHX ) ) allocate( CFTHX(0:NX2P,0:NZV-1,0:NY+1,1:N_TH) ) 
if ( allocated(FTHX) ) then
deallocate(FTHX, stat=s_1)
endif
if (s_1.NE.0) then
  write(I_OUT,*) "Error deallocating RTHX"
endif
 
endif

return
end




subroutine allocate_temps

 use Domain
 use Grid
 use variable_stat

 use IO,     only: I_OUT
 
 implicit none
!Passed Variables
 integer           :: ok

!Local Variables
 integer                    :: s_1
 ok=0
 s_1=0

! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! Variables for outputting statistics
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
allocate     (UBAR(0:NY+1),VBAR(0:NY+1),WBAR(0:NY+1), stat=s_1)
allocate     (URMS(0:NZ+1,0:NY+1),VRMS(0:NZ+1,0:NY+1), stat=s_1)
allocate     (WRMS(0:NZ+1,0:NY+1), stat=s_1)
allocate     (UV(0:NZ+1,0:NY+1),UW(0:NZ+1,0:NY+1), stat=s_1)
allocate     (WV(0:NZ+1,0:NY+1),PV(0:NZ+1,0:NY+1), PU(0:NZ+1,0:NY+1), stat=s_1)
allocate     (DWDY(0:NZ+1,0:NY+1), DUDY(0:NZ+1,0:NY+1), stat=s_1)
allocate     (DWDZ(0:NZ+1,0:NY+1), DUDZ(0:NZ+1,0:NY+1), stat=s_1)

      
allocate     (SHEAR(0:NZ+1,0:NY+1), stat=s_1)
allocate     (OMEGA_X(0:NZ+1,0:NY+1),OMEGA_Y(0:NZ+1,0:NY+1), stat=s_1)
allocate     (OMEGA_Z(0:NZ+1,0:NY+1),TKE(0:NZ+1,0:NY+1), stat=s_1)

if (s_1.NE.0) then
  write(I_OUT,*) "Error Allocating statistics Variables" 
endif
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! Variables needed for SAVE_STATS_TH
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
allocate     (THBAR(0:NY+1,1:N_TH), PE_DISS(0:NY+1,1:N_TH), stat=s_1)
allocate     (THRMS(0:NZ+1,0:NY+1,1:N_TH),Rig(0:NZ+1,0:NY+1), stat=s_1)
allocate     (th_wh(0:NZ+1,0:NY+1,1:N_TH+1), stat=s_1)
allocate     (tim(0:NZ+1,0:NY+1), stat=s_1)
allocate     (THV(0:NZ+1,0:NY+1,1:N_TH), stat=s_1)
allocate     (THW(0:NZ+1,0:NY+1,1:N_TH), stat=s_1)
allocate     (DTHDY(0:NZ+1,0:NY+1,1:N_TH), stat=s_1)
allocate     (DTHDZ(0:NZ+1,0:NY+1,1:N_TH), stat=s_1)
allocate     (DELTA_N2(0:NZ+1,0:NY+1,1:N_TH), stat=s_1)
if (s_1.NE.0) then
  write(I_OUT,*) "Error Allocating statistics Variables for theta" 
endif
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&      
! Variables for tkebudget
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
allocate     (EPSILON(0:NZ+1,0:NY+1), stat=s_1)
allocate     (tke_mean(0:NZ+1,0:NY+1),tke_mean_old(0:NZ+1,0:NY+1), stat=s_1)
allocate     (energy_mean(0:NZ+1,0:NY+1),energy_mean_old(0:NZ+1,0:NY+1), stat=s_1)
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

! allocate tke variavles 

      allocate (tke_1(0:NZ+1,0:NY+1), stat=s_1)
      allocate (tke_2(0:NZ+1,0:NY+1), stat=s_1)
      allocate (tke_2_1(0:NZ+1,0:NY+1), stat=s_1)
      allocate (tke_2_2(0:NZ+1,0:NY+1), stat=s_1)
      allocate (tke_3(0:NZ+1,0:NY+1), stat=s_1)
      allocate (tke_3_1(0:NZ+1,0:NY+1), stat=s_1)
      allocate (tke_3_2(0:NZ+1,0:NY+1), stat=s_1)
      allocate (tke_3_3(0:NZ+1,0:NY+1), stat=s_1)
      allocate (tke_4(0:NZ+1,0:NY+1), stat=s_1)
      allocate (tke_5_1(0:NZ+1,0:NY+1), stat=s_1)
      allocate (tke_5(0:NZ+1,0:NY+1,1:N_TH), stat=s_1)
      allocate (tke_6_1(0:NZ+1,0:NY+1), stat=s_1)
      allocate (tke_6_1_1(0:NZ+1,0:NY+1), stat=s_1)
      allocate (tke_6_1_2(0:NZ+1,0:NY+1), stat=s_1)
      ALLOCATE (tke_6_1_3(0:NZ+1,0:NY+1), stat=s_1)
      allocate (tke_6_2(0:NZ+1,0:NY+1), stat=s_1)
      allocate (tke_6_2_1(0:NZ+1,0:NY+1), stat=s_1)
      allocate (tke_6_2_2(0:NZ+1,0:NY+1), stat=s_1)
      allocate (tke_6_2_3(0:NZ+1,0:NY+1), stat=s_1)
      allocate (tke_7(0:NZ+1,0:NY+1), stat=s_1)
      allocate (S1_mean(0:NZ+1,0:NY+1), stat=s_1)
      allocate (p_mean(0:NZ+1,0:NY+1), stat=s_1)
      allocate (transport(0:NZ+1,0:NY+1), stat=s_1)
if (s_1.NE.0) then
  write(I_OUT,*) "Error Allocating statistics Variables for tke" 
endif      


 if (s_1.NE.0) write(I_OUT,'(a,i5)') "ERROR ALLOCATING TEMPS: stat=",s_1
  ok=s_1

 return
end subroutine allocate_temps

subroutine deallocate_temps
 use IO,     only: I_OUT
use variable_stat


implicit none
!Passed Variables
 integer       :: ok

!Local Variables
 integer                   :: s_1
 s_1=0
 ok=0

! Deallocate Variables for outputting statistics
     deallocate (UBAR,VBAR,WBAR)
     deallocate (URMS,VRMS, &
            WRMS, UV,UW, WV,PU,PV, &
            DWDY, DUDY,DWDZ, DUDZ, &
            SHEAR, OMEGA_X,OMEGA_Y,OMEGA_Z,TKE)
! deallocate variables needed for SAVE_STATS_TH
     
      deallocate (THBAR, PE_DISS, Rig)
      deallocate (THRMS,THV, THW, DTHDY, DTHDZ)

! Variables for tkebudget
      deallocate (epsilon, DELTA_N2, &
              tke_mean,tke_mean_old,energy_mean,energy_mean_old)     

! allocate tke variavles
	DEALLOCATE (tke_1,tke_2,tke_2_1, &
            tke_2_2,tke_3,tke_3_1,tke_3_2,tke_3_3,tke_4,tke_5, &
            tke_6_1,tke_6_1_1,tke_6_1_2,tke_6_1_3,tke_6_2,     &
            tke_6_2_1,tke_6_2_2,tke_6_2_3,tke_7,S1_mean,p_mean, &
            transport,tke_5_1)



 if (s_1.NE.0) write(I_OUT,'(a31,i5)') "ERROR DEALLOCATING TEMPS: stat=",s_1
 ok=s_1

return
end subroutine deallocate_temps

subroutine allocation_les_var
use ntypes
use Domain
use run_variable, only : NU_T, KAPPA_T
use les_chan_var
use IO,     only: I_OUT


implicit none

!Passed Variables
 integer        :: ok

!Local Variables
 integer           :: s_1
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!  variables
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

allocate     (tke_sgs_t(0:NZ+1,0:NY+1), stat=s_1)
allocate     (tke_sgs_p(0:NZ+1,0:NY+1), stat=s_1)
allocate     (tke_sgs_diss(0:NZ+1,0:NY+1), stat=s_1)
allocate     (Sij_mean(0:NZ+1,0:NY+1,1:6), stat=s_1)
allocate     (TAU_mean(0:NZ+1,0:NY+1,1:6), stat=s_1)
allocate     (NU_T_mean(0:NZ+1,1:NY+1), stat=s_1)
allocate     (NU_T(0:NXP,0:NZ+1,0:NY+1),  stat=s_1)
allocate     (KAPPA_T(0:NXP,0:NZ+1,0:NY+1,1), stat=s_1)

allocate     (U_2BAR(0:NXP,0:NZV-1,0:NY+1,1:3), stat=s_1)
allocate     (U_BAR_TIL(0:NXP,0:NZV-1,0:NY+1,1:3), stat=s_1)
allocate     (U_4BAR(0:NXP,0:NZV-1,0:NY+1,1:3), stat=s_1)
allocate     (S_2BAR(0:NXP,0:NZV-1,0:NY+1), stat=s_1)
allocate     (U1_bar_les(0:NZ+1,0:NY+1), stat=s_1)
allocate     (U2_bar_les(0:NZ+1,0:NY+1), stat=s_1)
allocate     (U3_bar_les(0:NZ+1,0:NY+1), stat=s_1)
allocate     (DELTA_Y(0:NZ+1,0:NY+1), stat=s_1)
allocate     (C_DYN(0:NZ+1,0:NY+1), stat=s_1)

allocate     (C_DYN_TH(0:NZ+1,0:NY+1,1:N_TH), stat=s_1)
allocate     (KAPPA_T_mean(0:NZ+1,0:NY+1,1:N_TH), stat=s_1)

allocate     (GMAT_11_z(0:NZ+2,0:NY+2,1:2), stat=s_1)
allocate     (GMAT_11_y(0:NZ+2,0:NY+2,1:2), stat=s_1)
allocate     (GMAT_22_z(0:NZ+2,0:NY+2,1:2), stat=s_1)
allocate     (GMAT_22_y(0:NZ+2,0:NY+2,1:2), stat=s_1)
allocate     (GMAT_12_z(0:NZ+2,0:NY+2,1:2), stat=s_1)
allocate     (GMAT_12_y(0:NZ+2,0:NY+2,1:2), stat=s_1)

allocate     (ST_rij(0:NXP,0:NZ+1,0:NY+1,1:6), stat=s_1)


if (s_1.NE.0) then
  write(I_OUT,*) "Error Allocating LES Variables "
  goto 1000
endif

write(I_OUT,'(a)') "ALLOCATION FOR LES COMPLETED"
return
1000 continue
 ok = s_1
 write(I_OUT,'(a)') "ALLOCATION LES VAR FAILED"

return
end


subroutine allocate_les_tmp

 use Domain
 use Grid
 use les_chan_var

 use IO,     only: I_OUT

 implicit none
!Passed Variables
 integer           :: ok

!Local Variables
 integer                    :: s_1
 ok=0
 s_1=0

! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! Variables for les  c_dyn
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
allocate (numerator(0:NXP,0:NZ+1,0:NY+1), stat=s_1 )
allocate (denominator(0:NXP,0:NZ+1,0:NY+1), stat=s_1 )
allocate (denominator_sum(0:NZ+1,0:NY+1), stat=s_1 )
allocate (numerator_sum(0:NZ+1,0:NY+1), stat=s_1 )

end subroutine allocate_les_tmp


subroutine deallocate_les_tmp

 use Domain
 use Grid
 use les_chan_var

 use IO,     only: I_OUT

 implicit none
!Passed Variables
 integer           :: ok

!Local Variables
 integer                    :: s_1
 ok=0
 s_1=0

! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! Variables for les  c_dyn
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
deallocate (numerator, stat=s_1 )
deallocate (denominator, stat=s_1 )
deallocate (denominator_sum, stat=s_1 )
deallocate (numerator_sum, stat=s_1 )

end subroutine deallocate_les_tmp
