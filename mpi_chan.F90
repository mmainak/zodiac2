      SUBROUTINE GHOST_CHAN_MPI_working
! This subroutine is part of the MPI package for the channel flow
! Diablo package.
! Here, we define a set of ghost cells on each process
! the ghost cells contain information from the neighboring nodes
! and allow us to compute finite differences over the local gridpoints.
! We need to update the contents of the ghost cells at the start of
! each Runge-Kutta substep




      RETURN
      END


      SUBROUTINE GHOST_CHAN_MPI

      RETURN
      END


      SUBROUTINE GHOST_ZERO_CHAN_MPI_previous

      RETURN
      END



      SUBROUTINE GHOST_ZERO_CHAN_MPI


      RETURN
      END

      SUBROUTINE GHOST_GRID_MPI
! This subroutine is part of the MPI package for the channel flow
! Diablo package.
! Here, we define a set of ghost cells on each process
! the ghost cells contain information from the neighboring nodes
! and allow us to compute finite differences over the local gridpoints.
! We need to update the contents of the ghost cells at the start of
! each Runge-Kutta substep

      RETURN
      END

! C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE THOMAS_FORWARD_REAL_MPI(A,B,C,G,NY,NX)


      RETURN
      END



! C----*|--.---------.---------.---------.---------.---------.---------.-|------
      SUBROUTINE THOMAS_FORWARD_COMPLEX_MPI(A,B,C,G,NY,NX)
! C----*|--.---------.---------.---------.---------.---------.---------.-|-----|

      RETURN
      END




! C----*|--.---------.---------.---------.---------.---------.---------.-|------
      SUBROUTINE THOMAS_BACKWARD_REAL_MPI(A,B,C,G,NY,NX)
! C----*|--.---------.---------.---------.---------.---------.---------.-|------


      RETURN
      END

! C----*|--.---------.---------.---------.---------.---------.---------.-|------
      SUBROUTINE THOMAS_BACKWARD_COMPLEX_MPI(A,B,C,G,NY,NX)
! C----*|--.---------.---------.---------.---------.---------.---------.-|-----|


      RETURN
      END


! C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE THOMAS_FORWARD_REAL_MPI_old(A,B,C,G,NY,NX)


      RETURN
      END



! C----*|--.---------.---------.---------.---------.---------.---------.-|------
      SUBROUTINE THOMAS_FORWARD_COMPLEX_MPI_old(A,B,C,G,NY,NX)


      RETURN
      END




! C----*|--.---------.---------.---------.---------.---------.---------.-|------
      SUBROUTINE THOMAS_BACKWARD_REAL_MPI_old(A,B,C,G,NY,NX)

      RETURN
      END

! C----*|--.---------.---------.---------.---------.---------.---------.-|------
      SUBROUTINE THOMAS_BACKWARD_COMPLEX_MPI_old(A,B,C,G,NY,NX)
! C----*|--.---------.---------.---------.---------.---------.---------.-|-----|


      RETURN
      END


!----*|--.---------.---------.---------.---------.---------.---------.-|------
      SUBROUTINE INIT_MPI
! C----*|--.---------.---------.---------.---------.---------.---------.-|-----|


      RETURN
      END


!----*|--.---------.---------.---------.---------.---------.---------.-|------
      SUBROUTINE INIT_CHAN_MPI


      RETURN
      END

      SUBROUTINE APPLY_BC_TH_MPI(MATL,MATD,MATU,VEC,N)

      
      RETURN
      END


      SUBROUTINE APPLY_BC_U2_MPI(MATL,MATD,MATU,VEC,K)

      RETURN
      END


      SUBROUTINE APPLY_BC_U1_MPI(MATL,MATD,MATU,VEC,k)

      RETURN
      END

      SUBROUTINE APPLY_BC_U1_MPI_C(MATL_C,MATD_C,MATU_C,VEC_C)

      RETURN
      END


      SUBROUTINE APPLY_BC_U3_MPI(MATL,MATD,MATU,VEC,k)

      RETURN
      END

      SUBROUTINE APPLY_BC_U3_MPI_C(MATL_C,MATD_C,MATU_C,VEC_C)

      RETURN
      END


      SUBROUTINE APPLY_BC_REM_DIV_MPI(MATL_C,MATD_C,MATU_C,VEC_C,K)

        RETURN
        END

      SUBROUTINE APPLY_BC_POISSON_MPI(MATL_C,MATD_C,MATU_C,VEC_C,K)
! This subroutine applies the boundary conditions for the Poisson Eq.
! Note, MATL, MATD, etc. are dimensioned in header


      RETURN
      END

      SUBROUTINE APPLY_WALL_MODEL_MPI
! This subroutine applies the boundary conditions for the Poisson Eq.
! Note, MATL, MATD, etc. are dimensioned in header
      INCLUDE 'header_duct'
      INCLUDE "mpif.h"
      INCLUDE 'header_mpi'


      RETURN
      END

      SUBROUTINE APPLY_BC_VEL_MPI
! This subroutine applies the boundary conditions for the Poisson Eq.
! Note, MATL, MATD, etc. are dimensioned in header
      INCLUDE 'header_duct'
      INCLUDE "mpif.h"
      INCLUDE 'header_mpi'
      INTEGER  J1i,J2e  

 
      RETURN
      END
 


      SUBROUTINE APPLY_BC_ANWM_MPI
! This subroutine applies the boundary conditions for the Poisson Eq.
! Note, MATL, MATD, etc. are dimensioned in header


       RETURN
       END 


      SUBROUTINE APPLY_NUT_NWM_MPI_UPPER
 
       RETURN
       END



      SUBROUTINE APPLY_NUT_NWM_MPI_LOWER


       RETURN
       END

! C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE THOMAS_FORWARD_REAL_MPI_BACK(A,B,C,G,NY,NX)


      RETURN
      END


! C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE THOMAS_BACKWARD_REAL_MPI_BACK(A,B,C,G,NY,NX)


      RETURN
      END





