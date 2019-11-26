      subroutine MULTIGRID(P1,RHS1,xind)
use ntypes
use Domain
use Grid
use Fft_var
use TIME_STEP_VAR
use run_variable
use mg_vari 
implicit none
! c-----------------------------------------------------------------------
! c     this is an example of a main program using mgd9v
! c-----------------------------------------------------------------------
! c
!       integer xind,i,j,k
! 
! c
! c     parameter( nm=   774, nxf= 65, nyf= 33 )
! c     parameter( nm=  2919, nxf=129, nyf= 65 )
! c      parameter( nm= 11304, nxf=257, nyf=129 )
! c
! c      double precision
! c     +      v(nxf*nyf), vc(nm), vb(nxf*nyf), vbc(nm),val(nxf,nyf),
! c     +      rhs(nxf*nyf), rhsc(nm), a(nxf*nyf*9), ac(9*nm),
! c     +      ldu(3*nxf*nyf), lduc(3*nm), work(nxf*12),
! c     +      wa(nxf*nyf), wac(nm), wb(nxf*nyf), wbc(nm), tol, resno,
! c     +      gyf(0:nyf+1),gzf(0:nxf+1),hx,hy,x,y,pi

      integer xind,i,j,k
      double precision  &
             hx,hy,x,y     
      real*8 P1((NY+2)*(NZ+2)),RHS1((NY+2)*(NZ+2)),bv(4)
!      real*8 xpoint(NZ+2,NY+2), ypoint(NZ+2,NY+2)

! c     user data statements,
! c
! c     data nxc,nyc,levels/5,3,5/
! c     data nxc,nyc,levels/5,3,6/
! c      data nxc,nyc,levels/5,3,7/
! c
! c      data maxit,istart,iprep/30,1,0/
! c      data iout/1,0,0,2,1,0/
 !     data iout/0,0,0,0,0,0/
         iout(:) = 0
! c     
!       data iprep/0/
! c      data nout/6/
! c
! c     open(unit=nout,file='output')
! c
! c     problem set up
! c
! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! c     boundary condition      
      bc(1) = 2
      bc(2) = 1
      bc(3) = 1
      bc(4) = 2
      
      bv(:) = 0.d0
! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCcc

      v(:,xind)   = P1
      rhs(:,xind) = RHS1


      IF ( INIT_FLAG .EQ. .TRUE. ) THEN

      
! c      write(nout,1)
! c    1 format('1')
      IF (IC_TYPE .eq. 7) THEN
       call siam_curv( NZ+2,NY+2,bc,KX2P(xind),bv,a(:,xind), &
       GMAT_11,GMAT_22, GMAT_12,INT_JACOB)
      ELSE
       call siam
      ENDIF  
! c      
! c
! c     solution of the linear system
! c
! c
      iprep = 0
      tol=1.0d-5
      ifail=111
      istart=1
      maxit=0
! c
! c      open(40,file='flow_field.txt',form='formatted',status='unknown')
! 
! c      do  j = 1, NY
! c        y=(j-1)*hy
! c        write(40,*)gyf(j)
! c      enddo
! c      do  i = 1, NZ
! c        x=(i-1)*hx
! c        write(40,*)gzf(i)
! c      enddo
! 
! c      call bc_con(NZ,NY,bc,bv,rhs(:,xind))
!      
      
      call mgd9v(levels, nzc, nyc, NZ+2, NY+2, nm, &
                iout, istart, iprep, maxit, tol,   &
                rhs(:,xind), rhsc(:,xind), a(:,xind), &
                ac(:,xind), v(:,xind), vc(:,xind),    &
                vb(:,xind), vbc(:,xind), work(:,xind), &
                ldu(:,xind), lduc(:,xind), wa(:,xind),  &
                wac(:,xind), wb(:,xind), wbc(:,xind),  &
                resno, ndid, ifail)

    
! c      write(6,*) "*MG grids have been initialized FOR*kx=",XIND,
! c     +    KX2(xind)

      ELSE
      tol=1.0d-8
      ifail=1
      iprep=1
      istart=0
      maxit=200

!       
! c      do j=1,NY
! c       do i=1,NZ
! c         k = (j-1)*NZ + i
! c         v(k,xind)= 0.0 !(-1)**k*0.5
! c          v(k,xind)=dsin(2.d0*pi*gzf(i))*dsin(2.d0*pi*gyf(j))
! c          v(k) = (gzf(i)-0.5)*(gzf(i)+0.5)*(gyf(j)-0.5)*(gyf(j)+0.5)
! c          write(40,*)rhs(k,xind)
! c       enddo
! c      enddo
!       

      call mgd9v(levels, nzc, nyc, NZ+2, NY+2, nm,       &
                iout, istart, iprep, maxit, tol,        &
                rhs(:,xind), rhsc(:,xind), a(:,xind),   &
                ac(:,xind), v(:,xind), vc(:,xind),      &
                vb(:,xind), vbc(:,xind), work(:,xind),  &
                ldu(:,xind), lduc(:,xind), wa(:,xind),  &
                wac(:,xind), wb(:,xind), wbc(:,xind),    &
                resno, ndid, ifail)

      P1 = v(:,xind)

! c      call SOR (v(:,xind),rhs(:,xind))
!       
! c      do j=1,NY+2
! c       do i=1,NZ+2
! c        k = (j-1)*(NZ+2) + i
! c        write(40,*)v(k,xind)
! c       enddo
! c      enddo
! c      close(40)
!       
!       
! c      open(21,file='GRID_XZ.dat',form='formatted',status='old') 
! 
! c      DO J=1,NY+2
! c       DO K=1,NZ+2
! c         READ(21,*)xpoint(K,J),ypoint(K,J)
! c       ENDDO
! c      ENDDO
! c      close(21) 
!       
! c111   format(3f12.8)
! 
! c      open(10,file='sample.dat',form='formatted',status='unknown')
! c      write(10,*) 'title = "sample mesh" '
! c      write(10,*) 'variables = "x", "y" "A" '
! c      write(10,*) 'zone i=',NZ+2 ,', j=', NY+2,'DATAPACKING=POINT'
!      
! c      do j=1,NY+2
! c       do i=1,NZ+2
! c         k = (j-1)*(NZ+2) + i
! c         write(10,111) xpoint(i,j),ypoint(i,j), v(k,xind)  
! c       enddo
! c      enddo
! c      close(10)
! c      stop

      ENDIF


      return
      end

      subroutine siam
! c-----------------------------------------------------------------------
      use mpi_var , only : RANK
      implicit none

      return
      end

      subroutine siam_curv( nx,ny,bc,KX2,bv,l,GMAT_11,GMAT_22,GMAT_12, &
                INT_JACOB)
! c-----------------------------------------------------------------------
use mpi_var, only:rank        
      implicit none
! c
! c     testproblem to be found in
! c     p.m. de zeeuw and e.j. van asselt
! c     the convergence rate of multi-level algorithms applied to the
! c     convection-diffusion equation
! c     siam j. sci. stat. comput. 6(2) april 1985, p499.
! c
! c     ( u  + u  ) + a u = f(x,y)
! c        xx   yy
! c
! c     version dd911009
! c-----------------------------------------------------------------------
      integer nx, ny, nout
      integer bc(4)
      double precision &
          eps, l(nx,ny,9), f(nx,ny),bv(4),pi 
      integer i, j, k,nxm
      double precision  &
          x, y, hx, hy, mu, mux, muy, a, axy, ah, b, bxy, bh, dir, r, &
          KX2, GMAT_11(NX+1,NY+1,2),GMAT_12(NX+1,NY+1,2),             &
          GMAT_22(NX+1,NY+1,2),INT_JACOB(NX+1,NY+1)
      parameter ( pi=3.141592654) 



      call zeros(nx*ny*9,l)
      call zeros(nx*ny  ,f)
! 
      do 20 j = 1, ny-1
        do 10 i = 1, nx-1
          l(i+1,j  ,4)= 0.25*GMAT_12(i+1,j+1,1)-0.25*GMAT_12(i+1,j,1) &
                         - GMAT_11(i+1,j,2)	  
          l(i+1,j+1,4)= 0.25*GMAT_12(i+1,j+2,1)-0.25*GMAT_12(i+1,j+1,1) &
                         - GMAT_11(i+1,j+1,2)
          l(i  ,j  ,6)= -0.25*GMAT_12(i,j+1,1)+0.25*GMAT_12(i,j,1)    &
                         - GMAT_11(i+1,j,2)
          l(i  ,j+1,6)= -0.25*GMAT_12(i,j+2,1)+0.25*GMAT_12(i,j+1,1)  &
                         - GMAT_11(i+1,j+1,2)

        
          l(i  ,j  ,8)= -0.25*GMAT_12(i+1,j,2)+0.25*GMAT_12(i,j,2)  &
                         - GMAT_22(i,j+1,1)
          l(i+1,j  ,8)= -0.25*GMAT_12(i+2,j,2)+0.25*GMAT_12(i+1,j,2) &
                         - GMAT_22(i+1,j+1,1)
          l(i  ,j+1,2)= +0.25*GMAT_12(i+1,j+1,2)-0.25*GMAT_12(i,j+1,2) &
                         - GMAT_22(i,j+1,1)
          l(i+1,j+1,2)= +0.25*GMAT_12(i+2,j+1,2)-0.25*GMAT_12(i+1,j+1,2) &
                         - GMAT_22(i+1,j+1,1)
     
         	  
          l(i+1,j+1,1)= -0.25*GMAT_12(i+1,j+1,2)-0.25*GMAT_12(i+1,j+1,1)
                       
	  l(i,j+1,3)  =  +0.25*GMAT_12(i+1,j+1,2)+0.25*GMAT_12(i,j+1,1)
	  
          l(i+1,j,7)  =  +0.25*GMAT_12(i+1,j,2)+0.25*GMAT_12(i+1,j+1,1)
         
	  l(i,j,9)    = -0.25*GMAT_12(i+1,j,2) -0.25*GMAT_12(i,j+1,1) 

   10   continue
   20 continue

      do 40 j = 1, ny
        do 30 i = 1, nx
           l(i,j,5) = GMAT_11(i+1,j,2)+GMAT_11(i,j,2) &
               +   GMAT_22(i,j+1,1)+GMAT_22(i,j,1)	   & 
               + KX2*INT_JACOB(i,j)
   30   continue
   40 continue

      
       
      IF (KX2 .LE. 0.d0 ) THEN
      write(6,*)'#################################################'
      write(6,*)'boundary condition'
      write(6,*) 'Bottom wall', bc(1), 'Right wall', bc(2),'Top wall', &
                bc(3),'Left wall', bc(4)
      write(6,*)'#################################################'
      write(6,*)'boundary value', 'Bottom wall', bv(1), 'Right wall', &
                bv(2),'Top wall',  bv(3),'Left wall', bv(4)
      ENDIF
      
      



! c     Boundary conditions
       do j = 1,ny
! c      Dirichlet/periodic left face
        if( bc(4).eq.1 ) then
         l(1,j,:)   =  0.0d0
         l(1,j,5)   =  1.d0*10000000
	 l(1,j,6)   =  1.d0*10000000
         f(1,j)     =  2.0*bv(4)*10000000
        else
         l(1,j,:)   =  0.0d0
         l(1,j,5)   = 1.0*1000000000
         l(1,j,6)   =-1.0*1000000000
         f(1,j)     = -bv(4)*1000000000
        endif

! c      Dirichlet/periodic right face
        if( bc(2).eq.1 ) then
         l(nx,j,:) =  0.0d0
         l(nx,j,5) =  1.0d0*1000000000
	 l(nx,j,4) =  1.0d0*1000000000
         f(nx,j)   =  2.0*bv(2)*1000000000
        else
         l(nx,j,:) =  0.0d0
         l(nx,j,5)   = 1.0*1000000000
         l(nx,j,4)   = -1.0*1000000000
         f(nx,j)     = bv(2)*1000000000
        endif
       enddo

       do i = 1,nx
! C      Dirichlet/periodic top face
        if( bc(3).eq.1 ) then
         l(i,ny,:) =  0.0d0
         l(i,ny,5) =  1.0d0*10000000
	 l(i,ny,2) =  1.0d0*10000000
         f(i,ny)   =  2.d0*bv(3)*10000000
        else
         l(i,ny,:) =  0.0d0
         l(i,ny,5) =  1.0*1000000000
         l(i,ny,2) = -1.0*1000000000
         f(i,ny)   =  bv(3)*1000000000  
        endif   
       enddo
       do i = 1,nx
! c      Dirichlet/periodic bottom face
        if( bc(1).eq.1 ) then
         l(i,1,:)   =  0.0d0
         l(i,1,5)   =  1.0d0*10000000
	 l(i,1,8)   =  1.0d0*10000000
         f(i,1)     =  2.0*bv(1)*10000000
        else
         l(i,1,:)   =  0.0d0
         l(i,1,5)   =  1.0*1000000000
         l(i,1,8)   = -1.0*1000000000
         f(i,1)     = -bv(1)*1000000000
        endif
       enddo

       if ((bc(1)==2) .and. (bc(2)==2).and. &
          (bc(3)==2) .and. (bc(4)==2) )then
        l(nx,ny,:) =  0.0d0
        l(nx,ny,5) =  1.0d0*10000000
! c        l(nx,ny,2) =  1.0d0*10000000

! c        l(nx,ny,5) =  1.0d0*10000000
! c        l(nx,ny,4) =  1.0d0*10000000
        f(nx,ny)   =  2.d0*bv(3)*10000000 
        write(6,*)'one point is zero'
       endif 


      return
      end
      
      
      subroutine bc_con(nx,ny,bc,bv,f)
      implicit none

      integer i,j,k,nx, ny,bc(4)
      real*8 f(nx,ny),bv(4)

! c      do j = 1,ny
! c       do i = 1,nx
! c       f(i,j)=1
! c       enddo
! c      enddo
! 
!      
! c     Boundary conditions
       do j = 1,ny
! c      Dirichlet/periodic left face
        if( bc(4).eq.1 ) then
         f(1,j)     =  bv(4)
        else
        endif

! c      Dirichlet/periodic right face
        if( bc(2).eq.1 ) then
         f(nx,j)   =  bv(2)
        else
        endif
       enddo

       do i = 1,nx
! c      Dirichlet/periodic bottom face
        if( bc(1).eq.1 ) then
         f(i,1)     =  bv(1)
        else
        endif

! c      Dirichlet/periodic top face
        if( bc(3).eq.1 ) then
         f(i,ny)   =  bv(3)
        else
        endif
       enddo 


      return
      end
 
      
      
      double precision function mu(eps,ah)
      double precision eps, ah

           if( ah.gt.eps )    then
             mu=    eps/(2.0d0*ah)
      else if( ah.lt.(-eps) ) then
             mu=1.0d0+eps/(2.0d0*ah)
      else
             mu=0.5d0
      end if
      return
      end
