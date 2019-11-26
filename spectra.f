      Program main
      INCLUDE 'header_fft'

      integer ni,nj,nk,nk_st,n_t, st,en,tstep
      integer nip,njp,NXP,NXV,NXP_L,NX_PAR,rank,NP
      parameter (ni=NX, nj=NY+2,nk_st=176,nk=245,NP=6)
      parameter (n_t = nk-nk_st-1)
      PARAMETER (NXV=ceiling(real(NX+2)/real(NP))*NP)
      PARAMETER (NXP=NXV/NP-1)
      Real*8 x(ni,nj),y(ni,nj)
      Real*8 u(ni,nj),v(ni,nj),w(ni,nj),p(ni,nj),th(ni,nj),
     &       th_wh(ni,nj)

      real*8   RI,U_0,omega,rho_0,nu,dt,dx,xloc
      PARAMETER( RI = 14.93d0, U_0 = 0.0395d0, omega = 1.d0,
     &           rho_0 = 1.d0, nu=1.0d-7 ) 
      integer  i,j,imax,jmax,kmax,k,kk,i_md,ind_lt,ind_rt
      integer  zindex,zindex2,xindex,xindex2,index_Mfac
      parameter (zindex=201,zindex2=240,xindex=341,xindex2=107)
      
      REAL*8     U_data(0:n_t+1,0:NZ+1,0:NY+1), t(0:n_t+1)	!(nk_st:nk)
      COMPLEX*16 CU_data(0:n_t/2,0:NZ+1,0:NY+1)
      EQUIVALENCE (U_data,CU_data)

      integer  debug,ier,itot
      integer  tecini,tecdat,teczne,tecnod,tecfil,tecend
      integer  visdouble,disdouble
      character*1 nulchar

	  CHARACTER*4 PID
	  CHARACTER*4 file_num
      CHARACTER*34 file_name_1,file_time
      CHARACTER*6 base_name
 	  character(len=18) :: title

      LOGICAL VEL_FIELD
!       PARAMETER (NKXV=ceiling(real(NKX+1)/real(NP))*NP) 
!       PARAMETER (NKXP=NKXV/NP-1)
!       PARAMETER (NX2V=ceiling(real(NXV/2)/real(NP))*NP)
!       PARAMETER (NX2P=NX2V/NP-1)    

      base_name = 'span2_'
	  print *, 'folder name=', base_name
      
      PI  = ATAN(1.0)*4.0

      write(6,*) 'U_0', U_0, 'RI', RI

      VEL_FIELD = .TRUE.
      

      nulchar = char(0)
      debug   = 0
      visdouble = 0
      disdouble = 1
      imax = ni
      jmax = nj
      kmax = 1        
      kk = 0


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       file_time = trim(folder_name)//'/time_bulk.txt'
!       open(20,file=file_time,status='old',form='formatted')

!       do k=0,nk_st-1
!       	read(20,*) time, dt,u_bulk,v_bulk,u_bc,v_bc
!       enddo 
      
!       time_st= time  

      
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!       file_time = trim(folder_name)//
!      &        '/fft_data.dat'
!       open(34,file=file_time,status='unknown',form='formatted')
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    
!       open(11,file='line_evol_w.dat',status='unknown',form='formatted')
!       write(11,*) 'title = "Grid" '
!       write(11,*) 'variables = "time", "z", "p", "u", "w", "pu", "pw" '
!       write(11,*) 'zone i=',nj ,', j=', nk-nk_st+1,'DATAPACKING=POINT'
       
      kk = 0
       
	DO k=nk_st,nk
	  kk= kk+1
		
      file_num = CHAR(MOD(k,10000)/1000+48) 
     &		//CHAR(MOD(k,1000)/100+48)     
     &		//CHAR(MOD(k,100)/10+48)       
     &		//CHAR(MOD(k,10)+48)
       
     	st=1
     
	 DO rank=0,NP-1
	   I=rank+1
       PID = '_'//CHAR(MOD(I,1000)/100+48)   
     &		//CHAR(MOD(I,100)/10+48)    
     &		//CHAR(MOD(I,10)+48)
       
       file_name_1 = 'yz_plane/'//base_name//file_num//PID//'.pln'

	   NXP_L = NX-(NXP+1)*rank-1
	   NX_PAR = min(NXP,NXP_L)
	   if(rank.EQ.0)  write(6,'(2i6,2a)') k,kk,' ',file_name_1

       open(22,file=file_name_1,status='old',form='unformatted')
!        print *, file_name_1
       read(22) time, dt, tstep, xloc, nip, njp
!        if(rank.EQ.0) then
!        	write(*,'(2f8.4,i6,f8.4,2i6)') time, dt, tstep, xloc, nip, njp
!        endif
       if(njp.NE.nj)	then 
       	print *, '!!! WARNING !!! Check NJ', nj,'NJP', njp
       	STOP
       endif
       en=st+nip-1
       read(22) y(1,:),dx,u(st:en,:), 
     &			w(st:en,:),v(st:en,:), 
     &			p(st:en,:),th(st:en,:)
       close(22)
       st=en+1
	 ENDDO
	 x(1,:) = 0.d0
	 DO i=2,ni
	 	y(i,:) = y(1,:)
	 	x(i,:) = x(i-1,:)+dx
	 ENDDO
! 	 print *, 'dx ', dx
! 	 print *, 'End index:', en, 'nx:', ni
 
!       read(20,*) time, dt,u_bulk,v_bulk,u_bc,v_bc 
!       write(6,*)k,kk,dt,time
      
!       t(kk) = time


      do i=1,ni       
       do j=nj,1,-1
        U_data(kk-1,i-1,j-1)   = u(i,j) 
       end do
      end do	! end of loop over 'i' index
!_____________________________________________________________________       
       
      IF(VEL_FIELD) THEN
      !
      ! Open the file and write the tecplot datafile header information.
      !
      file_name_1 = 'yz_plane/plt/'//base_name//file_num
     &				//'.plt'//nulchar
      print *, file_name_1

      write(title,'(a6,f12.6)') 'time =',TIME
       ier = tecini(title,
     &        'y,z,u,v,w,p,rho_d'//nulchar,
     &         file_name_1//nulchar,
     &             '.'//nulchar,
     &             debug,visdouble)


       ier = teczne(title,
     &             imax,jmax,kmax,
     &             'BLOCK'//nulchar,nulchar)

      ! Write out the field data.
      itot = imax*jmax*kmax
      ier = tecdat(itot,x,disdouble)
      ier = tecdat(itot,y,disdouble)
      ier = tecdat(itot,u,disdouble)
      ier = tecdat(itot,v,disdouble)
      ier = tecdat(itot,w,disdouble)
      ier = tecdat(itot,p,disdouble)
      ier = tecdat(itot,th,disdouble)

      ier = tecend()
      ENDIF
        
	ENDDO		! END OF LOOPING OVER FILES (i.e. end of TIME loop)
!_____________________________________________________________________       
!_____________________________________________________________________       

      write(6,*)'Total steps',  kk

    
! 
! 122   format(I8,9f10.5)  
! 123   format(I8,8f12.5)
! 124   format(I8,15f12.5)
! 125   format(I8,12f12.5)
! 126   format(I8,8f12.5)      
! 
!       time_end = time

!       DO k=0,NKX
! ! 		do i=1,ni
! ! 			do j=1,nj							
!         		U_data(k,:,:)   = dsin(2.d0*t(k))
!         		print *, 'k,NKX,time=',k,NKX,dsin(2.d0*t(k))
! !         	enddo
! !         enddo
!       ENDDO 

! 		STOP
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       FFT 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       LX = time_end-time_st   
!       
!       CALL INIT_FFT_TIME 
! 
!       CALL FFT_X_TO_FOURIER_TIME(U_data,CU_data,0,NY+1,0,NZ+1) 
! 
! !      CALL FFT_X_TO_FOURIER_TIME(U_test,CU_test,0,NY+1,0,NZ+1)
! 
! 
!       int_dissp = 0.0 
!       int_prod  = 0.0
! 
!       DO k=0,NKX
!       int_dissp = 0.0
!       do i=ind_lt,ind_rt
!         dx=x(i+1,j)-x(i,j)
!         do j=1,ind_j_y(i)
!          dz=y(i,j+1)-y(i,j)
!          int_dissp = int_dissp + CDABS(CU_data(k,i,j))*dx*dz/area_j
!         enddo
!       enddo
! 
! 
!        WRITE(34,222) KX(K),int_dissp,
!      &              CDABS(CU_data(K,i_md,18)),
!      &              CDABS(CU_data(K,i_md,200)),
!      &              CDABS(CU_data(K,i_md,300)),
!      &              SUM(CDABS(CU_data(K,i_md,2:20)))/19.0
!       enddo
! 
! 
! !      DO I=0,NKX
! !      WRITE(34,222) KX(I),CDABS(CU_data(I,ni,nj)),CDABS(CU_test(I,1,1)),
! !     &              CDABS(CU_test(I,1,2)),CDABS(CU_test(I,2,1)),
! !     &              CDABS(CU_test(I,2,2))
! !      ENDDO
!  
! 222   format(6f12.6)
! 
!       CLOSE(34)
!         stop
!         
        
        end

