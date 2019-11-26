      Program main
!      INCLUDE 'header_fft'

      IMPLICIT NONE
      INTEGER   NX, NY, NZ, N_TH 

      PARAMETER (NX=32)
      PARAMETER (NY=1153-2)
      PARAMETER (NZ=257-2)
      PARAMETER (N_TH=1)

      integer ni,nj,nf,nk,nk_st,n_t, st,en,tstep
      integer nkp,nip,njp,NXP,NXV,NXP_L,NX_PAR,rank,NP
      parameter (nf=NX, ni=NZ+2,nj=NY+2,nk_st=1990,nk=2990,NP=4)
      parameter (n_t = nk-nk_st-1)
      PARAMETER (NXV=ceiling(real(NX+2)/real(NP))*NP)
      PARAMETER (NXP=NXV/NP-1)
      Real*8 x(ni,nj),y(ni,nj)
      Real*8 u(nf,ni,nj),v(nf,ni,nj),w(nf,ni,nj),th(nf,ni,nj),
     &       th_wh(nf,ni,nj)

      real*8   pi,time,RI,U_0,omega,rho_0,nu,dt,dx,xloc
      integer  i,j,imax,jmax,kmax,k,kk,i_md,ind_lt,ind_rt
      integer  zindex,zindex2,xindex,xindex2,index_Mfac
      parameter (zindex=201,zindex2=240,xindex=341,xindex2=107)
      
      integer  debug,ier,itot
      integer  tecini,tecdat,teczne,tecnod,tecfil,tecend
      integer  visdouble,disdouble
      character*1 nulchar

      CHARACTER*4 PID
      CHARACTER*4 file_num
      CHARACTER*34 file_name_1,file_time,file_name_2
      CHARACTER*6 base_name
      character(len=18) :: title

      LOGICAL VEL_FIELD
!       PARAMETER (NKXV=ceiling(real(NKX+1)/real(NP))*NP) 
!       PARAMETER (NKXP=NKXV/NP-1)
!       PARAMETER (NX2V=ceiling(real(NXV/2)/real(NP))*NP)
!       PARAMETER (NX2P=NX2V/NP-1)    

      base_name = 'data_'
      print *, 'folder name=', base_name
      
      PI  = ATAN(1.0)*4.0

      

!      nulchar = char(0)
!      debug   = 0
!      visdouble = 0
!      disdouble = 1
!      imax = ni
!      jmax = nj
!      kmax = 1        
!      kk = 0

       
      kk = 0
       
      DO k=nk_st,nk
          kk= kk+1

      file_num = CHAR(MOD(k,10000)/1000+48) 
     &          //CHAR(MOD(k,1000)/100+48)     
     &          //CHAR(MOD(k,100)/10+48)       
     &          //CHAR(MOD(k,10)+48)
       
      st=1
     
      DO rank=0,NP-1
        I=rank+1
       PID = '_'//CHAR(MOD(I,1000)/100+48)   
     &          //CHAR(MOD(I,100)/10+48)    
     &          //CHAR(MOD(I,10)+48)
       
      file_name_1=trim('plane_3D/'//trim(base_name) 
     &                  //file_num//PID//'.pln')
      file_name_2='plane_3D_test/'//trim(base_name)//file_num//'.pln'

	   NXP_L = NX-(NXP+1)*rank-1
	   NX_PAR = min(NXP,NXP_L)
	   if(rank.EQ.0)  write(6,'(2i6,2a)') k,kk,' ',file_name_1

       open(22,file=file_name_1,status='old',form='unformatted')
       print *, file_name_1
       read(22) time, dt, tstep, nip,nkp, njp
!        if(rank.EQ.0) then
!        	write(*,'(2f8.4,i6,f8.4,2i6)') time, dt, tstep, xloc, nip, njp
!        endif
       if((njp.NE.nj).or.((nkp.NE.ni)))then 
       	print *, '!!! WARNING !!! Check NJ', ni,'NJP', nkp
        print *, '!!! WARNING !!! Check NJ', nj,'NJP', njp
       	STOP
       endif

       en=st+nip-1
       read(22) dx,u(st:en,:,:), 
     &          w(st:en,:,:),v(st:en,:,:), 
     &          th(st:en,:,:)
       close(22)
       st=en+1
	 ENDDO


        print *, file_name_2
        open(23,file=file_name_2,status='unknown',form='unformatted')
        write(23) dx,u(:,:,:),w(:,:,:),v(:,:,:),
     &          th(:,:,:)
       close(23)
!_____________________________________________________________________       
       
        
	ENDDO		! END OF LOOPING OVER FILES (i.e. end of TIME loop)
!_____________________________________________________________________       
!_____________________________________________________________________       

      write(6,*)'Total steps',  kk

    
        
        end

