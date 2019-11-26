      Program main
      integer n,i, IOUT
      parameter (n=64)
      real*8,allocatable,dimension(:) :: a
      real*8 b(n),c(n),&
       d(n)
 
!$OMP PARALLEL PRIVATE
      allocate (a(n))

      DO i=1,n
      a(i)=i
      b(i)=2.0
      write(IOUT,*)a(i)
      write(*,*)a(i)
      enddo
        
      write(IOUT,*) real(34)/real(4),ceiling (real(34)/real(4))    
!$OMP END PARALLEL                  



      stop
      end
