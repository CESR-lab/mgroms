      subroutine smooth(lev,nite)
      ! apply nite*2 smoothings at level lev
      ! smoothing operator is stored in grid(lev)%A

      use mg_variables 
!     provides the array of structure 'grid'
!     each entry of grid contains all the variables of a given level

      implicit none

      integer:: lev,nite

      ! local
      integer:: k

      do k = 1,nite
         ! we assume halo width >=2 and so smooth twice before updating
         ! the halo
         call smooth_twice(grid(lev)%nx,grid(lev)%ny,grid(lev)%nz,
     $                     grid(lev)%A,grid(lev)%x,grid(lev)%b)
         call fill_halo(lev,grid(lev)%x)
      enddo

      end subroutine
      

      subroutine smooth_twice(nx,ny,nz,A,x,b)

      implicit none

      integer:: nx,ny,nz
      real*8,dimensions(8,nx,ny,nz)::A
      real*8,dimensions(nx,ny,nz)::x,b

      ! local
      integer:: i,j,k
      real*8  ::coef


      ! do we have ghostpoints in the vertical? => likely no ...
      !
      do k=1,nz
         do j=2,ny-1
            do i=2,nx-1
               coef = 1/abs(A(8,i,j,k))

               x(i,j,k)=x(i,j,k) + coef*( -b(i,j,k)
     $              +A(1,i,j,k)*x(i-1,j-1,k)
     $              +A(2,i,j,k)*x(i-1,j  ,k-1)
     $              +A(3,i,j,k)*x(i  ,j-1,k)
! ...
     $              +A(8,i,j,k)*x(i  ,j  ,k) )
              
            enddo
         enddo
      enddo

      end subroutine
    
      
