      module mg_simpleop

      implicit none

      use mg_grids

      contains 

      !----------------------------------------
      subroutine add_to(lev,x,y)
      ! x += y
      !
      ! todo: openmp the loop
      do k=1,nz 
         do j=1,ny
            do i=1,nx
               x(i,j,k) = x(i,j,k)+y(i,j,k)
            enddo
         enddo
      enddo

      !----------------------------------------
      subroutine set_to_zero(lev,x)
      ! x = 0.
      !
      ! todo: openmp the loop
      do k=1,nz 
         do j=1,ny
            do i=1,nx
               x(i,j,k) = 0.
            enddo
         enddo
      enddo

      end module
