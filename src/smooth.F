      subroutine smooth(lev,grid,nite)
      
      implicit none

      integer:: lev,nite
      ! mg_grid is a structure containing the variables of a given level
      ! lev
      mg_grid:: grid

      ! local
      integer:: k

      do k = 1,nite
         ! we assume halo width >=2 and so smooth twice before updating
         ! the halo
         call smooth_twice(grid%nx,grid%ny,grid%nz,grid%A,grid%x,grid%b)
         call fill_halo(lev,grid,grid%x)
      enddo

      end
      
