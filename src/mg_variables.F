      module mg_variables

      implicit none

      type grid_type
      real*8,dimension(:,:,:)  ,allocatable :: x,b,r
      real*8,dimension(:,:,:,:),allocatable :: A
      integer:: nx,ny,nz
      end type grid_type


      type(grid_type),dimension(:),allocatable :: grid
      integer:: nlevs ! index of the coarsest level (1 is the finest)


      contains


      subroutine allocate_variables

      integer:: lev,nx,ny,nz

      do lev=1,nlevs
         nx = grid(lev)%nx
         ny = grid(lev)%ny
         nz = grid(lev)%nz
         allocate( grid(lev)%x(nx,ny,nz) )
         allocate( grid(lev)%b(nx,ny,nz) )
         allocate( grid(lev)%r(nx,ny,nz) )
         allocate( grid(lev)%A(8,nx,ny,nz) )
      end subroutine

      end module
