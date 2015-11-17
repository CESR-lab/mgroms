      module mg_interpolation
      ! where restriction and interpolations are defined
      !
      ! there are two drivers: finetocoarse and coarsetofine
      ! switching to optimized routines
      ! we use linear interpolations (not nearest value)

      implicit none

      use mg_grids
      use mg_mpi
      use mg_optimized

      !----------------------------------------
      subroutine finetocoarse(l1,l2,x,y)
      
      integer:: l1,l2 ! l1 is the fine grid index, l2=l1+1

      real,dimension(grid(l1)%nx,grid(l1)%ny,grid(l1)%nz) :: x
      real,dimension(grid(l2)%nx,grid(l2)%ny,grid(l2)%nz) :: y

      select case (grid(l1)%gather)
         case(0)
            call restrict_generic(x,y) ! regular coarsening
         case(1)
            call gather(l1,x,grid(l1)%dummy) ! gather
            call restrict_generic(grid(l1)%dummy,y) ! then coarsen
         case(2)
            call restrict_generic(x,grid(l1)%dummy) ! coarsen
            call gather(grid(l1)%dummy,y) ! then gather
      end select

      end subroutine

      !----------------------------------------
      subroutine restrict_generic(x,y)

      ! watch out: the dimensions of x or y do not correspond to a grid level
      ! if a gather is introduced

      select case(grid(l1)%coarsening_method)
      case(1)
         ! regular 3D coarsening
         call restrict_xyz(l1,l2,x,y)
         call fill_halo(l2,y)
      case(2)
         ! 2D coarsening
         call restrict_xy(l1,l2,x,y)
         call fill_halo(l2,y)
      case(3)
         ! simple vertical coarsening
         call restrict_z(l1,l2,x,y)
         ! no fill halo because no coupling in the horizontal!
      case(4)
         ! double vertical coarsening
         call restrict_zz(l1,l2,x,y)
      case(5)
         ! triple vertical coarsening
         call restrict_zzz(l1,l2,x,y)
      end select
      
      end subroutine



      end module
      
