      module mg_smoother
      ! where smoother and direct operators are defined
      !

      implicit none

      use mg_grids
      use mg_optimized


      !----------------------------------------
      subroutine define_smoothers()
      
      ! the finest grid has to be defined in relation with the outside
      ! todo: think about it
      
      do lev=1,nlevs-1

         select case(grid(lev)%coarsening_method)

      case(1)
         ! regular 3D coarsening
         call restrict_matrix_xyz(lev)


      case(2)
         ! 2D coarsening

      case(3)
         ! simple vertical coarsening

         ! no fill halo because no coupling in the horizontal!
      case(4)
         ! double vertical coarsening

      case(5)
         ! triple vertical coarsening
         call restrict_matrix_zzz(lev) ! this gonna be tricky with the 22 points stencil!...

      end select
         
      enddo

      end subroutine

      !----------------------------------------
      subroutine smooth(lev,x,b,nite)
      ! driver because several smoothers are available

      select case(grid(lev)%smoothing_method)

      case (0) ! SOR 15 points stencil on the finest grid
         call smooth_sor_15(lev,x,b,nite)

      case (1) ! implicit vertical with the 15 points stencil
         call smooth_implicit_15(lev,x,b,nite)

      case (2) ! SOR 5 points stencil for 2D horizontal grids
         call smooth_sor_5(lev,x,b,nite)

      end select

      end subroutine

      !----------------------------------------
      subroutine residual(lev,x,b,r)

      select case(grid(lev)%smoothing_method)

      case (0,1) ! 15 points stencil
         call residual_15(lev,x,b,r)

      case (2) ! 5 points stencil (purely horizontal)
         call residual_5(lev,x,b,r)
         ! etc

      end select
 
      end subroutine


      end module
      
