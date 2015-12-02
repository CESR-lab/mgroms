module mg_intergrids
  ! define the interface to grid restriction/interpolation
  use mg_mpi 

  use mg_grids
  use mg_restrict
  use mg_gather

  implicit none

  contains

      !----------------------------------------
      subroutine finetocoarse(l1,l2,x,y)
      
      integer(kind=is):: l1,l2 ! l1 is the fine grid index, l2=l1+1

      real(kind=rl),dimension(:,:,:),intent(in) :: x
      real(kind=rl),dimension(:,:,:),intent(out) :: y

      select case (grid(l1)%gather)
         case(0)
            call restrict_generic(l1,l2,x,y) ! regular coarsening
         case(1)
            call restrict_generic(l1,l2,x,grid(l1)%dummy) ! coarsen
            call gather(l2,grid(l1)%dummy,y) ! then gather
      end select

      end subroutine

      !----------------------------------------
      subroutine restrict_generic(l1,l2,x,y)

      ! watch out: the dimensions of x or y do not correspond to a grid level
      ! if a gather is introduced

      integer(kind=is):: l1,l2 ! l1 is the fine grid index, l2=l1+1

      real(kind=rl),dimension(:,:,:),intent(in) :: x
      real(kind=rl),dimension(:,:,:),intent(out) :: y

      select case(grid(l1)%coarsening_method)
      case(1)
         ! regular 3D coarsening
         call restrict_xyz(l1,l2,x,y)
         call fill_halo(l2,y)
      case(2)
         ! 2D coarsening
         call restrict_xy(l1,l2,x,y)
         call fill_halo(l2,y)
!!$      case(3)
!!$         ! simple vertical coarsening
!!$         call restrict_z(l1,l2,x,y)
!!$         ! no fill halo because no coupling in the horizontal!
!!$      case(4)
!!$         ! double vertical coarsening
!!$         call restrict_zz(l1,l2,x,y)
      case(5)
         ! triple vertical coarsening
         call restrict_zzz(l1,l2,x,y)
      end select
      
      end subroutine

  end module mg_intergrids
  
  
