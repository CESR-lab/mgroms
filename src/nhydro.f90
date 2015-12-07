module nhydro

  use mg_mpi
  use mg_grids
  use mg_define_matrix
  use mg_intergrids
  use mg_relax
  use mg_solvers

  implicit none

contains

  !--------------------------------------------------------------
  subroutine nhydro_init(npxg, npyg, neighb, dx, dy, zr, zw)

    integer(kind=4), dimension(4) , intent(in) :: neighb ! S, E, N, W
    integer(kind=4)               , intent(in) :: npxg, npyg
    real(kind=8), dimension(:,:)  , intent(in) :: dx, dy
    real(kind=8), dimension(:,:,:), intent(in) :: zr, zw

    integer(kind=4) :: nz, ny, nx
    integer(kind=4) :: ierr, lev

    nz = size(zr,dim=1)
    ny = size(zr,dim=2)
    nx = size(zr,dim=3)

    call mg_mpi_init()

    call define_grids(npxg, npyg, nx, ny, nz)

    call define_neighbours(neighb)

    call MPI_Barrier( MPI_COMM_WORLD ,ierr)
    if (myrank.eq.0)then
       do lev=1,nlevs
          if (grid(lev)%gather.eq.0)then
             write(*,100)"lev=",lev,": ", &
                  grid(lev)%nx,' x',grid(lev)%ny,' x',grid(lev)%nz, &
                  " on ",grid(lev)%npx,' x',grid(lev)%npy," procs"
          else
             write(*,100)"lev=",lev,": ", &
                  grid(lev)%nx,' x',grid(lev)%ny,' x',grid(lev)%nz, &
                  " on ",grid(lev)%npx,' x',grid(lev)%npy," procs / gather"
          endif
       enddo
    endif
100 format (A4,I2,A,I3,A,I3,A,I3,A,I3,A,I3,A)

    call define_matrices()

  end subroutine nhydro_init

  !--------------------------------------------------------------
  subroutine nhydro_solve(u,v,w)
    real(kind=8), dimension(:,:,:), allocatable, intent(inout) :: u,v,w

    real(kind=8) :: tol = 1.e-6

    integer(kind=4) :: nx, ny, nz

    integer(kind=4) :: maxite =2

    nz = size(u,dim=1)
    ny = size(u,dim=2)
    nx = size(v,dim=3)

    !- we need a MPI update 
    !grid(1)%b = rhs(u,v,w) ! div of u,v,w

    grid(1)%b(1:nz,1:ny,1:nx) = 0._8
    call random_number(grid(1)%p)
    
!!        &
!!         u(:,:,2:nx+1) - u(:,:,1:nx) + &
!!         v(:,2:ny+1,:) - v(:,1:ny,:) + &
!!         w(2:nz+1,:,:) - w(1:nz,:,:)

    call solve(tol,maxite)

  end subroutine nhydro_solve

  !--------------------------------------------------------------
  subroutine nhydro_clean()
    
  end subroutine nhydro_clean

end module nhydro
