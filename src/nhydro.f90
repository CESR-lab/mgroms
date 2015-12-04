module nhydro

  use mg_mpi
  use mg_grids
  use mg_define_matrix
  use mg_intergrids
  use mg_relax
  !TODO use mg_solvers

contains

  !--------------------------------------------------------------
  subroutine nhydro_init(npxg, npyg, neighb, dx, dy, zr, zw)

    integer(kind=4), dimension(4) , intent(in) :: neighb ! S, E, N, W
    integer(kind=4)               , intent(in) :: npxg, npyg
    real(kind=8), dimension(:,:)  , intent(in) :: dx, dy
    real(kind=8), dimension(:,:,:), intent(in) :: zr, zw

    integer(kind=4) :: nz, ny, nx

    nz = size(zr,dim=1)
    ny = size(zr,dim=2)
    nx = size(zr,dim=3)

    call mg_mpi_init()

    call define_grids(npxg, npyg, nz,ny,nx)

    call define_neighbours(neighb)

    call define_matrices()

  end subroutine nhydro_init

  !--------------------------------------------------------------
  subroutine nhydro_solve(u,v,w)
    real(kind=8), dimension(:,:,:), intent(inout) :: u,v,w

    !TODO call mg_solve()

  end subroutine nhydro_solve

  !--------------------------------------------------------------
  subroutine nhydro_clean()
    
  end subroutine nhydro_clean

end module nhydro
