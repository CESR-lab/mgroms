module nhydro

  use mg_grids
  use mg_define_matrix
  use mg_restrict
  use mg_relax

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

    call define_grids(npxg, npyg, nz,ny,nx)

    call define_neighbours(neighb)

    call define_matrices()

  end subroutine nhydro_init

  !--------------------------------------------------------------
  subroutine nhydro_solve(u,v,w)
    real(kind=8), dimension(:,:,:), intent(inout) :: u,v,w

    real(kind=8) :: res

    do lev=1,nlevs-1
       !- discuss with Guillaume !!!
       call compute_residual(lev,res)
       call restrict(lev)
    enddo

    !- discuss with Guillaume !!!
    do lev=1,nlevs

       res0=0.

       do it=1, nit
          call relax(lev,nsweeps)
          call compute_residual(lev,res)

          call prolong(lev) ! interpolation

          conv = log(res0/res)/log(10.)
          res0=res

       enddo

    enddo

  end subroutine nhydro_solve

  !--------------------------------------------------------------
  subroutine nhydro_clean()
    
  end subroutine nhydro_clean

end module nhydro
