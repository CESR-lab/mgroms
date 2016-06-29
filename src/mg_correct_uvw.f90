module mg_correct_uvw

  use mg_mpi
  use mg_tictoc
  use mg_namelist
  use mg_compute_rhs
  use mg_grids
  use mg_mpi_exchange
  use mg_netcdf_out

  implicit none

contains
  !-------------------------------------------------------------------------     
  subroutine correct_uvw(u,v,w)

    real(kind=rp), dimension(:,:,:), pointer, intent(inout) :: u,v,w

    integer(kind=ip):: k, j, i
    integer(kind=ip):: nx, ny, nz

    real(kind=rp), dimension(:,:)  , pointer :: dx,dy
    real(kind=rp), dimension(:,:,:), pointer :: zr,zw
    real(kind=rp), dimension(:,:)  , pointer :: dxu,dyv
    real(kind=rp), dimension(:,:,:), pointer :: dz,dzw
    real(kind=rp), dimension(:,:,:), pointer :: p

    !NG comment: constants in a mg_cst.f90 file ?
    real(kind=rp), parameter :: two  = 2._rp
    real(kind=rp), parameter :: one  = 1._rp
    real(kind=rp), parameter :: hlf  = 0.5_rp
    real(kind=rp), parameter :: qrt  = 0.25_rp
    real(kind=rp), parameter :: zero  = 0._rp

    nx = grid(1)%nx
    ny = grid(1)%ny
    nz = grid(1)%nz

    dx => grid(1)%dx
    dy => grid(1)%dy
    zr => grid(1)%zr
    zw => grid(1)%zw

    if (myrank==0) write(*,*)'- correct u,v,w:'

    !! Cell heights
    allocate(dzw(nz+1,0:ny+1,0:nx+1))
    do i = 0,nx+1
       do j = 0,ny+1
          dzw(1,j,i) = zr(1,j,i)-zw(1,j,i) !!
          do k = 2,nz
             dzw(k,j,i) = zr(k,j,i)-zr(k-1,j,i) !!  cell height at w-points
          enddo
          dzw(nz+1,j,i) = zw(nz+1,j,i)-zr(nz,j,i) !!
       enddo
    enddo

    !! Cell widths
    allocate(dxu(0:ny+1,nx+1))
    do i = 1,nx+1
       do j = 0,ny+1
          dxu(j,i) = hlf * (dx(j,i)+dx(j,i-1))
       enddo
    enddo
    allocate(dyv(ny+1,0:nx+1))
    do i = 0,nx+1
       do j = 1,ny+1
          dyv(j,i) = hlf * (dy(j,i)+dy(j-1,i))
       enddo
    enddo

    !! Correct
    p => grid(1)%p

    do i = 1,nx+1
       do j = 1,ny+1 
          do k = 1,nz
             u(k,j,i) = u(k,j,i) - one / dxu(j,i)*(p(k,j,i)-p(k,j,i-1))
             v(k,j,i) = v(k,j,i) - one / dyv(j,i)*(p(k,j,i)-p(k,j-1,i))
          enddo
       enddo
    enddo

    do i = 1,nx
       do j = 1,ny 
          do k = 2,nz !interior and upper levels
             w(k,j,i) = w(k,j,i) - one / dzw(k,j,i)*(p(k,j,i)-p(k-1,j,i))
          enddo
          k = nz+1 !surface
          w(k,j,i) = w(k,j,i) - one / dzw(k,j,i)*(-p(k-1,j,i))
       enddo
    enddo

  end subroutine correct_uvw

  !-------------------------------------------------------------------------     
  subroutine check_correction(nx,ny,nz,ua,va,wa)

    integer(kind=ip), intent(in) :: nx, ny, nz
    real(kind=rp), dimension(nz  ,0:ny+1,0:nx+1), target, intent(in) :: ua
    real(kind=rp), dimension(nz  ,0:ny+1,0:nx+1), target, intent(in) :: va
    real(kind=rp), dimension(nz+1,0:ny+1,0:nx+1), target, intent(in) :: wa

    real(kind=rp), dimension(:,:,:), pointer :: u, v, w

    if (myrank==0) write(*,*)'- check correction:'

    u => ua
    v => va
    w => wa

    write(*,*)'check_correction -> Lbound(v):',lbound(v)

    call compute_rhs(u,v,w)

  end subroutine check_correction

end module mg_correct_uvw
