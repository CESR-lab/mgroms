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
  subroutine correct_uvw(dx,dy,rmask,zr,zw,u,v,w)

    real(kind=rp), dimension(:,:)  , pointer, intent(in) :: dx,dy,rmask
    real(kind=rp), dimension(:,:,:), pointer, intent(in) :: zr,zw
    real(kind=rp), dimension(:,:,:), pointer, intent(inout) :: u,v,w

    integer(kind=ip):: k, j, i
    integer(kind=ip):: nx, ny, nz

    real(kind=rp), dimension(:,:)  ,   pointer :: umask,vmask
    real(kind=rp), dimension(:,:)  ,   pointer :: dxu,dyv
    real(kind=rp), dimension(:,:,:),   pointer :: dz,dzw
    real(kind=rp), dimension(:,:,:),   pointer :: p

    nx = grid(1)%nx
    ny = grid(1)%ny
    nz = grid(1)%nz

    if (myrank==0) write(*,*)'- correct u,v,w:'

    !! umask and vmask
    allocate(umask(0:ny+1,0:nx+1))
    allocate(vmask(0:ny+1,0:nx+1))
!!$    be care
    do i = 1,nx+1
       do j = 1,ny+1
          umask(j,i) = rmask(j,i-1)*rmask(j,i)
          vmask(j,i) = rmask(j-1,i)*rmask(j,i)
       enddo
    enddo
!!$      call fill_halo(lev,umask)
!!$      call fill_halo(lev,vmask)

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
          dxu(j,i) = 0.5_8*(dx(j,i)+dx(j,i-1))
       enddo
    enddo
    allocate(dyv(ny+1,0:nx+1))
    do i = 0,nx+1
       do j = 1,ny+1
          dyv(j,i) = 0.5_8*(dy(j,i)+dy(j-1,i))
       enddo
    enddo

    !! Correct
    p => grid(1)%p

    do i = 1,nx+1
       do j = 1,ny+1 
          do k = 1,nz
             u(k,j,i) = u(k,j,i) - 1._8/dxu(j,i)*(p(k,j,i)-p(k,j,i-1)) *umask(j,i)
             v(k,j,i) = v(k,j,i) - 1._8/dyv(j,i)*(p(k,j,i)-p(k,j-1,i)) *vmask(j,i)
          enddo
       enddo
    enddo

    do i = 1,nx
       do j = 1,ny 
          do k = 2,nz !interior and upper levels
             w(k,j,i) = w(k,j,i) - 1._8/dzw(k,j,i)*(p(k,j,i)-p(k-1,j,i))
          enddo
          k = nz+1 !surface
          w(k,j,i) = w(k,j,i) - 1._8/dzw(k,j,i)*(-p(k-1,j,i))
       enddo
    enddo

  end subroutine correct_uvw

  !-------------------------------------------------------------------------     
  subroutine check_correction(dx,dy,rmask,zr,zw,u,v,w)

    real(kind=rp), dimension(:,:)  , pointer, intent(in) :: dx,dy,rmask
    real(kind=rp), dimension(:,:,:), pointer, intent(in) :: zr,zw
    real(kind=rp), dimension(:,:,:), pointer, intent(inout) :: u,v,w

    if (myrank==0) write(*,*)'- check correction:'

    call compute_rhs(dx,dy,rmask,zr,zw,u,v,w)

    if (netcdf_output) then
       call write_netcdf(grid(1)%b,vname='b',netcdf_file_name='check_correction.nc',rank=myrank)
    endif

  end subroutine check_correction

end module mg_correct_uvw
