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
  subroutine correct_uvw(rmask,u,v,w)

    real(kind=rp), dimension(:,:)  , pointer, intent(inout) :: rmask
    real(kind=rp), dimension(:,:,:), pointer, intent(inout) :: u,v,w

    integer(kind=ip):: k, j, i
    integer(kind=ip):: nx, ny, nz

    real(kind=rp), dimension(:,:)  , pointer :: umask,vmask
    real(kind=rp), dimension(:,:)  , pointer :: dx,dy
    real(kind=rp), dimension(:,:,:), pointer :: zr,zw
    real(kind=rp), dimension(:,:,:), pointer :: p
    real(kind=rp) :: dxu,dyv
    real(kind=rp) :: dzw

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

    !! umask and vmask
    allocate(umask(0:ny+1,0:nx+1))
    allocate(vmask(0:ny+1,0:nx+1))
    if (bmask) then
!!$    be care
       umask(:,:)=0._8
       vmask(:,:)=0._8
       do i = 1,nx+1
          do j = 0,ny+1
             umask(j,i) = rmask(j,i-1)*rmask(j,i)
          enddo
       enddo
       do i = 0,nx+1
          do j = 1,ny+1
             vmask(j,i) = rmask(j-1,i)*rmask(j,i)
          enddo
       enddo
    else
       umask(:,:)=1._8
       vmask(:,:)=1._8
    endif

    !! Correct
    p => grid(1)%p

    do i = 1,nx+1
       do j = 0,ny+1 
          do k = 1,nz

             dxu = hlf * (dx(j,i)+dx(j,i-1))

             u(i,j,k) = u(i,j,k) - one / dxu * (p(k,j,i)-p(k,j,i-1)) * umask(j,i)

          enddo
       enddo
    enddo

    do i = 0,nx+1
       do j = 1,ny+1 
          do k = 1,nz

             dyv = hlf * (dy(j,i)+dy(j-1,i))

             v(i,j,k) = v(i,j,k) - one / dyv * (p(k,j,i)-p(k,j-1,i  )) * vmask(j,i)

          enddo
       enddo
    enddo

    do i = 0,nx+1
       do j = 0,ny+1

          do k = 2,nz !interior and upper levels
             dzw = zr(k,j,i)-zr(k-1,j,i)
             w(i,j,k-1) = w(i,j,k-1) - one / dzw * (p(k,j,i)-p(k-1,j,i))
          enddo

          k = nz+1 !surface
          dzw = zw(nz+1,j,i)-zr(nz,j,i)
          w(i,j,k-1) = w(i,j,k-1) - one / dzw * (-p(k-1,j,i))

       enddo
    enddo

    deallocate(umask)
    deallocate(vmask)

  end subroutine correct_uvw

end module mg_correct_uvw
