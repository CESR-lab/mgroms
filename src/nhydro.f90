module nhydro

  use mg_mpi
  use mg_grids
  use mg_namelist
  use mg_define_matrix
  use mg_intergrids
  use mg_relax
  use mg_solvers
  use mg_mpi_exchange

  implicit none

contains

  !--------------------------------------------------------------
  subroutine nhydro_init(nx, ny, nz, npxg, npyg, neighb, dx, dy, zr, zw, umask, vmask)

    integer(kind=ip), dimension(4)              , intent(in) :: neighb ! S, E, N, W
    integer(kind=ip)                            , intent(in) :: nx, ny, nz
    integer(kind=ip)                            , intent(in) :: npxg, npyg
    real(kind=rp), dimension(:,:)  , allocatable, intent(in) :: dx, dy
    real(kind=rp), dimension(:,:)  , allocatable, intent(in) :: umask, vmask
    real(kind=rp), dimension(:,:,:), allocatable, intent(in) :: zr, zw

    integer(kind=ip) :: ierr, lev

    call mg_mpi_init()

    call define_grids(npxg, npyg, nx, ny, nz)

    call define_neighbours(neighb)

    call print_grids()

    call define_matrices(dx, dy, zr, zw, umask, vmask)

  end subroutine nhydro_init

  !--------------------------------------------------------------
  subroutine nhydro_solve(u,v,w)

    real(kind=rp), dimension(:,:,:), allocatable, intent(inout) :: u,v,w

    integer(kind=ip) :: nx, ny, nz
    real(kind=rp)    :: tol = 1.e-12
    integer(kind=ip) :: maxite = 30

    nz = size(u,dim=1)
    ny = size(u,dim=2)
    nx = size(v,dim=3)

    !- we need a MPI update 
    ! grid(1)%b = rhs(u,v,w) ! div of u,v,w

    grid(1)%p(:,:,:) = 0.

!ND
!    call relax(1,1)
    call solve(tol,maxite)

    !ND
    call write_netcdf(grid(1)%p,vname='p',netcdf_file_name='p.nc',rank=myrank)
    call write_netcdf(grid(1)%r,vname='r',netcdf_file_name='r.nc',rank=myrank)

  end subroutine nhydro_solve

  !--------------------------------------------------------------
  subroutine nhydro_clean()

    call grids_dealloc()
    
  end subroutine nhydro_clean

!!$  !-------------------------------------------------------------------------     
!!$  subroutine Kinetic_energy(u,v,w,uf,vf,wf);
!!$
!!$    real(kind=rp)   :: dimension(:,:,:),intent(in):: u,v,w
!!$    real(kind=rp)   :: dimension(:,:,:),intent(in):: uf,vf,wf
!!$
!!$    real(kind=rp)    :: Ek = 0.0_8
!!$    integer(kind=ip) :: i,j,k
!!$
!!$ 
!!$    do i = 1,nx
!!$       do j = 1,ny
!!$          do k = 1,nz
!!$              Ek = Ek + 
!!$   &            0.5*( u(k,j,i  )*uf(k,j,i  )*dxu(j,i  )*dyu(j,i  )*dzu(k,j,i  ) + 
!!$   &                  u(k,j,i+1)*uf(k,j,i+1)*dxu(j,i+1)*dyu(j,i+1)*dzu(k,j,i+1) ) + 
!!$   &            0.5*( v(k,j  ,i)*vf(k,j  ,i)*dxv(j  ,i)*dyv(j  ,i)*dzv(k,j  ,i) + 
!!$   &                  v(k,j+1,i)*vf(k,j+1,i)*dxv(j+1,i)*dyv(j+1,i)*dzv(k,j+1,i) ) + 
!!$   &            0.5*( w(k  ,j,i)*wf(k  ,j,i)*dx(j,i)*dy(j,i)*dzw(k  ,j,i) + 
!!$   &                  w(k+1,j,i)*wf(k+1,j,i)*dx(j,i)*dy(j,i)*dzw(k+1,j,i) ) + 
!!$          enddo
!!$       enddo
!!$    enddo
!!$
!!$  end subroutine Kinetic Energy
!!$!----------------------------------------
!!$  subroutine Momentum2flux(u,v,w,uf,vf,wf,zx)
!!$    !!
!!$    !!  uf =  um - wm*zx
!!$    !!  vf =  vm - wm*zy
!!$    !!  wf = -zx*um - zx*vm + (1+zx^2+ zy^2)*wm
!!$    !!
!!$    real(kind=rp)   :: dimension(:,:,:),intent(in) :: u,v,w
!!$    real(kind=rp)   :: dimension(:,:,:),intent(out):: uf,vf,wf
!!$    real(kind=rp)   :: dimension(:,:,:),intent(in) :: zx
!!$
!!$    !! Slopes are defined at rho-points
!!$    do i = 1,nx+1
!!$       do j = 1,ny
!!$          do k = 1,nz
!!$            uf(k,j,i) = um(k,j,i) - 0.25*(zx(k,j,i  )*wm(k,j,i  ) + zx(k,j,i  )*wm(k+1,j,i  )
!!$                                          zx(k,j,i-1)*wm(k,j,i-1) + zx(k,j,i-1)*wm(k+1,j,i-1) )
!!$            vf(k,j,i) = vm(k,j,i) - 0.25*(zy(k,j  ,i)*wm(k,j  ,i) + zx(k,j  ,i)*wm(k+1,j  ,i)
!!$                                          zy(k,j-1,i)*wm(k,j-1,i) + zx(k,j  ,i)*wm(k+1,j  ,i) )
!!$            wf(k,j,i) = vm(k,j,i) - 0.25*(zy(k,j  ,i)*wm(k,j  ,i) + zx(k,j  ,i)*wm(k+1,j  ,i)
!!$                                          zy(k,j-1,i)*wm(k,j-1,i) + zx(k,j  ,i)*wm(k+1,j  ,i) )
!!$       enddo
!!$       zw(nz+1,j,i) = 0.0
!!$    enddo
!!$
!!$  end subroutine Momentum2flux

end module nhydro
