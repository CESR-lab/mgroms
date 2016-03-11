module mg_solvers

  use mg_mpi
  use mg_tictoc
  use mg_namelist
  use mg_grids
  use mg_relax
  use mg_intergrids

  implicit none

contains

  !---------------------------------------------------------------------
  subroutine testgalerkin(lev)
    real(kind=8) :: norm_c,norm_f,dummy
    integer(kind=4) :: lev,nx,ny,nz,nh

    nx = grid(lev)%nx
    ny = grid(lev)%ny
    nz = grid(lev)%nz
    nh = grid(lev)%nh

    call random_number(grid(lev)%p)!
    grid(lev)%p(:,1,:)= 1._8
!    grid(lev)%p(2,:,:)=0._8
    call fill_halo(lev,grid(lev)%p)
    grid(lev)%b = 0._8
    call compute_residual(lev,dummy)    
    call norm(lev,grid(lev)%p,grid(lev)%r,nx,ny,nz,norm_c)

    grid(lev-1)%p = 0._8 
    call coarse2fine(lev-1) ! interpolate p to r and add r to p
    grid(lev-1)%b = 0._8
    call compute_residual(lev-1,dummy)

    nx = grid(lev-1)%nx
    ny = grid(lev-1)%ny
    nz = grid(lev-1)%nz
    nh = grid(lev-1)%nh
    call norm(lev-1,grid(lev-1)%p,grid(lev-1)%r,nx,ny,nz,norm_f)

    if (myrank==0) then
       write(*,*)"======== lev ",lev,"==========="
       write(*,*)"norm coarse = ",norm_c
       write(*,*)"norm fine   = ",norm_f/16._8
    endif

end subroutine testgalerkin


  !---------------------------------------------------------------------
subroutine norm(lev,x,y,nx,ny,nz,res)
  use mg_mpi_exchange
  integer(kind=4) :: lev,i,j,k
  integer(kind=4) :: nx,ny,nz
  real(kind=8) :: r,res
  real(kind=8),dimension(:,:,:)  , pointer :: x,y

  r=0._8
  do i=1,nx
     do j=1,ny
        do k=1,nz
           r=r+x(k,j,i)*y(k,j,i)
        enddo
     enddo
  enddo
  call global_sum(lev,r,res)
end subroutine norm

  !---------------------------------------------------------------------
  subroutine solve(tol,maxite)

    real(kind=rp)   , intent(in) :: tol
    integer(kind=ip), intent(in) :: maxite
 
    ! local
    real(kind=rp)    :: rnorm,bnorm,res0,conv,rnorm0
    integer(kind=ip) :: nite

    integer(kind=ip) :: nx,ny,nz,nh
    real(kind=rp), dimension(:,:,:), allocatable :: p0,b0
    real(kind=rp), dimension(:,:,:), pointer :: p,b,r

    real(kind = lg) :: tstart,tend,perf
    real(kind=rp) :: rnxg,rnyg,rnzg
    real(kind=rp) :: rnpxg,rnpyg

    p  => grid(1)%p
    b  => grid(1)%b
    r  => grid(1)%r

    nx = grid(1)%nx
    ny = grid(1)%ny
    nz = grid(1)%nz
    nh = grid(1)%nh

    allocate(p0(nz,1-nh:ny+nh,1-nh:nx+nh))
    allocate(b0(nz,1-nh:ny+nh,1-nh:nx+nh))

    p0 = p
    b0 = b

    call tic(1,'solve')
    call cpu_time(tstart)
    
!    bnorm = maxval(abs(grid(1)%b))
!    call global_max(bnorm)

    res0 = sum(grid(1)%b(1:nz,1:ny,1:nx)**2)
    call global_sum(1,res0,bnorm)

    call compute_residual(1,rnorm) ! residual returns both 'r' and its norm
    
    if (myrank == 0) write(*,*)'rnom:', rnorm,' bnorm:', bnorm

    res0 = sqrt(rnorm/bnorm)
    rnorm0 = res0

    nite=0

    do while ((nite < maxite).and.(res0 > tol))

       call Vcycle(1)
       !call relax(1,1)

       call compute_residual(1,rnorm)
       rnorm = sqrt(rnorm/bnorm)
       conv = res0/rnorm ! error reduction after this iteration
       res0 = rnorm

       nite = nite+1
       if (myrank == 0) write(*,10) nite, rnorm, conv
       if (myrank == 0) write(100,*) rnorm, conv
    enddo

    call cpu_time(tend)
    call toc(1,'solve')

    if (myrank == 0) then
       rnpxg=real(grid(1)%npx,kind=rp)
       rnpyg=real(grid(1)%npy,kind=rp)
       rnxg=real(grid(1)%nx,kind=rp)*rnpxg
       rnyg=real(grid(1)%ny,kind=rp)*rnpyg
       rnzg=real(grid(1)%nz,kind=rp)
       ! the rescaled time should be expressed in terms of error reduction,
       ! therefore the ratio rnorm/rnorm0 [the rnorm0 was missing prior Dec 11th]
       perf = (tend-tstart)*(rnpxg*rnpyg)/(-log(rnorm/rnorm0)/log(10._8))/(rnxg*rnyg*rnzg)
       write(*,*)'--- summary ---'
       write(*,'(A,F8.3,A)')"time spent to solve :",tend-tstart," s"
       write(*,'(A,E10.3)')"rescaled performance:",perf
       write(*,*)'---------------'
    end if

10  format("ite = ",I2,": res = ",E10.3," / conv = ",F7.1)

  end subroutine solve

  !---------------------------------------------------------------------
  subroutine Fcycle()

    integer(kind=ip):: lev

    do lev=1,nlevs-1
       grid(lev)%r=grid(lev)%b
       call fine2coarse(lev)
    enddo

    grid(nlevs)%p(:,:,:) = 0._8

    call relax(nlevs, ns_coarsest)

    do lev=nlevs-1,1,-1
       grid(lev)%p(:,:,:) = 0._8
       call coarse2fine(lev) !- fill_halo ?
       call Vcycle(lev)
    enddo

  end subroutine Fcycle

  !----------------------------------------
  subroutine Vcycle(lev1)

    integer(kind=ip), intent(in) :: lev1
    integer(kind=ip)             :: lev, nlevs0
    real(kind=rp)                :: rnorm

    nlevs0=nlevs

    do lev=lev1,nlevs0-1
       call relax(lev,ns_pre)
       call compute_residual(lev,rnorm)
       call fine2coarse(lev)
       grid(lev+1)%p(:,:,:) = 0._8
    enddo

    call relax(nlevs0,ns_coarsest)

    do lev=nlevs0-1,lev1,-1
       call coarse2fine(lev)
       call relax(lev,ns_post)
    enddo

  end subroutine Vcycle

end module mg_solvers
