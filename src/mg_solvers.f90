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
  subroutine solve(tol,maxite)
    real(kind=rp)   , intent(in) :: tol
    integer(kind=ip), intent(in) :: maxite
 
    ! local
    real(kind=rp)    :: rnorm,bnorm,res0,conv
    integer(kind=ip) :: nite

    bnorm = maxval(abs(grid(1)%b))
    call global_max(bnorm)

    call compute_residual(1,rnorm) ! residual returns both 'r' and its norm

    if (myrank == 0) write(*,*)' rnom:', rnorm,' bnorm:', bnorm

    res0 = rnorm/bnorm

    nite=0

    do while ((nite < maxite).and.(res0 > tol))

       call Fcycle()
       call compute_residual(1,rnorm)
       if (myrank == 0) write(*,*)' rnom:', rnorm,' bnorm:', bnorm
       rnorm = rnorm/bnorm
       conv=res0/rnorm ! error reduction after this iteration
       res0=rnorm
       nite=nite+1
       write(*,10) nite, rnorm, conv

    enddo

10  format("ite = ",I4,": res = ",G," / conv = ",G)

  end subroutine solve

  !---------------------------------------------------------------------
  subroutine Fcycle()

    integer(kind=ip):: lev

    do lev=1,nlevs-1
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

    integer(kind=ip),intent(in):: lev1

    integer(kind=ip):: lev
    real(kind=rp)   :: rnorm

    do lev=lev1,nlevs-1
       call relax(lev,ns_pre)
       call compute_residual(lev,rnorm)
       if (myrank == 0) write(*,*)' vcycle lev:', lev,' rnorm:', rnorm
       call fine2coarse(lev)
       grid(lev+1)%p(:,:,:) = 0._8
    enddo

    call relax(nlevs, ns_coarsest)

    do lev=nlevs-1,lev1,-1
       call coarse2fine(lev)
       call relax(lev,ns_post)
    enddo

  end subroutine Vcycle

end module mg_solvers