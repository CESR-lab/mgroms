module mg_solvers

  use mg_tictoc
  use mg_grids
  use mg_smoother
  use mg_interpolation     
  use mg_simpleop

  implicit none

contains

  !----------------------------------------
  subroutine mg_solve(x,b,tol,maxite,res,nite)
    !     input: x=first guess / b=RHS
    !     output: x(overwritten) = solution


    real(kind=8),dimension(:,:,:), intent(inout) :: x
    real(kind=8),dimension(:,:,:), intent(in)    :: b
    real(kind=8), intent(in):: tol
    real(kind=8), intent(out):: res   
    integer(kind=4), intent(in):: maxite
    integer(kind=4), intent(out)::nite

    ! local
    real*8:: rnorm,bnorm,res0,conv

    call norm(1,b,bnorm) ! norm of b on level=1

    call residual(1,x,b,grid(1)%b,rnorm) ! residual returns both 'r' and its norm
    res0 = rnorm/bnorm

    nite=0

    do while ((nite.lt.maxite).and.(res0.gt.tol))
       call mg_Fcycle(1)
       call add_to(1,x,grid(1)%x)
       call residual(1,x,b,grid(1)%b,rnorm)
       res = rnorm/bnorm
       conv=res0/res ! error reduction after this iteration
       res0=res
       nite=nite+1
       write(*,10) nite,res,conv
    enddo

10  format("ite = ",I," / res = ",G," / conv = ",G)

  end subroutine mg_solve

  !----------------------------------------
  subroutine mg_Vcycle(lev1)

    integer:: lev1,lev
    real*8:: rnorm

    do lev=lev1,nlevs-1
       call smooth(lev,grid(lev)%x,grid(lev)%b,npre)
       call residual(lev,grid(lev)%x,grid(lev)%b,grid(lev)%r,rnorm)
       call finetocoarse(lev,lev+1,grid(lev)%r,grid(lev+1)%b)
       call set_to_zero(lev+1,grid(lev+1)%x)
    enddo

    lev=nlevs
    call smooth(lev,grid(lev)%x,grid(lev)%b,ndeepest)

    do lev=nlevs-1,lev1,-1
       call coarsetofine(lev+1,lev,grid(lev+1)%x,grid(lev)%r)
       call add_to(lev,grid(lev)%x,grid(lev)%r) ! add x to r
       call smooth(lev,grid(lev)%x,grid(lev)%b,npost)
    enddo

  end subroutine mg_Vcycle

  !----------------------------------------
  subroutine mg_Fcycle(lev1)

    integer:: lev1,lev

    do lev=lev1,nlevs-1
       call finetocoarse(lev,lev+1,grid(lev)%b,grid(lev+1)%b)
    enddo

    lev=nlevs
    call set_to_zero(lev,grid(lev)%x)
    call smooth(lev,grid(lev)%x,grid(lev)%b,ndeepest)

    do lev=nlevs-1,lev1,-1
       call coarsetofine(lev+1,lev,grid(lev+1)%x,grid(lev)%x)
       call mg_Vcycle(lev)
    enddo

  end subroutine mg_Fcycle

end module mg_solvers
