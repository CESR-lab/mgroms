module mg_intergrids
  !
  ! Collection of restriction subroutines
  !
  use mg_tictoc
  use mg_grids

  implicit none

contains
  !------------!
  !- RESTRICT -!
  !------------------------------------------------------------
  subroutine restrict(lev)

    integer(kind=is), intent(in) :: lev

    real(kind=8),dimension(:,:,:),pointer :: r
    real(kind=8),dimension(:,:,:),pointer :: b

    integer(kind=is) :: nx, ny, nz

    nx = grid(lev+1)%nx
    ny = grid(lev+1)%ny
    nz = grid(lev+1)%nz

    r => grid(lev)%r
    b => grid(lev+1)%b

    grid(lev+1)%p = 0._8

    if ((aggressive).and.(lev==1)) then
       call restrict_aggressive(r,b,nx,ny,nz)

    elseif (grid(lev)%nz == 1) then
       call restrict_2D(r,b,nx,ny)

    else
       call tic(lev,'restrict_3D')

       call restrict_3D(r,b,nx,ny,nz)

      call toc(lev,'restrict_3D')

    end if

  end subroutine restrict

  !----------------------------------------
  subroutine restrict_aggressive(x,y,nx,ny,nz)

    real(kind=rl),dimension(:,:,:), intent(in) :: x !fine
    real(kind=rl),dimension(:,:,:), intent(inout) :: y ! coarse
    integer(kind=is), intent(in) :: nx, ny, nz

    ! local
    integer(kind=is):: i,j,k,k2

    do k=1,nz
       k2=(k-1)/8+1
       if(mod(k,8).eq.1)then
          do j=1,ny
             do i=1,nx               
                y(i,j,k2) = x(i,j,k)*0.125_8
             enddo
          enddo
       else
          do j=1,ny
             do i=1,nx               
                y(i,j,k2) = y(i,j,k2)+x(i,j,k)*0.125_8
             enddo
          enddo
       endif
    enddo

  end subroutine restrict_aggressive

  !------------------------------------------------------------
  subroutine restrict_2D(x,y,nx,ny)
    real(kind=rl),dimension(:,:,:),pointer,intent(in) :: x
    real(kind=rl),dimension(:,:,:),pointer,intent(out) :: y
    integer(kind=is), intent(in) :: nx, ny

    !TODO
    integer(kind=is) ::idum ! line to remove
    idum = nx               ! line to remove
    idum = ny               ! line to remove
    y = x                   ! line to remove
    write(*,*)'Error: prolong_2D  not available yet !'
    stop -1
    !TODO

  end subroutine restrict_2D

 !----------------------------------------
  subroutine restrict_xy(l1,l2,x,y)

    integer:: l1,l2
    real*8,dimension(grid(l1)%nz,grid(l1)%ny,grid(l1)%nx) :: x
    real*8,dimension(grid(l2)%nz,grid(l2)%ny,grid(l2)%nx) :: y

    ! local
    integer:: i,j,i2,j2
    integer:: nx2,ny2

    nx2 = grid(l2)%nx
    ny2 = grid(l2)%ny

!!$<<<<<<< HEAD
!!$    ! indices (nh+1,nh+2) on fine grid are glued to (nh+1) on coarse grid
!!$    do j2=2,ny2-1
!!$       j=2*(j2-nhalo)+1       ! take into account the halo!!!
!!$       do i2=2,nx2-1
!!$          i=2*(i2-nhalo)+1
!!$          y(i2,j2,1) = (x(i,j,1)+x(i+1,j,1)+x(i,j+1,1)+x(i+1,j+1,1))*0.25
!!$=======
    do j2=1,ny2
       j=2*j2-1       
       do i2=1,nx2
          i=2*i2-1
          y(1,j2,i2) = (x(1,j,i)+x(1,j,i+1)+x(1,j+1,i)+x(1,j+1,i+1))*0.25
!!$>>>>>>> 5d76062d572541a52c0574dba607c2e2d63cb883
       enddo
    enddo

  end subroutine restrict_xy

  !------------------------------------------------------------
  subroutine restrict_3D(x,y,nx,ny,nz)
    !
    ! Restrict 'x' from fine level l1 to 'y' on coarse level l2=l1+1
    real(kind=rl),dimension(:,:,:),pointer,intent(in) :: x
    real(kind=rl),dimension(:,:,:),pointer,intent(out) :: y
    integer(kind=is), intent(in) :: nx, ny, nz
    ! local
    integer(kind=is) :: i,j,k,i2,j2,k2
    real(kind=8):: z

    ! 
    do i2=1,nx
       i=2*i2-1
       do j2=1,ny
          j=2*j2-1
          do k2=1,nz
             k=2*k2-1
             z = x(k,j,i)  +x(k,j,i+1)  +x(k,j+1,i)  +x(k,j+1,i+1) &
                  + x(k+1,j,i)+x(k+1,j,i+1)+x(k+1,j+1,i)+x(k+1,j+1,i+1)
             y(k2,j2,i2) = z * 0.125_8
          enddo
       enddo
    enddo

  end subroutine restrict_3D

  !-----------!
  !- PROLONG -!
  !------------------------------------------------------------
  subroutine prolong(lev)

    !- prolong from level lev+1 to level lev
    integer(kind=is), intent(in) :: lev

    real(kind=8),dimension(:,:,:),pointer :: pf
    real(kind=8),dimension(:,:,:),pointer :: pc

    integer(kind=is) :: nxc, nyc, nzc

    nxc = grid(lev+1)%nx
    nyc = grid(lev+1)%ny
    nzc = grid(lev+1)%nz

    pf => grid(lev)%p
    pc => grid(lev+1)%p

    if ((aggressive).and.(lev==1)) then
       call prolong_aggressive(pf,pc,nxc,nyc,nzc)

    elseif (grid(lev)%nz == 1) then
       call prolong_2D(pf,pc,nxc,nyc)

    else
       call prolong_3D(pf,pc,nxc,nyc,nzc)

    end if

  end subroutine prolong

  !------------------------------------------------------------
  subroutine prolong_aggressive(x,y,nx,ny,nz)
    real(kind=8),dimension(:,:,:),intent(in)  :: x
    real(kind=8),dimension(:,:,:),intent(out) :: y
    integer(kind=is),intent(in) :: nx, ny, nz

    !TODO
    integer(kind=is) ::idum ! line to remove
    idum = nx               ! line to remove
    idum = ny               ! line to remove
    idum = nz               ! line to remove
    y = x                   ! line to remove
    write(*,*)'Error:  prolong_aggressive not available yet !'
    stop -1
    !TODO

  end subroutine prolong_aggressive

  !------------------------------------------------------------
  subroutine prolong_2D(x,y,nx,ny)
    real(kind=8),dimension(:,:,:),intent(in)  :: x
    real(kind=8),dimension(:,:,:),intent(out) :: y
    integer(kind=is),intent(in) :: nx, ny

    !TODO
    integer(kind=is) ::idum ! line to remove
    idum = nx               ! line to remove
    idum = ny               ! line to remove
    y = x                   ! line to remove
    write(*,*)'Error: prolong_2D  not available yet !'
    stop -1
    !TODO

  end subroutine prolong_2D

  !------------------------------------------------------------
  subroutine prolong_3D(x,y,nx,ny,nz)
    real(kind=8),dimension(:,:,:),intent(in)  :: x
    real(kind=8),dimension(:,:,:),intent(out) :: y
    integer(kind=is),intent(in) :: nx, ny, nz

    ! local
    integer(kind=is) :: i,j,k,i2,j2,k2
    ! 
    do i2=1,nx
       i=2*i2-1
       do j2=1,ny
          j=2*j2-1
          do k2=1,nz
             k=2*k2-1
             y(k  ,j  ,i  ) = x(k2,j2,i2)
             y(k+1,j  ,i  ) = x(k2,j2,i2)
             y(k  ,j+1,i  ) = x(k2,j2,i2)
             y(k+1,j+1,i  ) = x(k2,j2,i2)
             y(k  ,j  ,i+1) = x(k2,j2,i2)
             y(k+1,j  ,i+1) = x(k2,j2,i2)
             y(k  ,j+1,i+1) = x(k2,j2,i2)
             y(k+1,j+1,i+1) = x(k2,j2,i2)
          enddo
       enddo
    enddo

  end subroutine prolong_3D

!!$  !----------------------------------------
!!$  subroutine interpolate_zzz(l2,l1,y,x)
!!$
!!$    integer(kind = 4), intent(in):: l1,l2
!!$    real*8,dimension(:,:,:), intent(out):: x
!!$    real*8,dimension(:,:,:), intent(in) :: y
!!$
!!$    integer(kind = 4):: i, j, k, k2
!!$    integer(kind = 4):: nx, ny, nz
!!$
!!$    do k=1,nz
!!$       k2=(k-1)/8+1
!!$       do j=1,ny
!!$          do i=1,nx               
!!$             x(i,j,k)=y(i,j,k2)
!!$          enddo
!!$       enddo
!!$    enddo
!!$
!!$  end subroutine interpolate_zzz

end module mg_intergrids
