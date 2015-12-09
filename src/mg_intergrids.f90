module mg_intergrids

  use mg_mpi
  use mg_tictoc
  use mg_namelist
  use mg_grids
  use mg_mpi_exchange
  use mg_gather

  implicit none

contains
  !---------------!
  !- FINE2COARSE -! fine to coarse grid
  !------------------------------------------------------------
  subroutine fine2coarse(lev)

    ! coarsen grid(lev)%r to grid(lev+1)%b

    integer(kind=ip), intent(in) :: lev

    real(kind=rp),dimension(:,:,:),pointer :: r
    real(kind=rp),dimension(:,:,:),pointer :: b

    integer(kind=ip) :: nx, ny, nz

    nx = grid(lev+1)%nx
    ny = grid(lev+1)%ny
    nz = grid(lev+1)%nz

    r => grid(lev)%r

    if (grid(lev+1)%gather == 1) then
       b => grid(lev+1)%dummy3
       nx = grid(lev+1)%nx / grid(lev+1)%ngx
       ny = grid(lev+1)%ny / grid(lev+1)%ngy
!       if(myrank == 0) write(*,*)"gather lev=",lev+1,"nx,ny,nz=",nx,ny,nz
    else
       b => grid(lev+1)%b
!       if(myrank == 0) write(*,*)"F2C   lev=",lev+1,"nx,ny,nz=",nx,ny,nz
    endif

    grid(lev+1)%p = 0._8

    if ((aggressive).and.(lev==1)) then
       call fine2coarse_aggressive(r,b,nx,ny,nz)

    elseif (grid(lev)%nz == 1) then
       call fine2coarse_2D(r,b,nx,ny)

    else
       call tic(lev,'fine2coarse_3D')

       call fine2coarse_3D(r,b,nx,ny,nz)

      call toc(lev,'fine2coarse_3D')

    end if

    if (grid(lev+1)%gather == 1) then
!       if(myrank == 0) write(*,*)" *** dummy3(1,1,1)=",b(1,1:ny,1:nx)
       r => grid(lev+1)%dummy3
       b => grid(lev+1)%b
       call gather(lev+1,r,b)
!       if(myrank == 0) write(*,*)" *** after gather **** b(1,1,1)=",grid(lev+1)%b(1,1,1)
    endif

    call fill_halo(lev+1,grid(lev+1)%b)

  end subroutine fine2coarse
  !----------------------------------------
  subroutine fine2coarse_aggressive(x,y,nx,ny,nz)

    real(kind=rp)   , dimension(:,:,:), intent(in)    :: x !fine
    real(kind=rp)   , dimension(:,:,:), intent(inout) :: y ! coarse
    integer(kind=ip)                  , intent(in)    :: nx, ny, nz

    ! local
    integer(kind=ip):: i,j,k,k2

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

  end subroutine fine2coarse_aggressive

  !------------------------------------------------------------
  subroutine fine2coarse_2D(x,y,nx,ny)
    real(kind=rp),dimension(:,:,:),pointer,intent(in) :: x
    real(kind=rp),dimension(:,:,:),pointer,intent(out) :: y
    integer(kind=ip), intent(in) :: nx, ny

    !TODO
    integer(kind=ip) ::idum ! line to remove
    idum = nx               ! line to remove
    idum = ny               ! line to remove
    y = x                   ! line to remove
    write(*,*)'Error: coarse2fine_2D  not available yet !'
    stop -1
    !TODO

  end subroutine fine2coarse_2D

 !----------------------------------------
  subroutine fine2coarse_xy(l1,l2,x,y)

    integer:: l1,l2
    real*8,dimension(grid(l1)%nz,grid(l1)%ny,grid(l1)%nx) :: x
    real*8,dimension(grid(l2)%nz,grid(l2)%ny,grid(l2)%nx) :: y

    ! local
    integer:: i,j,i2,j2
    integer:: nx2,ny2

    nx2 = grid(l2)%nx
    ny2 = grid(l2)%ny

    do j2=1,ny2
       j=2*j2-1       
       do i2=1,nx2
          i=2*i2-1
          y(1,j2,i2) = (x(1,j,i)+x(1,j,i+1)+x(1,j+1,i)+x(1,j+1,i+1))*0.25
       enddo
    enddo

  end subroutine fine2coarse_xy

  !------------------------------------------------------------
  subroutine fine2coarse_3D(x,y,nx,ny,nz)
    !
    ! Fine2coarse 'x' from fine level l1 to 'y' on coarse level l2=l1+1
    real(kind=rp),dimension(:,:,:),pointer,intent(in) :: x
    real(kind=rp),dimension(:,:,:),pointer,intent(out) :: y
    integer(kind=ip), intent(in) :: nx, ny, nz
    ! local
    integer(kind=ip) :: i,j,k,i2,j2,k2
    real(kind=rp):: z

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

  end subroutine fine2coarse_3D

  !---------------!
  !- COARSE2FINE -! coarse to fine grid
  !------------------------------------------------------------
  subroutine coarse2fine(lev)

    !- coarse2fine from level lev+1 to level lev
    integer(kind=ip), intent(in) :: lev

    real(kind=rp),dimension(:,:,:),pointer :: rf
    real(kind=rp),dimension(:,:,:),pointer :: pc

    integer(kind=ip) :: nxc, nyc, nzc

    nxc = grid(lev+1)%nx
    nyc = grid(lev+1)%ny
    nzc = grid(lev+1)%nz


    pc => grid(lev+1)%p

    if (grid(lev+1)%gather == 1) then
!       if(myrank==0)write(*,*)"SPLIT!!!"
       rf => grid(lev+1)%dummy3
       call split(lev+1,pc,rf)
       pc => grid(lev+1)%dummy3
       nxc = grid(lev+1)%nx / grid(lev+1)%ngx
       nyc = grid(lev+1)%ny / grid(lev+1)%ngy
       nzc = grid(lev+1)%nz       
    endif
    rf => grid(lev)%r


   if ((aggressive).and.(lev==1)) then
       call coarse2fine_aggressive(rf,pc,nxc,nyc,nzc)

    elseif (grid(lev)%nz == 1) then
       call coarse2fine_2D(rf,pc,nxc,nyc)

    else
       call coarse2fine_3D(rf,pc,nxc,nyc,nzc)
    end if

    call fill_halo(lev,grid(lev)%r)

    grid(lev)%p = grid(lev)%p + grid(lev)%r


  end subroutine coarse2fine

  !------------------------------------------------------------
  subroutine coarse2fine_aggressive(x,y,nx,ny,nz)
    real(kind=rp),dimension(:,:,:),pointer,intent(in)  :: x
    real(kind=rp),dimension(:,:,:),pointer,intent(out) :: y
    integer(kind=ip),intent(in) :: nx, ny, nz

    !TODO
    integer(kind=ip) ::idum ! line to remove
    idum = nx               ! line to remove
    idum = ny               ! line to remove
    idum = nz               ! line to remove
    y = x                   ! line to remove
    write(*,*)'Error:  coarse2fine_aggressive not available yet !'
    stop -1
    !TODO

  end subroutine coarse2fine_aggressive

  !------------------------------------------------------------
  subroutine coarse2fine_2D(x,y,nx,ny)
    real(kind=rp),dimension(:,:,:),pointer,intent(in)  :: x
    real(kind=rp),dimension(:,:,:),pointer,intent(out) :: y
    integer(kind=ip),intent(in) :: nx, ny

    !TODO
    integer(kind=ip) ::idum ! line to remove
    idum = nx               ! line to remove
    idum = ny               ! line to remove
    y = x                   ! line to remove
    write(*,*)'Error: coarse2fine_2D  not available yet !'
    stop -1
    !TODO

  end subroutine coarse2fine_2D

  !------------------------------------------------------------
  subroutine coarse2fine_3D(xf,xc,nx,ny,nz)
    real(kind=rp),dimension(:,:,:),pointer,intent(out) :: xf
    real(kind=rp),dimension(:,:,:),pointer,intent(in)  :: xc
    integer(kind=ip),intent(in) :: nx, ny, nz

    ! local
    integer(kind=ip) :: i,j,k,i2,j2,k2
    ! 
    do i2=1,nx
       i=2*i2-1
       do j2=1,ny
          j=2*j2-1
          do k2=1,nz
             k=2*k2-1
             xf(k  ,j  ,i  ) = xc(k2,j2,i2)
             xf(k+1,j  ,i  ) = xc(k2,j2,i2)
             xf(k  ,j+1,i  ) = xc(k2,j2,i2)
             xf(k+1,j+1,i  ) = xc(k2,j2,i2)
             xf(k  ,j  ,i+1) = xc(k2,j2,i2)
             xf(k+1,j  ,i+1) = xc(k2,j2,i2)
             xf(k  ,j+1,i+1) = xc(k2,j2,i2)
             xf(k+1,j+1,i+1) = xc(k2,j2,i2)
          enddo
       enddo
    enddo

  end subroutine coarse2fine_3D

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
