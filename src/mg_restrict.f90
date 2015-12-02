module mg_restrict
   !
   ! Collection of restriction subroutines
   !
   use mg_grids
   !      use mg_mpi
   implicit none



contains

   !------------------------------------------------------------
   subroutine restrict_xyz(l1,l2,x,y)
   !
   ! Restrict 'x' from fine level l1 to 'y' on coarse level l2=l1+1

   integer(kind=is) :: l1,l2
   real(kind=8),dimension(:,:,:),intent(in) :: x
   real(kind=8),dimension(:,:,:),intent(out) :: y
!!$   real(kind=rl),dimension( &
!!$         grid(l1)%nz,       &
!!$         1-grid(l1)%nh:grid(l1)%ny+grid(l1)%nh, &
!!$         1-grid(l1)%nh:grid(l1)%nx+grid(l1)%nh), intent(in) :: x
!!$   real(kind=rl),dimension( &
!!$         grid(l2)%nz,       &
!!$         1-grid(l2)%nh:grid(l2)%ny+grid(l2)%nh, &
!!$         1-grid(l2)%nh:grid(l2)%nx+grid(l2)%nh), intent(out) :: y


   ! local
   integer(kind=is) :: i,j,k,i2,j2,k2
   real(kind=8):: z

   ! 
   do i2=1,grid(l2)%nx
      i=2*i2-1
      do j2=1,grid(l2)%ny
         j=2*j2-1
         do k2=1,grid(l2)%nz
            k=2*k2-1
            z = x(k,j,i)  +x(k,j,i+1)  +x(k,j+1,i)  +x(k,j+1,i+1) &
              + x(k+1,j,i)+x(k+1,j,i+1)+x(k+1,j+1,i)+x(k+1,j+1,i+1)
            y(k2,j2,i2) = z * 0.125_8
          enddo
       enddo
    enddo

   end subroutine restrict_xyz

   !------------------------------------------------------------
   subroutine interp_xyz(l2,l1,x,y)
   !
   ! Transpose operation of restrict_xyz
   ! Interpolate 'x' from coarse level l2 to 'y' on fine level l1=l2-1

   integer(kind=is) :: l1,l2
   real(kind=8),dimension(:,:,:),intent(in) :: x
   real(kind=8),dimension(:,:,:),intent(out) :: y

   ! local
   integer(kind=is) :: i,j,k,i2,j2,k2

   ! 
   do i2=1,grid(l2)%nx
      i=2*i2-1
      do j2=1,grid(l2)%ny
         j=2*j2-1
         do k2=1,grid(l2)%nz
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

   end subroutine interp_xyz

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

    do j2=1,ny2
       j=2*j2-1       
       do i2=1,nx2
          i=2*i2-1
          y(1,j2,i2) = (x(1,j,i)+x(1,j,i+1)+x(1,j+1,i)+x(1,j+1,i+1))*0.25
       enddo
    enddo

  end subroutine restrict_xy

  !----------------------------------------
  subroutine restrict_zzz(l1,l2,x,y)

    integer:: l1,l2
    real*8,dimension(grid(l1)%nx,grid(l1)%ny,grid(l1)%nz) :: x
    real*8,dimension(grid(l2)%nx,grid(l2)%ny,grid(l2)%nz) :: y

    ! local
    integer:: i,j,k,k2
    integer:: nx,ny,nz

    nx = grid(l1)%nx
    ny = grid(l1)%ny
    nz = grid(l1)%nz      

    do k=1,nz
       k2=(k-1)/8+1
       if(mod(k,8).eq.1)then
          do j=1,ny
             do i=1,nx               
                y(i,j,k2) = x(i,j,k)*0.125
             enddo
          enddo
       else
          do j=1,ny
             do i=1,nx               
                y(i,j,k2) = y(i,j,k2)+x(i,j,k)*0.125
             enddo
          enddo
       endif
    enddo

  end subroutine restrict_zzz

  !----------------------------------------
  subroutine interpolate_zzz(l2,l1,y,x)

    integer(kind = 4), intent(in):: l1,l2
    real*8,dimension(:,:,:), intent(out):: x
    real*8,dimension(:,:,:), intent(in) :: y


    integer(kind = 4):: i, j, k, k2
    integer(kind = 4):: nx, ny, nz

    do k=1,nz
       k2=(k-1)/8+1
       do j=1,ny
          do i=1,nx               
             x(i,j,k)=y(i,j,k2)
          enddo
       enddo
    enddo

  end subroutine interpolate_zzz




end module mg_restrict
