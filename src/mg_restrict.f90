module mg_restrict
   implicit none
   !
   ! Collection of restriction subroutines
   !
   use mg_grids
   !      use mg_mpi



contains

   !------------------------------------------------------------
   subroutine restrict_xyz(l1,l2,x,y)
   !
   ! Restrict the residual from fine to coarse grid

   integer:: l1,l2
   real*8,dimension(grid(l1)%nx,grid(l1)%ny,grid(l1)%nz) :: x
   real*8,dimension(grid(l2)%nx,grid(l2)%ny,grid(l2)%nz) :: y

   ! local
   integer:: i,j,k,i2,j2,k2
   integer:: nx2,ny2,nz2
   real*8:: z

   nx2 = grid(l2)%nx
   ny2 = grid(l2)%ny
   nz2 = grid(l2)%nz      

   do k2=1,nz2
       k=2*k2-1
       ! indices (nh+1,nh+2) on fine grid are glued to (nh+1) on coarse grid
       do j2=2,ny2-1
          j=2*(j2-nhalo)+1 ! take into account the halo!!!
          do i2=2,nx2-1
             i=2*(i2-nhalo)+1
             z = x(i,j,k)+x(i+1,j,k)+x(i,j+1,k)+x(i+1,j+1,k) &
                        + x(i,j,k+1)+x(i+1,j,k+1)+x(i,j+1,k+1)+x(i+1,j+1,k+1)
             y(i2,j2,k2) = z * 0.125
          enddo
       enddo
    enddo

   end subroutine restrict_xyz

  !----------------------------------------
  subroutine restrict_xy(l1,l2,x,y)

    integer:: l1,l2
    real*8,dimension(grid(l1)%nx,grid(l1)%ny,grid(l1)%nz) :: x
    real*8,dimension(grid(l2)%nx,grid(l2)%ny,grid(l2)%nz) :: y

    ! local
    integer:: i,j,k,i2,j2,k2
    integer:: nx2,ny2,nz2

    nx2 = grid(l2)%nx
    ny2 = grid(l2)%ny
    nz2 = grid(l2)%nz      


    ! indices (nh+1,nh+2) on fine grid are glued to (nh+1) on coarse grid
    do j2=2,ny2-1
       j=2*(j2-nhalo)+1       ! take into account the halo!!!
       do i2=2,nx2-1
          i=2*(i2-nhalo)+1
          y(i2,j2,1) = (x(i,j,1)+x(i+1,j,1)+x(i,j+1,1)+x(i+1,j+1,1))*0.25
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

!!NG: 16 nov 2015 comment this -> #if defined FULLSET
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

  !----------------------------------------
  subroutine smooth_sor_5(lev,x,b,nite)

    integer(kind = 4), intent(in):: lev
    integer(kind = 4), intent(in):: nite
    real*8,dimension(:,:,:), intent(out):: x
    real*8,dimension(:,:,:), intent(in) :: b

    integer(kind = 4):: k, kt, j, i
    integer(kind = 4):: nx, ny
    real(kind=8)     :: c1, c2, c3

    c1 = omega ! relaxation parameter
    c2=1-c1
    k=1
    do kt=1,nite
       ! SOR is hard to multithread...
       do j=2,ny-1
          do i=2,nx-1
!!$             c3=c1/abs(A(2,i,j,k))
!!$             x(i,j,k) = x(i,j,k)*c2 + c3*(                    &
!!$             &                        A(0,i,j,k)  *x(i,j-1,k) &
!!$             &                        A(1,i,j,k)  *x(i-1,j,k) &
!!$             &                        A(1,i+1,j,k)*x(i+1,j,k) &
!!$             &                        A(1,i,j+1,k)*x(i,j+1,k) &
!!$             &                       -b(i,j,k) )
          enddo
       enddo
       ! improve that, frequency of call depends on the halo width
!!$       call fill_halo(lev,x)
    enddo

  end subroutine smooth_sor_5

  !----------------------------------------
  subroutine residual_15(lev,x,b,r)
    ! WATCH OUT: THERE IS NO GHOSTPOINT IN Z 
    !
    ! remember the Gmeiner 2015 paper, if a line has too many
    ! instructions there is a chance to saturate the registers
    !

    integer(kind = 4), intent(in) :: lev
    real(kind = 8), dimension(:,:,:), intent(in)  :: x
    real(kind = 8), dimension(:,:,:), intent(in)  :: b
    real(kind = 8), dimension(:,:,:), intent(out) :: r

    integer(kind = 4):: i, j, k, km, kp, nx, ny, nz

!!$    ! bottom level
!!$    do k=1,1
!!$       kp=k+1
!!$       do j=2,ny-1
!!$          do i=2,nx-1
!!$             r(i,j,k) = b(i,j,k)                    &
!!$             &              - A(5,i,j,k)*x(i,j-1,k) &
!!$             &              - A(6,i,j,k)*x(i-1,j,k) &
!!$             &              - A(7,i,j,k)*x(i,j,k)   &
!!$             r(i,j,k) = r(i,j,k)                         &
!!$             &              - A(6,i+1,j,k)*x(i+1,j,k)    &
!!$             &              - A(5,i,j+1,k)*x(i,j+1,k)    &
!!$             &              - A(4,i,j-1,k+1)*x(i,j-1,kp) &
!!$             &              - A(3,i-1,j,k+1)*x(i-1,j,kp) &
!!$             &              - A(2,i,j,k+1)*x(i,j,kp)     &
!!$             &              - A(1,i+1,j,k+1)*x(i+1,j,kp) & 
!!$             &              - A(0,i,j+1,k+1)*x(i,j+1,kp)
!!$          enddo
!!$       enddo
!!$    enddo
!!$
!!$    ! interior levels
!!$    do k=2,nz-1
!!$       km=k-1
!!$       kp=k+1
!!$       do j=2,ny-1
!!$          do i=2,nx-1
!!$             r(i,j,k) = b(i,j,k)                     &
!!$             &              - A(0,i,j,k)*x(i,j-1,km) &
!!$             &              - A(1,i,j,k)*x(i-1,j,km) &
!!$             &              - A(2,i,j,k)*x(i,j,km)   &
!!$             &              - A(3,i,j,k)*x(i+1,j,km) &
!!$             &              - A(4,i,j,k)*x(i,j+1,km) &
!!$             &              - A(5,i,j,k)*x(i,j-1,k)  &
!!$             &              - A(6,i,j,k)*x(i-1,j,k)  &
!!$             &              - A(7,i,j,k)*x(i,j,k)
!!$             r(i,j,k) = r(i,j,k)                     &
!!$             &              - A(6,i,j,k)*x(i+1,j,k)  &
!!$             &              - A(5,i,j,k)*x(i,j+1,k)  &
!!$             &              - A(4,i,j,k)*x(i,j-1,kp) &
!!$             &              - A(3,i,j,k)*x(i-1,j,kp) &
!!$             &              - A(2,i,j,k)*x(i,j,kp)   &
!!$             &              - A(1,i,j,k)*x(i+1,j,kp) &
!!$             &              - A(0,i,j,k)*x(i,j+1,kp) 
!!$          enddo
!!$       enddo
!!$    enddo
!!$
!!$    ! top level
!!$    do k=nz,nz
!!$       km=k-1
!!$       do j=2,ny-1
!!$          do i=2,nx-1
!!$             r(i,j,k) = b(i,j,k)                     &
!!$             &              - A(0,i,j,k)*x(i,j-1,km) &
!!$             &              - A(1,i,j,k)*x(i-1,j,km) &
!!$             &              - A(2,i,j,k)*x(i,j,km)   &
!!$             &              - A(3,i,j,k)*x(i+1,j,km) &
!!$             &              - A(4,i,j,k)*x(i,j+1,km) &
!!$             &              - A(5,i,j,k)*x(i,j-1,k)  &
!!$             &              - A(6,i,j,k)*x(i-1,j,k)  &
!!$             &              - A(7,i,j,k)*x(i,j,k)    &
!!$             r(i,j,k) = r(i,j,k)                     &
!!$             &              - A(6,i,j,k)*x(i+1,j,k)  &
!!$             &              - A(5,i,j,k)*x(i,j+1,k)  
!!$          enddo
!!$       enddo
!!$    enddo

!!$    call fill_halo(lev,r)

  end subroutine residual_15

  !----------------------------------------
  subroutine restrict_matrix_xyz(lev) ! define A(lev+1) from A(lev)

    ! do the matrix multiplication A(lev+1)=R*A(lev)*I

    ! WATCH OUT: THERE IS NO GHOSTPOINT IN Z

    integer(kind = 4), intent(in) :: lev


    ! local
    integer(kind=4):: i,j,k,i2,j2,k2,m
    integer(kind=4):: l
    integer(kind=4):: kx, ky, kz
    integer(kind=4):: nx2,ny2,nz2
    real(kind=8) :: a
    integer(kind=4) :: jx, jy, jz

!!$    do k2=1,nz2
!!$       do j2=2,ny2-1
!!$          do i2=2,nx2-1
!!$             do l=0,7     
!!$                a = 0.
!!$                do k=0,7
!!$                   kx=mod(k,2)
!!$                   ky=mod(k/2,2)
!!$                   kz=mod(k/4,2)
!!$                   do m=0,7
!!$                      jx=ix15(m)
!!$                      jy=iy15(m)
!!$                      jz=iz15(m)                        
!!$                      a = a+cff(kx)*cff(ky)*cff(kz)                    &
!!$                      &                       *cff(jx)*cff(jy)*cff(jz) &
!!$                      &                       *grid(lev)%A(i,j,k,m)
!!$                   enddo
!!$                enddo
!!$                grid(lev+1)%A(i2,j2,k2,l) = a
!!$             enddo
!!$          enddo
!!$       enddo
!!$    enddo

!!$    call fill_halo_matrix(lev+1,grid(lev+1)%A)

  end subroutine restrict_matrix_xyz
!!NG: 16 nov 2015 comment this -> #endif

end module mg_optimized
