module mg_define_rhs

  use mg_grids

  implicit none

   real(kind=rl):: dx, dy, dz

contains

  subroutine define_rhs(nxg, nyg, nzg, npxg, npyg)

    integer(kind=is), intent(in) :: nxg, nyg, nzg
    integer(kind=is), intent(in) :: npxg, npyg

    real(kind=rl), dimension(:,:,:), pointer :: rhs

    integer(kind=is) :: nx, ny, nz
    integer(kind=is) :: k, j, i
    integer(kind=is) :: pi, pj

    real(kind=rl):: sizex, sizey, sizez
    real(kind=rl):: x0, y0, z0
    real(kind=rl):: x1, y1, z1

    real(kind=rl):: x,y,z

    real(kind=rl):: bet

    integer(kind=is) :: lev

    sizex=1.d4
    sizey=1.d4
    sizez=4.0d3

    bet = 600._8 / (sizex*sizex)

    x0 = sizex * 0.65_8
    y0 = sizey * 0._8
    z0 = sizez * (0.75_8 - 1._8)

    x1 = sizex * 0.75_8
    y1 = sizey * 0._8
    z1 = sizez * (0.65_8 - 1._8)

    lev = 1

    nx = grid(lev)%nx
    ny = grid(lev)%ny
    nz = grid(lev)%nz

    rhs => grid(lev)%b

    dx = sizex / nxg
    dy = sizey / nyg
    dz = sizez / nz

    pj = myrank/npxg   
    pi = myrank-pj*npxg 

    do i = 1,nx
       ! nx is not necessarily the same for all processes
       x =( real(i+(pi*nx),kind=rl)- 0.5_8) * dx  

       do j = 1,ny
          !y =real(j,kind=rl)

          do k = 1,nz
             z = -(real(k,kind=rl)- 0.5_8) * dz

             rhs(k,j,i) = &
                  exp(-bet * ((x-x0)**2 + (z-z0)**2)) - &
                  exp(-bet * ((x-x1)**2 + (z-z1)**2))

          enddo
       enddo
    enddo

    write(*,*)'myrank - sum(rhs):', myrank, sum(rhs(:,1:ny,1:nx))

  end subroutine define_rhs

end module mg_define_rhs
