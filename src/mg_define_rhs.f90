module mg_define_rhs

  use mg_grids
  use mg_mpi_exchange

  implicit none

   real(kind=rp):: dx, dy, dz

contains

  subroutine define_rhs(nxg, nyg, npxg)

    integer(kind=ip), intent(in) :: nxg,nyg
    integer(kind=ip), intent(in) :: npxg

    real(kind=rp), dimension(:,:,:), pointer :: rhs

    integer(kind=ip) :: nx, ny, nz
    integer(kind=ip) :: k, j, i
    integer(kind=ip) :: pi, pj

    real(kind=rp):: sizex, sizey, sizez
    real(kind=rp):: x0, y0, z0
    real(kind=rp):: x1, y1, z1

    real(kind=rp):: x,z

    real(kind=rp):: bet,sumglo

    integer(kind=ip) :: lev

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
       x =( real(i+(pi*nx),kind=rp)- 0.5_8) * dx  

       do j = 1,ny
          !y =real(j,kind=rp)

          do k = 1,nz

             z = -(real(nz-k+1,kind=rp)-0.5_8) * dz

             rhs(k,j,i) = &
                  exp(-bet * ((x-x0)**2 + (z-z0)**2)) - &
                  exp(-bet * ((x-x1)**2 + (z-z1)**2))

          enddo
       enddo
    enddo

    z = sum(rhs(:,1:ny,1:nx))
    call global_sum(1,z,sumglo)

    sumglo = sumglo / (nxg*nyg*nz)

    rhs(:,1:ny,1:nx) = rhs(:,1:ny,1:nx) - sumglo

    z = sum(rhs(:,1:ny,1:nx))
    call global_sum(1,z,sumglo)

    !write(*,*)'myrank - sum(rhs):', myrank, sum(rhs(:,1:ny,1:nx))
    if (myrank.eq.0)then
       write(*,*)'- sum(rhs):',sumglo
    endif

  end subroutine define_rhs

end module mg_define_rhs
